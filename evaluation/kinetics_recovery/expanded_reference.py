#!/usr/bin/env python3
"""Build isolated expanded-reference summaries; never change canonical scores.

Without --artifact-root, regenerate compact summaries from committed scores.
With --artifact-root, also validate frozen sources, rebuild the gene scope,
and recompute TSvelo organoid alpha using the shared kinetics processing.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

MODULE = Path(__file__).resolve().parent
RESULTS = MODULE / "results/expanded_reference"
CANONICAL = MODULE / "results/central_results"
KEYS = ["Method", "Dataset", "Branch", "Parameter", "Gene"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def candidate_genes(scope: pd.DataFrame, dataset: str) -> list[str]:
    rows = scope[scope.Dataset.eq(dataset) & scope.Included]
    return rows.sort_values("Candidate_order").Gene.tolist()


def load_scores(results: Path = RESULTS) -> tuple[pd.DataFrame, pd.DataFrame]:
    manifest = json.loads((results / "manifest.json").read_text())
    if sha256(CANONICAL / "gene_metrics.csv") != manifest["canonical_gene_metrics_sha256"]:
        raise ValueError("Canonical scores changed: re-audit the expanded reference scope")
    scope = pd.read_csv(results / "gene_scope.csv")
    canonical = pd.read_csv(CANONICAL / "gene_metrics.csv")
    extra = pd.read_csv(results / "tsvelo_organoid_gene_metrics.csv")
    if not extra.Method.eq("TSvelo").all() or not extra.Dataset.eq("organoid").all():
        raise ValueError("The isolated extra table must contain only organoid TSvelo")
    if not extra.Parameter.eq("alpha").all():
        raise ValueError("TSvelo expanded recovery supports alpha only")
    if extra["Branch"].isna().any() or not extra["Branch"].isin([1, 2]).all():
        raise ValueError("TSvelo expanded scores require branches 1 and 2")
    frame = pd.concat([canonical, extra], ignore_index=True)
    if frame.duplicated(KEYS).any():
        raise ValueError("Expanded scores overlap canonical scoring keys")
    selected = pd.concat([
        frame[frame.Dataset.eq(dataset) & frame.Gene.isin(candidate_genes(scope, dataset))]
        for dataset in ("rpe1", "organoid")
    ], ignore_index=True)
    return scope, selected


def summaries(scope: pd.DataFrame, scores: pd.DataFrame) -> dict[str, pd.DataFrame]:
    coverage, common, membership = [], [], []
    for (dataset, branch, parameter), frame in scores.groupby(
        ["Dataset", "Branch", "Parameter"], dropna=False, sort=True
    ):
        genes = candidate_genes(scope, dataset)
        panel = f"{dataset}_{'all' if pd.isna(branch) else 'branch'+str(int(branch))}_{parameter}"
        valid_sets = {
            method: set(group.loc[np.isfinite(group.Spearman_Corr), "Gene"])
            for method, group in frame.groupby("Method")
        }
        conventional = {m: g for m, g in valid_sets.items() if not m.startswith("Dynamo")}
        groups = {
            "core_conventional": {m: g for m, g in conventional.items() if m != "TSvelo"},
            "all_conventional": conventional,
            "all_displayed_methods": valid_sets,
        }
        intersections = {name: set.intersection(*group.values()) for name, group in groups.items()}
        for name, intersection in intersections.items():
            membership.extend({"Panel": panel, "Comparison_scope": name, "Gene": gene}
                              for gene in genes if gene in intersection)
        for method, group in frame.groupby("Method"):
            values = group.loc[np.isfinite(group.Spearman_Corr), "Spearman_Corr"]
            coverage.append({
                "Panel": panel, "Dataset": dataset, "Branch": branch,
                "Parameter": parameter, "Method": method, "Candidate_N": len(genes),
                "Present_N": len(group), "Valid_N": len(values),
                "Undefined_N": len(group)-len(values), "Unavailable_N": len(genes)-len(group),
                "Mean": values.mean(), "Median": values.median(), "SD_across_genes": values.std(),
                "Role": "label_data_control" if method.startswith("Dynamo") else "conventional_us",
                "Provenance": "expanded_tsvelo_organoid" if dataset == "organoid" and method == "TSvelo"
                else "canonical_central_results",
            })
            for name, intersection in intersections.items():
                if method not in groups[name]:
                    continue
                chosen = group.loc[group.Gene.isin(intersection), "Spearman_Corr"]
                common.append({
                    "Panel": panel, "Comparison_scope": name, "Common_genes_N": len(intersection),
                    "Method": method, "Method_valid_N": len(chosen),
                    "Mean": chosen.mean(), "Median": chosen.median(),
                })
    return {"coverage_summary": pd.DataFrame(coverage), "common_gene_summary": pd.DataFrame(common),
            "common_gene_membership": pd.DataFrame(membership)}


def frozen_inputs(artifact_root: Path) -> dict[str, Path]:
    from .metrics import source_path
    inputs = json.loads((RESULTS / "manifest.json").read_text())["external_inputs"]
    resolved = {}
    for key, source in inputs.items():
        path = source_path(artifact_root, source["path"])
        if sha256(path) != source["sha256"]:
            raise ValueError(f"Frozen input hash mismatch: {key}")
        resolved[key] = path
    return resolved


def build_scope(inputs: dict[str, Path]) -> pd.DataFrame:
    """Intersect published reference genes with frozen conventional score genes."""
    scores = pd.read_csv(CANONICAL / "gene_metrics.csv")
    rpe = pd.read_csv(inputs["table_s1"], skiprows=1)
    org = pd.read_csv(inputs["table_s4"], skiprows=1)
    table_s3 = pd.read_csv(inputs["table_s3"], skiprows=1)
    selected = set(table_s3.loc[table_s3.selected, "geneName"])
    rpe_genes = list(pd.read_csv(CANONICAL / "matrices/rpe1/reference_alpha.csv", index_col=0).index)
    org_scores = scores[scores.Dataset.eq("organoid") & scores.Method.isin(
        ["cellDancer", "DeepVelo", "DeepKINET"])]
    sets = [set(g.Gene) for _, g in org_scores.groupby(["Method", "Branch", "Parameter"])]
    if not all(g == sets[0] for g in sets) or len(sets[0]) != 270 or len(selected) != 295:
        raise ValueError("Conventional organoid scope no longer matches the frozen 270/295 contract")
    org_genes = [g for g in org.gene_symbol if g in sets[0]]
    if len(rpe_genes) != 250 or not set(org_genes) <= selected:
        raise ValueError("Unexpected expanded candidate universe")
    rows = []
    for dataset, reference, candidates in (("rpe1", rpe, rpe_genes), ("organoid", org, org_genes)):
        order = {gene: n+1 for n, gene in enumerate(candidates)}
        for row in reference.to_dict("records"):
            gene = row["gene_symbol"]
            rows.append({
                "Dataset": dataset, "Gene": gene, "Reference_order": int(row["position"]),
                "Candidate_order": order.get(gene), "Included": gene in order,
                "Source_selected": gene in selected if dataset == "organoid" else True,
                "Strategy_group": row["strategy_group"],
                "Exclusion_reason": "" if gene in order else
                "housekeeping_control" if dataset == "organoid" and gene not in selected else
                "not_in_common_conventional_input",
            })
    return pd.DataFrame(rows)


def score_tsvelo(inputs: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Validate mouse gene IDs and merge overlapping lineage alpha by cell ID."""
    import anndata as ad
    from .metrics import (calculate_metrics, csv_roundtrip, layer_frame,
                          organoid_raw_from_layer, process_organoid,
                          read_organoid_reference, select_organoid_branch)
    settings = json.loads((MODULE / "config.json").read_text())["datasets"]["organoid"]
    reference = pd.read_csv(inputs["table_s4"], skiprows=1).set_index("gene_symbol")
    folded = {gene.casefold(): gene for gene in reference.index}
    if len(folded) != len(reference):
        raise ValueError("Reference symbols have case-insensitive collisions")
    source = ad.read_h5ad(inputs["organoid_source"], backed="r")
    ids, cells = source.var["ID"].copy(), source.obs_names.copy()
    source.file.close()
    frames, metadata, mappings = [], [], []
    obs_keys = [settings[k] for k in ("branch_key", "cluster_key", "pseudotime_key")]
    for key in ("tsvelo_lineage0", "tsvelo_lineage1"):
        result = ad.read_h5ad(inputs[key])
        if result.obs_names.has_duplicates or result.var_names.has_duplicates:
            raise ValueError("TSvelo cell and gene indices must be unique")
        if not result.obs_names.isin(cells).all():
            raise ValueError("TSvelo cells are absent from the frozen organoid input")
        names = [g for g in result.var_names if g.casefold() in folded]
        rename = {g: folded[g.casefold()] for g in names}
        if len(set(rename.values())) != len(names):
            raise ValueError("TSvelo symbols have case-insensitive collisions")
        mapping = pd.DataFrame([
            {"model_symbol": g, "reference_symbol": rename[g], "gene_id": result.var.loc[g, "ID"],
             "original_id_match": result.var.loc[g, "ID"] == ids.loc[rename[g]]}
            for g in names
        ])
        if not mapping.original_id_match.all():
            raise ValueError("Case normalization failed original mouse gene-ID validation")
        mappings.append(mapping)
        alpha = layer_frame(result[:, names], "alpha").rename(columns=rename)
        if not np.isfinite(alpha.to_numpy()).all():
            raise ValueError("TSvelo alpha contains non-finite values")
        frames.append(alpha)
        metadata.append(result.obs[obs_keys].copy())
    pd.testing.assert_frame_equal(mappings[0], mappings[1])
    meta = pd.concat(metadata)
    duplicates = meta[meta.index.duplicated(keep=False)].groupby(level=0).nunique(dropna=False)
    if (duplicates > 1).any().any():
        raise ValueError("Overlapping TSvelo cells have inconsistent branch metadata")
    meta = meta[~meta.index.duplicated(keep="first")]
    alpha = pd.concat(frames).groupby(level=0).mean()
    meta = meta.loc[alpha.index].copy()
    # Exported categorical values may be integer-valued floats.
    for key in obs_keys[:2]:
        meta[key] = meta[key].astype(str).str.replace(r"\.0$", "", regex=True)
    if len(alpha) != settings["n_obs"]:
        raise ValueError("Merged TSvelo cells do not match the frozen organoid input")
    merged = ad.AnnData(obs=meta, var=pd.DataFrame(index=alpha.columns))
    merged.layers["alpha"] = alpha.to_numpy()
    rows = []
    for branch, branch_settings in settings["branches"].items():
        selected, aligned_time = select_organoid_branch(merged, settings, branch)
        raw = csv_roundtrip(organoid_raw_from_layer(selected, aligned_time, "alpha", settings))
        prediction = csv_roundtrip(process_organoid(raw, settings))
        truth = csv_roundtrip(read_organoid_reference(
            inputs["table_s4"], branch, "alpha", branch_settings["reference_reverse_time"]
        ).loc[prediction.index])
        scores = calculate_metrics(prediction, truth)
        reasons = []
        for gene in scores.Gene:
            x, y = prediction.loc[gene].to_numpy(float), truth.loc[gene].to_numpy(float)
            valid = np.isfinite(x) & np.isfinite(y)
            reasons.append("fewer_than_two_pairs" if valid.sum() < 2 else
                           "constant_prediction" if np.ptp(x[valid]) == 0 else
                           "constant_reference" if np.ptp(y[valid]) == 0 else "")
        scores["Undefined_reason"] = reasons
        scores["Status"] = np.where(scores.Spearman_Corr.notna(), "valid", "undefined")
        scores["Method"], scores["Dataset"] = "TSvelo", "organoid"
        scores["Branch"], scores["Parameter"] = int(branch), "alpha"
        rows.append(scores[KEYS + ["Spearman_Corr", "Cosine_Similarity", "Valid_Samples", "Status", "Undefined_reason"]])
    return pd.concat(rows, ignore_index=True), mappings[0]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--artifact-root", type=Path)
    args = parser.parse_args()
    from .metrics import ensure_output
    ensure_output(args.output)
    scope, scores = load_scores()
    if args.artifact_root:
        inputs = frozen_inputs(args.artifact_root)
        scope = build_scope(inputs)
        extra, mapping = score_tsvelo(inputs)
        extra.to_csv(args.output / "tsvelo_organoid_gene_metrics.csv", index=False)
        mapping.to_csv(args.output / "tsvelo_organoid_gene_mapping.csv", index=False)
        scope.to_csv(args.output / "gene_scope.csv", index=False)
        scores = pd.concat([scores[~(scores.Dataset.eq("organoid") & scores.Method.eq("TSvelo"))], extra])
        scores = pd.concat([scores[scores.Dataset.eq(d) & scores.Gene.isin(candidate_genes(scope, d))]
                            for d in ("rpe1", "organoid")])
    for name, frame in summaries(scope, scores).items():
        frame.to_csv(args.output / f"{name}.csv", index=False)


if __name__ == "__main__":
    main()
