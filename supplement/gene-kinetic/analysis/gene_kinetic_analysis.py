#!/usr/bin/env python3
"""Reproducible gene-kinetic gene/cell edge analysis.

Inputs are the 15 gene-kinetic manifests, the fixed reference-30 gene edge
table, the Diaph3 supplement, and the 2026-09-20 canonical/gene-matched
cell-level CBDir/BATC edge tables.  The script uses gene-matched cell output
for DeepVelo, UniTVelo, and veloVI; all other methods use the canonical table.
It writes only compact CSV/Markdown analysis products.
"""

import argparse
import csv
import glob
import hashlib
import json
import math
import os
import re
import statistics
from collections import defaultdict
from pathlib import Path


GENE_METHODS = [
    "cellDancer", "DeepVelo", "DeepVelo_scv", "Graphvelo_scv", "MultiVelo",
    "scVelo_dyn", "TIvelo", "UniTVelo", "veloVI",
]

# display method -> (cell table, method id in that table)
CELL_SOURCE = {
    "cellDancer": ("canonical", "celldancer"),
    "DeepVelo": ("gene_matched", "deepvelo"),
    "DeepVelo_scv": ("canonical", "deepvelo2022"),
    "Graphvelo_scv": ("canonical", "graphvelo_scv"),
    "MultiVelo": ("canonical", "multivelo"),
    "scVelo_dyn": ("canonical", "scvelo_dynamic"),
    "TIvelo": ("canonical", "TIvelo"),
    "UniTVelo": ("gene_matched", "unitvelo"),
    "veloVI": ("gene_matched", "veloVI"),
}

ALIASES = {"scvelodynamic": "scVelo_dyn", "scvelodyn": "scVelo_dyn"}


def method_key(value):
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def canonical_method(value):
    key = method_key(value)
    if key in ALIASES:
        return ALIASES[key]
    for method in GENE_METHODS:
        if method_key(method) == key:
            return method
    return value


def number(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fields):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_text(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)


def dataset_for_case(case):
    prefixes = {
        "capo_": "capo_pre", "dentate_": "DentateGyrus_lamanno",
        "gastrulation_": "gastrulation", "mousebrain_": "mousebrain",
        "pancreas_": "pancreas", "retina_": "retina",
    }
    for prefix, dataset in prefixes.items():
        if case.startswith(prefix):
            return dataset
    raise ValueError(f"unknown manifest case: {case}")


def cell_dataset(dataset):
    return "capofibro" if dataset == "capo_pre" else dataset


def parse_edge_text(text):
    left, score = text.rsplit(":", 1)
    source, target = left.split("→", 1)
    return source.strip(), target.strip(), number(score.strip())


def average_rank(values, reverse=False):
    indexed = sorted(enumerate(values), key=lambda item: item[1], reverse=reverse)
    ranks = [0.0] * len(values)
    i = 0
    while i < len(indexed):
        j = i + 1
        while j < len(indexed) and indexed[j][1] == indexed[i][1]:
            j += 1
        rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[indexed[k][0]] = rank
        i = j
    return ranks


def pearson(x, y):
    if len(x) < 2:
        return None
    mx, my = statistics.fmean(x), statistics.fmean(y)
    dx = [v - mx for v in x]
    dy = [v - my for v in y]
    sx = math.sqrt(sum(v * v for v in dx))
    sy = math.sqrt(sum(v * v for v in dy))
    if sx == 0 or sy == 0:
        return None
    return sum(a * b for a, b in zip(dx, dy)) / (sx * sy)


def spearman(x, y):
    return pearson(average_rank(x), average_rank(y))


def correlation(x, y):
    return pearson(x, y), spearman(x, y)


def stats(values):
    values = [v for v in values if v is not None]
    if not values:
        return {"n": 0, "mean": None, "median": None, "sd": None, "min": None, "max": None,
                "n_positive": 0, "n_negative": 0, "positive_fraction": None,
                "negative_fraction": None, "mean_abs": None}
    n = len(values)
    return {"n": n, "mean": statistics.fmean(values), "median": statistics.median(values),
            "sd": statistics.pstdev(values), "min": min(values), "max": max(values),
            "n_positive": sum(v > 0 for v in values), "n_negative": sum(v < 0 for v in values),
            "positive_fraction": sum(v > 0 for v in values) / n,
            "negative_fraction": sum(v < 0 for v in values) / n,
            "mean_abs": statistics.fmean(abs(v) for v in values)}


def load_manifests(repo_root):
    manifests = []
    paths = sorted(glob.glob(os.path.join(repo_root, "supplement/gene-kinetic/figures/*_manifest.json")))
    if len(paths) != 15:
        raise ValueError(f"expected 15 manifests, found {len(paths)}")
    for path in paths:
        with open(path, encoding="utf-8") as handle:
            data = json.load(handle)
        data["manifest_path"] = os.path.relpath(path, repo_root)
        data["dataset"] = dataset_for_case(data["case"])
        data["edge_source_by_method"] = dict(zip(data["methods"], data.get("source_paths", [])))
        manifests.append(data)
    return manifests


def load_gene_records(repo_root):
    records = {}
    main = os.path.join(repo_root, "evaluation/gene/results/gene_cbdir_edge_details_reference30_all.csv")
    supplement = os.path.join(repo_root, "supplement/gene-kinetic/results/diaph3_nipc_to_nbl1.csv")
    for row in read_csv(main):
        key = (row["dataset"], canonical_method(row["method"]), row["gene"], row["source"], row["target"])
        records[key] = {"gene_cbdir": number(row.get("new_edge_score")),
                        "gene_evaluable_source_cells": row.get("evaluable_source_cells", ""),
                        "gene_status": "ok" if number(row.get("new_edge_score")) is not None else "missing",
                        "gene_source_table": "evaluation/gene/results/gene_cbdir_edge_details_reference30_all.csv"}
    for row in read_csv(supplement):
        key = (row["dataset"], canonical_method(row["method"]), row["gene"], row["source"], row["target"])
        if key not in records:
            records[key] = {"gene_cbdir": number(row.get("edge_score")),
                            "gene_evaluable_source_cells": row.get("evaluable_source_cells", ""),
                            "gene_status": row.get("status", "missing"),
                            "gene_source_table": "supplement/gene-kinetic/results/diaph3_nipc_to_nbl1.csv"}
    return records


def load_cell_tables(repo_root):
    paths = {
        "canonical": os.path.join(repo_root, "supplement/gene-kinetic/results/cell_level_cbdir_batc_rerun_20260920_v2/results/cell_metrics_edge_details_extended.csv"),
        "gene_matched": os.path.join(repo_root, "supplement/gene-kinetic/results/cell_level_cbdir_batc_gene_matched_20260920/results/cell_metrics_edge_details_extended.csv"),
    }
    records = {}
    all_rows = {}
    for table, path in paths.items():
        rows = read_csv(path)
        all_rows[table] = rows
        for row in rows:
            key = (row["dataset"], method_key(row["method"]), row["source"], row["target"], row["metric"])
            records[(table,) + key] = {
                "score": number(row.get("score")),
                "evaluable_source_cells": row.get("evaluable_source_cells", ""),
                "raw_nonfinite_source_cells": row.get("raw_nonfinite_source_cells", ""),
                "branch_winner_source_cells": number(row.get("branch_winner_source_cells")),
                "method_id": row["method"],
            }
    return paths, records, all_rows


def load_audits(repo_root):
    paths = {
        "canonical": os.path.join(repo_root, "supplement/gene-kinetic/results/cell_level_cbdir_batc_rerun_20260920_v2/results/cell_metrics_audit.csv"),
        "gene_matched": os.path.join(repo_root, "supplement/gene-kinetic/results/cell_level_cbdir_batc_gene_matched_20260920/results/cell_metrics_audit.csv"),
    }
    audits = {}
    for table, path in paths.items():
        for row in read_csv(path):
            audits[(table, row["dataset"], method_key(row["method"]))] = row
    return audits


def normalized_locator(value):
    value = value.replace("artifact://velocity/", "").replace("velocity:", "")
    return value.lstrip("/")


def build_source_mapping(manifests, audits):
    rows = []
    for manifest in manifests:
        for raw_method in manifest["methods"]:
            method = canonical_method(raw_method)
            table, source_method = CELL_SOURCE[method]
            c_dataset = cell_dataset(manifest["dataset"])
            audit = audits.get((table, c_dataset, method_key(source_method)), {})
            manifest_path = manifest["edge_source_by_method"].get(raw_method, "")
            result_path = audit.get("result_path", "")
            phase_layer_s = manifest.get("phase_contract", {}).get("layer_v_s", "")
            phase_layer_u = manifest.get("phase_contract", {}).get("layer_v_u", "")
            audit_layer = audit.get("velocity_layer_key", "")
            if phase_layer_s == audit_layer:
                layer_note = "exact_name"
            elif phase_layer_s == "velocity" and audit_layer in {"velocity_reconstructed", "velocity_gv"}:
                layer_note = "phase_velocity_alias_to_audit_layer"
            else:
                layer_note = "layer_name_mismatch"
            rows.append({
                "case": manifest["case"], "manifest_path": manifest["manifest_path"],
                "gene_dataset": manifest["dataset"], "cell_dataset": c_dataset,
                "gene": manifest["gene"], "method": method, "cell_table": table,
                "cell_method_id": source_method, "manifest_source_path": manifest_path,
                "audit_result_path": result_path,
                "source_path_match": bool(manifest_path and result_path and normalized_locator(manifest_path) == normalized_locator(result_path)),
                "phase_contract_velocity_s_key": phase_layer_s,
                "phase_contract_velocity_u_key": phase_layer_u,
                "audit_velocity_layer_key": audit_layer,
                "velocity_layer_contract_note": layer_note,
                "audit_status": audit.get("status", "missing_audit"),
                "result_sha256": audit.get("result_sha256", ""),
            })
    return rows


def parsed_manifest_edges(manifest):
    primary = tuple(manifest["edge"])
    parsed = {}
    for raw_method, edge_scores in manifest.get("edge_scores", {}).items():
        method = canonical_method(raw_method)
        for edge_text in edge_scores:
            source, target, score = parse_edge_text(edge_text)
            parsed[(method, source, target)] = score
    edge_keys = {primary}
    edge_keys.update((source, target) for _, source, target in parsed)
    return primary, parsed, edge_keys


def make_branch_denominators(all_rows):
    denominators = {}
    target_counts = {}
    for table, rows in all_rows.items():
        unique = {}
        for row in rows:
            if row["metric"] != "batc_candidate_source_mean":
                continue
            key = (row["dataset"], method_key(row["method"]), row["source"], row["target"])
            unique[key] = row
        by_source = defaultdict(list)
        for (dataset, method, source, target), row in unique.items():
            by_source[(dataset, method, source)].append(row)
        for key, rows_for_source in by_source.items():
            counts = [number(row.get("branch_winner_source_cells")) for row in rows_for_source]
            counts = [v for v in counts if v is not None]
            denominators[(table,) + key] = sum(counts) if counts else None
            target_counts[(table,) + key] = len(rows_for_source)
    return denominators, target_counts


def make_observations(manifests, gene_records, cell_records, source_mapping, denominators, target_counts):
    mapping_by_case_method = {(r["case"], r["method"]): r for r in source_mapping}
    rows = []
    for manifest in manifests:
        primary, manifest_scores, edge_keys = parsed_manifest_edges(manifest)
        for source, target in sorted(edge_keys):
            edge_role = "primary" if (source, target) == primary else "secondary_control"
            for raw_method in manifest["methods"]:
                method = canonical_method(raw_method)
                gene_key = (manifest["dataset"], method, manifest["gene"], source, target)
                gene = gene_records.get(gene_key, {})
                mapping = mapping_by_case_method[(manifest["case"], method)]
                table, source_method = CELL_SOURCE[method]
                c_dataset = mapping["cell_dataset"]
                base = (table, c_dataset, method_key(source_method), source, target)
                cbdir = cell_records.get(base + ("cbdir",), {})
                batc = cell_records.get(base + ("batc_candidate_source_mean",), {})
                denom_key = (table, c_dataset, method_key(source_method), source)
                winner = batc.get("branch_winner_source_cells")
                denom = denominators.get(denom_key)
                n_targets = target_counts.get(denom_key)
                rows.append({
                    "case": manifest["case"], "manifest_path": manifest["manifest_path"],
                    "dataset": manifest["dataset"], "cell_dataset": c_dataset,
                    "gene": manifest["gene"], "source": source, "target": target,
                    "edge_role": edge_role, "method": method,
                    "manifest_edge_score": manifest_scores.get((method, source, target)),
                    "gene_cbdir": gene.get("gene_cbdir"),
                    "gene_evaluable_source_cells": gene.get("gene_evaluable_source_cells", ""),
                    "gene_status": gene.get("gene_status", "missing"),
                    "gene_source_table": gene.get("gene_source_table", ""),
                    "cell_table": table, "cell_method_id": source_method,
                    "cell_source_mapping_status": "audit_ok_path_match" if mapping["audit_status"] == "ok" and mapping["source_path_match"] else "audit_or_path_check_failed",
                    "cell_result_path": mapping["audit_result_path"],
                    "cell_velocity_layer_key": mapping["audit_velocity_layer_key"],
                    "phase_contract_velocity_s_key": mapping["phase_contract_velocity_s_key"],
                    "velocity_layer_contract_note": mapping["velocity_layer_contract_note"],
                    "cell_cbdir": cbdir.get("score"),
                    "cell_batc_candidate_source_mean": batc.get("score"),
                    "cell_evaluable_source_cells": cbdir.get("evaluable_source_cells", "") or batc.get("evaluable_source_cells", ""),
                    "cell_raw_nonfinite_source_cells": cbdir.get("raw_nonfinite_source_cells", "") or batc.get("raw_nonfinite_source_cells", ""),
                    "branch_winner_source_cells": winner,
                    "branch_winner_total_source_cells": denom,
                    "branch_source_target_count": n_targets,
                    "branching_source": bool(n_targets and n_targets > 1),
                    "branch_winner_share": winner / denom if winner is not None and denom else None,
                })
    return rows


def aggregate_gene_edge_method(observations):
    groups = defaultdict(list)
    for row in observations:
        groups[(row["dataset"], row["cell_dataset"], row["source"], row["target"], row["edge_role"], row["method"])].append(row)
    out = []
    for key, rows in sorted(groups.items()):
        dataset, c_dataset, source, target, edge_role, method = key
        genes = [r["gene_cbdir"] for r in rows if r["gene_cbdir"] is not None]
        cell_scores = {}
        for field in ("cell_cbdir", "cell_batc_candidate_source_mean", "branch_winner_share"):
            vals = [r[field] for r in rows if r[field] is not None]
            cell_scores[field] = vals[0] if vals else None
            cell_scores[field + "_n_unique"] = len({round(v, 14) for v in vals})
        out.append({"dataset": dataset, "cell_dataset": c_dataset, "source": source, "target": target,
                    "edge_role": edge_role, "method": method, "n_genes": len(genes),
                    "genes": ";".join(sorted({r["gene"] for r in rows})),
                    "gene_mean": statistics.fmean(genes) if genes else None,
                    "gene_median": statistics.median(genes) if genes else None,
                    "gene_sd": statistics.pstdev(genes) if len(genes) > 1 else 0.0 if genes else None,
                    "gene_min": min(genes) if genes else None, "gene_max": max(genes) if genes else None,
                    "cell_cbdir": cell_scores["cell_cbdir"], "cell_cbdir_n_unique": cell_scores["cell_cbdir_n_unique"],
                    "cell_batc_candidate_source_mean": cell_scores["cell_batc_candidate_source_mean"],
                    "cell_batc_n_unique": cell_scores["cell_batc_candidate_source_mean_n_unique"],
                    "branch_winner_share": cell_scores["branch_winner_share"],
                    "branch_winner_share_n_unique": cell_scores["branch_winner_share_n_unique"],
                    "cell_evidence_status": "matched" if cell_scores["cell_cbdir"] is not None and cell_scores["cell_batc_candidate_source_mean"] is not None else "missing_edge",
                    "branching_source": rows[0]["branching_source"]})
    return out


def sign_match(a, b):
    return None if a is None or b is None else (a == 0 or b == 0 or (a > 0) == (b > 0))


def correlation_summary(observations, aggregated):
    result = []
    fields = [("cbdir", "cell_cbdir"), ("batc_candidate_source_mean", "cell_batc_candidate_source_mean")]
    scopes = [("primary", lambda r: r["edge_role"] == "primary"),
              ("secondary_control", lambda r: r["edge_role"] == "secondary_control"),
              ("all_displayed", lambda r: True)]
    for level, source_rows in (("raw_gene_method_edge", observations), ("gene_aggregated_edge_method", aggregated)):
        for scope, predicate in scopes:
            selected = [r for r in source_rows if predicate(r)]
            for metric, field in fields:
                gene_field = "gene_cbdir" if level.startswith("raw") else "gene_mean"
                pairs = [(r[gene_field], r[field]) for r in selected if r[gene_field] is not None and r[field] is not None]
                x, y = [a for a, _ in pairs], [b for _, b in pairs]
                p, s = correlation(x, y)
                result.append({"level": level, "scope": scope, "metric": metric, "centering": "none",
                               "n_pairs": len(pairs), "n_unique_edges": len({(r["dataset"], r["source"], r["target"]) for r in selected if r[gene_field] is not None and r[field] is not None}),
                               "n_unique_genes": len({r["gene"] for r in selected if r[gene_field] is not None and r[field] is not None}) if level.startswith("raw") else "",
                               "pearson": p, "spearman": s,
                               "sign_match_count": sum(sign_match(a, b) for a, b in pairs) if pairs else 0,
                               "sign_match_fraction": sum(sign_match(a, b) for a, b in pairs) / len(pairs) if pairs else None,
                               "interpretation": "exploratory_small_n" if len(pairs) < 30 else "descriptive_candidate_subset"})
                if level.startswith("gene_aggregated"):
                    groups = defaultdict(list)
                    for r in selected:
                        if r[gene_field] is not None and r[field] is not None:
                            groups[(r["dataset"], r["source"], r["target"])].append(r)
                    centered = []
                    for group_rows in groups.values():
                        gx = statistics.fmean(r[gene_field] for r in group_rows)
                        gy = statistics.fmean(r[field] for r in group_rows)
                        centered.extend((r[gene_field] - gx, r[field] - gy) for r in group_rows)
                    cx, cy = zip(*centered) if centered else ([], [])
                    p, s = correlation(list(cx), list(cy))
                    result.append({"level": level, "scope": scope, "metric": metric, "centering": "within_edge",
                                   "n_pairs": len(centered), "n_unique_edges": len(groups), "n_unique_genes": "",
                                   "pearson": p, "spearman": s, "sign_match_count": "", "sign_match_fraction": "",
                                   "interpretation": "exploratory_small_n" if len(groups) < 10 else "descriptive_candidate_subset"})
    return result


def correlation_by_method(aggregated):
    result = []
    fields = [("cbdir", "cell_cbdir"), ("batc_candidate_source_mean", "cell_batc_candidate_source_mean")]
    scopes = [("primary", lambda r: r["edge_role"] == "primary"),
              ("secondary_control", lambda r: r["edge_role"] == "secondary_control"),
              ("all_displayed", lambda r: True)]
    for scope, predicate in scopes:
        for method in sorted({r["method"] for r in aggregated}):
            selected = [r for r in aggregated if predicate(r) and r["method"] == method]
            for metric, field in fields:
                pairs = [(r["gene_mean"], r[field]) for r in selected if r["gene_mean"] is not None and r[field] is not None]
                x, y = [a for a, _ in pairs], [b for _, b in pairs]
                p, s = correlation(x, y)
                result.append({"scope": scope, "method": method, "metric": metric, "n_edges": len(pairs),
                               "pearson": p, "spearman": s,
                               "sign_match_count": sum(sign_match(a, b) for a, b in pairs) if pairs else 0,
                               "sign_match_fraction": sum(sign_match(a, b) for a, b in pairs) / len(pairs) if pairs else None,
                               "interpretation": "exploratory_small_n" if len(pairs) < 10 else "descriptive_candidate_subset"})
    return result


def edge_method_ranks(aggregated):
    out = []
    for scope, predicate in (("primary", lambda r: r["edge_role"] == "primary"),
                             ("secondary_control", lambda r: r["edge_role"] == "secondary_control")):
        selected = [r for r in aggregated if predicate(r)]
        groups = defaultdict(list)
        for row in selected:
            groups[(row["dataset"], row["cell_dataset"], row["source"], row["target"], row["edge_role"])].append(row)
        for key, rows in sorted(groups.items()):
            gene_rows = [r for r in rows if r["gene_mean"] is not None]
            gene_rank = dict(zip([id(r) for r in gene_rows], average_rank([r["gene_mean"] for r in gene_rows], reverse=True)))
            for metric, field in (("cbdir", "cell_cbdir"), ("batc_candidate_source_mean", "cell_batc_candidate_source_mean")):
                cell_rows = [r for r in rows if r["gene_mean"] is not None and r[field] is not None]
                cell_rank = dict(zip([id(r) for r in cell_rows], average_rank([r[field] for r in cell_rows], reverse=True)))
                for row in rows:
                    out.append({"scope": scope, "dataset": key[0], "cell_dataset": key[1], "source": key[2], "target": key[3],
                                "edge_role": key[4], "metric": metric, "method": row["method"], "genes": row["genes"],
                                "n_genes": row["n_genes"], "gene_mean": row["gene_mean"], "cell_score": row[field],
                                "gene_rank_high_is_1": gene_rank.get(id(row)), "cell_rank_high_is_1": cell_rank.get(id(row)),
                                "rank_delta_cell_minus_gene": (cell_rank.get(id(row)) - gene_rank.get(id(row))) if id(row) in cell_rank and id(row) in gene_rank else None,
                                "sign_match": sign_match(row["gene_mean"], row[field]),
                                "branch_winner_share": row["branch_winner_share"], "branching_source": row["branching_source"]})
    return out


def edge_correlation_summary(aggregated):
    groups = defaultdict(list)
    for row in aggregated:
        groups[(row["dataset"], row["cell_dataset"], row["source"], row["target"], row["edge_role"])].append(row)
    out = []
    for key, rows in sorted(groups.items()):
        for metric, field in (("cbdir", "cell_cbdir"), ("batc_candidate_source_mean", "cell_batc_candidate_source_mean")):
            pairs = [(r["gene_mean"], r[field]) for r in rows if r["gene_mean"] is not None and r[field] is not None]
            x, y = [a for a, _ in pairs], [b for _, b in pairs]
            p, s = correlation(x, y)
            out.append({"dataset": key[0], "cell_dataset": key[1], "source": key[2], "target": key[3], "edge_role": key[4],
                        "metric": metric, "n_methods": len(pairs), "pearson": p, "spearman": s,
                        "sign_match_count": sum(sign_match(a, b) for a, b in pairs) if pairs else 0,
                        "sign_match_fraction": sum(sign_match(a, b) for a, b in pairs) / len(pairs) if pairs else None,
                        "interpretation": "exploratory_small_n" if len(pairs) < 10 else "descriptive_candidate_subset"})
    return out


def edge_rank_summary(rank_rows):
    groups = defaultdict(list)
    for row in rank_rows:
        groups[(row["scope"], row["dataset"], row["cell_dataset"], row["source"], row["target"], row["metric"])].append(row)
    result = []
    for key, rows in sorted(groups.items()):
        valid = [r for r in rows if r["sign_match"] != "" and r["sign_match"] is not None]
        result.append({"scope": key[0], "dataset": key[1], "cell_dataset": key[2], "source": key[3], "target": key[4], "metric": key[5],
                       "method_rank_order_gene": ">".join(r["method"] for r in sorted(rows, key=lambda r: (r["gene_rank_high_is_1"] is None, r["gene_rank_high_is_1"] or 999))),
                       "method_rank_order_cell": ">".join(r["method"] for r in sorted(rows, key=lambda r: (r["cell_rank_high_is_1"] is None, r["cell_rank_high_is_1"] or 999))),
                       "n_methods_with_both": len(valid), "sign_match_count": sum(bool(r["sign_match"]) for r in valid),
                       "sign_match_fraction": sum(bool(r["sign_match"]) for r in valid) / len(valid) if valid else None,
                       "branching_source": rows[0]["branching_source"]})
    return result


def method_overall_summary(observations):
    result = []
    for scope, predicate in (("primary", lambda r: r["edge_role"] == "primary"),
                             ("secondary_control", lambda r: r["edge_role"] == "secondary_control"),
                             ("all_displayed", lambda r: True)):
        groups = defaultdict(list)
        for row in observations:
            if predicate(row):
                groups[row["method"]].append(row["gene_cbdir"])
        for method, values in sorted(groups.items()):
            summary = stats(values)
            result.append({"scope": scope, "method": method, **summary,
                           "candidate_scope_note": "small purposeful manifest set; descriptive only"})
    return result


def candidate_coverage(observations):
    groups = defaultdict(list)
    for row in observations:
        if row["edge_role"] == "primary":
            groups[(row["case"], row["dataset"], row["gene"], row["source"], row["target"])].append(row)
    out = []
    for key, rows in sorted(groups.items()):
        missing = sorted(r["method"] for r in rows if r["gene_cbdir"] is None)
        out.append({"case": key[0], "dataset": key[1], "gene": key[2], "source": key[3], "target": key[4],
                    "expected_methods": len(rows), "valid_gene_methods": sum(r["gene_cbdir"] is not None for r in rows),
                    "missing_gene_methods": ";".join(missing),
                    "cell_matched_methods": ";".join(sorted(r["method"] for r in rows if r["cell_cbdir"] is not None and r["cell_batc_candidate_source_mean"] is not None)),
                    "cell_match_count": sum(r["cell_cbdir"] is not None and r["cell_batc_candidate_source_mean"] is not None for r in rows)})
    return out


def render_report(output_dir, manifests, observations, aggregated, correlations, mappings):
    valid_gene = sum(r["gene_cbdir"] is not None for r in observations)
    matched_both = sum(r["cell_cbdir"] is not None and r["cell_batc_candidate_source_mean"] is not None for r in observations)
    path_ok = sum(r["source_path_match"] and r["audit_status"] == "ok" for r in mappings)
    layer_alias = sum(r["velocity_layer_contract_note"] == "phase_velocity_alias_to_audit_layer" for r in mappings)
    primary = [r for r in correlations if r["level"] == "gene_aggregated_edge_method" and r["scope"] == "primary" and r["centering"] == "none"]
    ptext = []
    for r in primary:
        ptext.append("- %s: n=%s, Pearson=%s, Spearman=%s, sign match=%s." % (
            r["metric"], r["n_pairs"],
            f"{r['pearson']:.3f}" if r["pearson"] is not None else "NA",
            f"{r['spearman']:.3f}" if r["spearman"] is not None else "NA",
            f"{r['sign_match_fraction']:.3f}" if r["sign_match_fraction"] is not None else "NA"))
    def edge_summary(dataset, source, target):
        rows = [r for r in aggregated if r["dataset"] == dataset and r["source"] == source and r["target"] == target and r["edge_role"] == "primary" and r["gene_mean"] is not None]
        result = {}
        for metric, field in (("cbdir", "cell_cbdir"), ("batc", "cell_batc_candidate_source_mean")):
            pairs = [(r["gene_mean"], r[field]) for r in rows if r[field] is not None]
            x, y = [a for a, _ in pairs], [b for _, b in pairs]
            p, s = correlation(x, y)
            result[metric] = {"n": len(pairs), "pearson": p, "spearman": s,
                              "sign": sum(sign_match(a, b) for a, b in pairs),
                              "rows": rows}
        return result
    gastr = edge_summary("gastrulation", "Erythroid2", "Erythroid3")
    capo = edge_summary("capo_pre", "G1", "S")
    pancreas = edge_summary("pancreas", "Pre-endocrine", "Beta")
    retina = edge_summary("retina", "Neuroblast", "RGC")
    ncald = edge_summary("DentateGyrus_lamanno", "ImmGranule1", "ImmGranule2")
    def fmt(value):
        return f"{value:.3f}" if value is not None else "NA"
    def method_triplet(summary, method):
        row = next(r for r in summary["cbdir"]["rows"] if r["method"] == method)
        return f"{method} (gene={fmt(row['gene_mean'])}, cell-CBDir={fmt(row['cell_cbdir'])}, BATC-candidate={fmt(row['cell_batc_candidate_source_mean'])})"
    primary_sign = {r["metric"]: r for r in correlations if r["level"] == "gene_aggregated_edge_method" and r["scope"] == "primary" and r["centering"] == "none"}
    control_sign = {r["metric"]: r for r in correlations if r["level"] == "gene_aggregated_edge_method" and r["scope"] == "secondary_control" and r["centering"] == "none"}
    pancreas_negative = ", ".join(r["method"] for r in pancreas["cbdir"]["rows"] if r["gene_mean"] < 0)
    retina_negative = ", ".join(r["method"] for r in retina["cbdir"]["rows"] if r["gene_mean"] < 0)
    text = f"""# Exploratory gene–cell edge analysis

This exploratory report is generated by `gene_kinetic_analysis.py` from the 15 phase-portrait manifests, the current reference-30 Gene edge table, the Diaph3 supplement, and the retained canonical/gene-matched Cell edge tables. It compares gene- and cell-level directionality on the displayed transitions. Gene scores include veloVI's coordinate-unit restoration; the phase portraits retain native coordinates. Input and output identities are in `provenance.json`.

## Source pairing

- 15 manifests produce 28 displayed gene edges and {len(observations)} gene-method-edge rows.
- Valid gene-CBDir values: {valid_gene}/{len(observations)}.
- Exact cell pairing for both native CBDir and `batc_candidate_source_mean`: {matched_both}/{len(observations)} rows.
- Source audit/path agreement: {path_ok}/{len(mappings)} case×method mappings.
- The phase contract names the spliced velocity layer `velocity`; the mapping table records {layer_alias} method-specific aliases (`velocity_reconstructed` or `velocity_gv`).
- DeepVelo, UniTVelo, and veloVI use gene-matched H5AD outputs. The other six methods use the retained canonical Cell outputs. The full mapping is in `source_mapping_audit.csv`.

## Primary-edge association after gene aggregation

To avoid pseudo-replication, genes sharing the same dataset/source/target/method are first averaged into one edge-method row. The values below compare gene-CBDir with the corresponding cell-level directionality score.

""" + "\n".join(ptext) + f"""

Within-edge-centered correlations are reported in `correlation_summary.csv`. They use method-centered residuals and should be treated as exploratory when the number of independent edges is small.
Method-specific correlations across edges are in `correlation_by_method.csv`; most have fewer than ten independent edges and are therefore exploratory.
The per-edge nine-method correlations used below are recorded in `edge_correlation_summary.csv`.

## Concrete edge-level findings (descriptive)

- Gastrulation Erythroid2→Erythroid3: the seven genes are aggregated by method before comparing the nine methods. Gene vs cell-CBDir has Spearman {fmt(gastr['cbdir']['spearman'])}, Pearson {fmt(gastr['cbdir']['pearson'])}; gene vs BATC candidate has Spearman {fmt(gastr['batc']['spearman'])}, Pearson {fmt(gastr['batc']['pearson'])}. Sign agreement is {gastr['cbdir']['sign']}/{gastr['cbdir']['n']} for cell-CBDir and {gastr['batc']['sign']}/{gastr['batc']['n']} for BATC candidate. Positive examples are {method_triplet(gastr, 'TIvelo')} and {method_triplet(gastr, 'cellDancer')}; negative examples are {method_triplet(gastr, 'DeepVelo_scv')} and {method_triplet(gastr, 'MultiVelo')}.
- capo G1→S: after aggregating the two genes, cell-CBDir sign agreement is {capo['cbdir']['sign']}/{capo['cbdir']['n']}. The gene vs cell-CBDir association is Pearson {fmt(capo['cbdir']['pearson'])}/Spearman {fmt(capo['cbdir']['spearman'])}; for BATC candidate it is Pearson {fmt(capo['batc']['pearson'])}/Spearman {fmt(capo['batc']['spearman'])}. These are descriptive associations from one edge and nine methods.
- Pancreas Pre-endocrine→Beta: negative Gene scores occur for {pancreas_negative}. Cell-CBDir spans {fmt(min(r['cell_cbdir'] for r in pancreas['cbdir']['rows']))} to {fmt(max(r['cell_cbdir'] for r in pancreas['cbdir']['rows']))}; BATC candidate means span {fmt(min(r['cell_batc_candidate_source_mean'] for r in pancreas['batc']['rows']))} to {fmt(max(r['cell_batc_candidate_source_mean'] for r in pancreas['batc']['rows']))}. This is direction-level decoupling on the displayed edge.
- Retina Neuroblast→RGC: negative Gene scores occur for {retina_negative}. Sign agreement is {retina['cbdir']['sign']}/{retina['cbdir']['n']} for cell-CBDir and {retina['batc']['sign']}/{retina['batc']['n']} for BATC candidate.
- Ncald ImmGranule1→ImmGranule2: {method_triplet(ncald, 'cellDancer')} has a positive Gene score and negative Cell scores on the displayed transition.
- After gene aggregation, primary-edge sign agreement is {primary_sign['cbdir']['sign_match_fraction']:.3f} for cell-CBDir and {primary_sign['batc_candidate_source_mean']['sign_match_fraction']:.3f} for BATC candidate, versus {control_sign['cbdir']['sign_match_fraction']:.3f} and {control_sign['batc_candidate_source_mean']['sign_match_fraction']:.3f} on secondary/control edges. The primary/control difference is descriptive and based on only seven and six independent displayed edges, respectively.

## How to read BATC candidate values

`batc_candidate_source_mean` is the mean candidate cosine before per-source branch maximization. `branch_winner_share` is the fraction of source cells for which the edge won among the candidate targets in the full extended table. Branching-source diagnostics are in `branching_winner_share.csv`.

## Cross-level findings

- Paired data summarize associations between gene-CBDir and native cell CBDir or pre-max candidate means on the predefined edges.
- Method rankings and sign agreement vary by edge; the per-edge rank and sign tables expose this heterogeneity.
- Branching edges can separate candidate directional similarity from branch selection, because winner share is reported independently.

"""
    write_text(os.path.join(output_dir, "gene_cell_edge_analysis.md"), text.rstrip() + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", default=os.path.abspath(os.path.join(os.path.dirname(__file__), "../../..")))
    parser.add_argument("--output-dir", default=None)
    args = parser.parse_args()
    repo_root = os.path.abspath(args.repo_root)
    output_dir = os.path.abspath(args.output_dir or os.path.join(repo_root, "supplement/gene-kinetic/analysis"))
    manifests = load_manifests(repo_root)
    genes = load_gene_records(repo_root)
    _, cells, all_cell_rows = load_cell_tables(repo_root)
    audits = load_audits(repo_root)
    source_mapping = build_source_mapping(manifests, audits)
    denominators, target_counts = make_branch_denominators(all_cell_rows)
    observations = make_observations(manifests, genes, cells, source_mapping, denominators, target_counts)
    aggregated = aggregate_gene_edge_method(observations)
    correlations = correlation_summary(observations, aggregated)
    rank_rows = edge_method_ranks(aggregated)
    rank_summary = edge_rank_summary(rank_rows)

    obs_fields = ["case", "manifest_path", "dataset", "cell_dataset", "gene", "source", "target", "edge_role", "method",
                  "manifest_edge_score", "gene_cbdir", "gene_evaluable_source_cells", "gene_status", "gene_source_table",
                  "cell_table", "cell_method_id", "cell_source_mapping_status", "cell_result_path", "cell_velocity_layer_key",
                  "phase_contract_velocity_s_key", "velocity_layer_contract_note",
                  "cell_cbdir", "cell_batc_candidate_source_mean", "cell_evaluable_source_cells", "cell_raw_nonfinite_source_cells",
                  "branch_winner_source_cells", "branch_winner_total_source_cells", "branch_source_target_count", "branching_source", "branch_winner_share"]
    write_csv(os.path.join(output_dir, "manifest_edge_observations.csv"), observations, obs_fields)
    write_csv(os.path.join(output_dir, "primary_edge_observations.csv"), [r for r in observations if r["edge_role"] == "primary"], obs_fields)
    mapping_fields = ["case", "manifest_path", "gene_dataset", "cell_dataset", "gene", "method", "cell_table", "cell_method_id",
                      "manifest_source_path", "audit_result_path", "source_path_match", "phase_contract_velocity_s_key",
                      "phase_contract_velocity_u_key", "audit_velocity_layer_key", "velocity_layer_contract_note", "audit_status", "result_sha256"]
    write_csv(os.path.join(output_dir, "source_mapping_audit.csv"), source_mapping, mapping_fields)
    agg_fields = ["dataset", "cell_dataset", "source", "target", "edge_role", "method", "n_genes", "genes", "gene_mean", "gene_median", "gene_sd",
                  "gene_min", "gene_max", "cell_cbdir", "cell_cbdir_n_unique", "cell_batc_candidate_source_mean", "cell_batc_n_unique",
                  "branch_winner_share", "branch_winner_share_n_unique", "cell_evidence_status", "branching_source"]
    write_csv(os.path.join(output_dir, "aggregated_edge_method_pairs.csv"), aggregated, agg_fields)
    method_fields = ["scope", "method", "n", "mean", "median", "sd", "min", "max", "n_positive", "n_negative",
                     "positive_fraction", "negative_fraction", "mean_abs", "candidate_scope_note"]
    write_csv(os.path.join(output_dir, "method_overall_summary.csv"), method_overall_summary(observations), method_fields)
    coverage_fields = ["case", "dataset", "gene", "source", "target", "expected_methods", "valid_gene_methods",
                       "missing_gene_methods", "cell_matched_methods", "cell_match_count"]
    write_csv(os.path.join(output_dir, "candidate_coverage_summary.csv"), candidate_coverage(observations), coverage_fields)
    corr_fields = ["level", "scope", "metric", "centering", "n_pairs", "n_unique_edges", "n_unique_genes", "pearson", "spearman",
                   "sign_match_count", "sign_match_fraction", "interpretation"]
    write_csv(os.path.join(output_dir, "correlation_summary.csv"), correlations, corr_fields)
    by_method_fields = ["scope", "method", "metric", "n_edges", "pearson", "spearman", "sign_match_count", "sign_match_fraction", "interpretation"]
    write_csv(os.path.join(output_dir, "correlation_by_method.csv"), correlation_by_method(aggregated), by_method_fields)
    edge_corr_fields = ["dataset", "cell_dataset", "source", "target", "edge_role", "metric", "n_methods", "pearson", "spearman",
                        "sign_match_count", "sign_match_fraction", "interpretation"]
    write_csv(os.path.join(output_dir, "edge_correlation_summary.csv"), edge_correlation_summary(aggregated), edge_corr_fields)
    rank_fields = ["scope", "dataset", "cell_dataset", "source", "target", "edge_role", "metric", "method", "genes", "n_genes", "gene_mean",
                   "cell_score", "gene_rank_high_is_1", "cell_rank_high_is_1", "rank_delta_cell_minus_gene", "sign_match", "branch_winner_share", "branching_source"]
    write_csv(os.path.join(output_dir, "edge_method_rank_sign.csv"), rank_rows, rank_fields)
    rank_summary_fields = ["scope", "dataset", "cell_dataset", "source", "target", "metric", "method_rank_order_gene", "method_rank_order_cell",
                           "n_methods_with_both", "sign_match_count", "sign_match_fraction", "branching_source"]
    write_csv(os.path.join(output_dir, "edge_rank_summary.csv"), rank_summary, rank_summary_fields)
    branch_rows = [r for r in observations if r["branching_source"]]
    branch_fields = obs_fields
    write_csv(os.path.join(output_dir, "branching_winner_share.csv"), branch_rows, branch_fields)
    primary_pairs = [r for r in observations if r["edge_role"] == "primary" and r["gene_cbdir"] is not None and r["cell_cbdir"] is not None and r["cell_batc_candidate_source_mean"] is not None]
    write_csv(os.path.join(output_dir, "cell_match_pairs.csv"), primary_pairs, obs_fields)
    render_report(output_dir, manifests, observations, aggregated, correlations, source_mapping)
    inputs = [os.path.join(repo_root, manifest["manifest_path"]) for manifest in manifests]
    inputs += [os.path.join(repo_root, "evaluation/gene/results/gene_cbdir_edge_details_reference30_all.csv"),
               os.path.join(repo_root, "supplement/gene-kinetic/results/diaph3_nipc_to_nbl1.csv")]
    for folder in ["cell_level_cbdir_batc_rerun_20260920_v2", "cell_level_cbdir_batc_gene_matched_20260920"]:
        inputs += [os.path.join(repo_root, "supplement/gene-kinetic/results", folder, "results", name)
                   for name in ["cell_metrics_edge_details_extended.csv", "cell_metrics_audit.csv"]]
    inputs.append(os.path.abspath(__file__))
    outputs = sorted(glob.glob(os.path.join(output_dir, "*.csv"))) + [os.path.join(output_dir, "gene_cell_edge_analysis.md")]
    provenance = {"analysis": "exploratory gene–cell directionality cases",
                  "inputs_sha256": {os.path.relpath(path, repo_root): hashlib.sha256(Path(path).read_bytes()).hexdigest() for path in inputs},
                  "outputs_sha256": {os.path.basename(path): hashlib.sha256(Path(path).read_bytes()).hexdigest() for path in outputs}}
    write_text(os.path.join(output_dir, "provenance.json"), json.dumps(provenance, indent=2) + "\n")
    print(json.dumps({"manifests": len(manifests), "displayed_edges": len({(r['dataset'], r['source'], r['target']) for r in observations}),
                      "observations": len(observations), "valid_gene": sum(r['gene_cbdir'] is not None for r in observations),
                      "matched_both_cell_metrics": sum(r['cell_cbdir'] is not None and r['cell_batc_candidate_source_mean'] is not None for r in observations),
                      "source_mappings": len(source_mapping), "output_dir": output_dir}, indent=2))


if __name__ == "__main__":
    main()
