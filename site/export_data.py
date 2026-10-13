#!/usr/bin/env python3
"""Export existing benchmark scores for static browsing; never recompute metrics."""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
CELL = "evaluation/cell/results/cell_metrics.csv"
GENE = "evaluation/gene/results/gene_cbdir_reference30_all.csv"
COSINE = "evaluation/gene/results/simulation/gene_cbdir_cossim_reference.csv"
KINETICS = "evaluation/kinetics_recovery/results/central_results/"
PLOT = "evaluation/kinetics_recovery/figures/fig6/plot_data/"
EXPANDED = "evaluation/kinetics_recovery/results/expanded_reference/"
CONTENT = "site/content/website_content.json"
# Existing Gene result identifiers -> method registry identities. No display names here.
GENE_IDS = {
    "scVelo_dyn": "C03", "MultiVelo": "C04", "UniTVelo": "C05",
    "TIvelo": "C08", "veloVI": "C12", "cellDancer": "C15",
    "DeepVelo": "C16", "DeepVelo_scv": "C18", "Graphvelo_scv": "C19",
}
SIM_IDS = {"murk": "sim_only_murk", "branching": "sim_only_branching", "all": "sim_hybrid"}


def numeric(value):
    """Preserve finite source precision; missing/nonfinite values stay missing."""
    if value in (None, "", "NA", "NaN", "nan"):
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def keyed(rows, columns):
    result = {}
    for row in rows:
        key = tuple(row[column] for column in columns)
        if key in result:
            raise ValueError(f"Duplicate key {key}")
        result[key] = row
    return result


def validate_visuals(root, directory):
    """Reject stale source data or damaged generated displays before deployment."""
    index = json.loads((directory / 'index.json').read_text())
    if index['schema_version'] != 1:
        raise ValueError('Unsupported visual data version')
    for source in index['sources']:
        if hashlib.sha256((root / source['path']).read_bytes()).hexdigest() != source['sha256']:
            raise ValueError(f"Visual source changed: {source['path']}")
    for entry in index['files']:
        path = (directory / entry['path']).resolve()
        if not path.is_relative_to(directory.resolve()):
            raise ValueError('Visual asset path escapes the display directory')
        if hashlib.sha256(path.read_bytes()).hexdigest() != entry['sha256']:
            raise ValueError(f"Visual asset changed: {entry['path']}")
    return index


def export_data(root=ROOT, output=None):
    root = Path(root).resolve()
    output = Path(output).resolve() if output else root / "site/public/data"
    # Validate everything before writing the generated display projection.
    if output == root or output in root.parents:
        raise ValueError("Output must not be the repository or one of its parents")
    if output.is_relative_to(root) and not output.is_relative_to(root / "site"):
        raise ValueError("Repository outputs must stay under site/")
    sources = []
    media_copies = []

    def register_file(relative, download=True):
        existing = next((s for s in sources if s["path"] == relative), None)
        if existing is not None:
            return existing
        path = root / relative
        source = {"path": relative, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                  "bytes": path.stat().st_size,
                  "download": f"data/files/{relative}" if download else None}
        sources.append(source)
        return source

    def read(relative, download=True):
        path = root / relative
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t" if path.suffix == ".tsv" else ",")
            rows = list(reader)
        register_file(relative, download).update(rows=len(rows), columns=reader.fieldnames)
        return rows

    def media_asset(relative, expected_sha256):
        if not (root / relative).resolve().is_relative_to(root):
            raise ValueError("Analysis image path escapes the repository")
        source = register_file(relative, download=False)
        if source["sha256"] != expected_sha256:
            raise ValueError(f"Analysis image changed: {relative}")
        target = "media/" + relative
        media_copies.append((root / relative, target))
        return {"path": "data/" + target, "source": relative, "sha256": source["sha256"]}

    methods = read("methods/configurations.tsv")
    datasets = read("data/datasets.tsv")
    by_method = {key[0]: row for key, row in keyed(methods, ["config_id"]).items()}
    by_dataset = {key[0]: row for key, row in keyed(datasets, ["dataset_id"]).items()}
    inputs = read("evaluation/cell/inputs.tsv")
    input_keys = keyed(inputs, ["evaluation_dataset", "operational_id"])
    supplementary = keyed(read("evaluation/cell/results/cell_metrics_supplementary.csv"), ["dataset", "method"])
    cell_rows = read(CELL)
    if set(keyed(cell_rows, ["dataset", "method"])) != set(input_keys):
        raise ValueError("Cell scores and frozen inputs have different keys")
    if set(supplementary) != set(input_keys):
        raise ValueError("Cell supplementary scores and frozen inputs have different keys")
    cell = []
    for row in cell_rows:
        binding = input_keys[(row["dataset"], row["method"])]
        by_method[binding["config_id"]]
        by_dataset[binding["dataset_id"]]
        cell.append({"dataset": binding["dataset_id"], "method": binding["config_id"],
                     "source": CELL, "n_cells": int(binding["n_common"]), "n_reference": int(binding["n_reference"]),
                     **{field: numeric(row[field]) for field in ["consistency", "iccoh", "cbdir", "transcosine", "batc"]},
                     "details": supplementary[(row["dataset"], row["method"])]})

    gene_rows = read(GENE)
    gene_summary_path = 'evaluation/gene/results/gene_cbdir_summary_by_dataset_method.csv'
    gene_summary = [{'dataset': SIM_IDS.get(r['dataset'], r['dataset'].lower()),
                     'method': GENE_IDS[r['method']], 'cbdir': numeric(r['new_mean']),
                     'n_genes': int(r['n_genes']), 'source': gene_summary_path}
                    for r in read(gene_summary_path)]
    gene_keys = keyed(gene_rows, ["dataset", "method", "gene"])
    cosine_rows = read(COSINE)
    cosine = keyed(cosine_rows, ["Dataset", "Method", "Gene"])
    expected_cosine = {key for key in gene_keys if key[0] in SIM_IDS}
    if set(cosine) != expected_cosine:
        raise ValueError("Simulation CosSim and reference-30 Gene scores have different keys")
    gene = []
    for row in gene_rows:
        dataset = SIM_IDS.get(row["dataset"], row["dataset"].lower())
        method = GENE_IDS[row["method"]]
        by_dataset[dataset]
        by_method[method]
        cs = cosine.get((row["dataset"], row["method"], row["gene"]))
        gene.append({"dataset": dataset, "method": method, "gene": row["gene"],
                     "cbdir": numeric(row["new_cbdir"]), "cossim": numeric(cs["CosSim"]) if cs else None,
                     "n_edges": numeric(row["n_evaluable_edges"]), "n_cells": numeric(row["n_evaluable_source_cells"]),
                     "status": row["status"], "source": GENE})

    kinetics_methods = {name: identifier for name, identifier in GENE_IDS.items()
                        if name in {"cellDancer", "DeepVelo"}}
    kinetics_methods.update({"DeepKINET": "C13", "TSvelo": "C14",
                            "Dynamo": "dynamo_labelled", "Dynamo(label data)": "dynamo_labelled"})
    # Labelled Dynamo is a distinct control, never silently joined to C07 RNA_seq.
    kinetics = []
    full = read(KINETICS + "gene_metrics.csv")
    keyed(full, ["Method", "Dataset", "Branch", "Parameter", "Gene"])
    for row in full:
        kinetics.append({"dataset": "sceu_" + row["Dataset"], "method": kinetics_methods[row["Method"]],
                         "branch": row["Branch"], "parameter": row["Parameter"], "gene": row["Gene"],
                         "score": numeric(row["Spearman_Corr"]), "n_samples": int(row["Valid_Samples"]),
                         "scope": "full", "source": KINETICS + "gene_metrics.csv"})
    for parameter in ["alpha", "gamma"]:
        source = PLOT + f"rpe1_{parameter}_heatmap.csv"
        for row in read(source):
            for gene_name, score in row.items():
                if gene_name == "Method":
                    continue
                kinetics.append({"dataset": "sceu_rpe1", "method": kinetics_methods[row["Method"]],
                                 "branch": "", "parameter": parameter, "gene": gene_name,
                                 "score": numeric(score), "n_samples": None, "scope": "selected", "source": source})
    for branch in ["1", "2"]:
        source = PLOT + f"organoid_branch{branch}_gene_metrics.csv"
        for row in read(source):
            method_key = row["Method"].removesuffix(f"_branch{branch}")
            kinetics.append({"dataset": "sceu_organoid", "method": kinetics_methods[method_key],
                             "branch": branch, "parameter": row["Parameter"].lower(), "gene": row["Gene"],
                             "score": numeric(row["Correlation"]), "n_samples": None,
                             "scope": "selected", "source": source})
    keyed(kinetics, ["scope", "dataset", "branch", "parameter", "method", "gene"])
    components = read(KINETICS + "component_means.csv")
    keyed(components, ["Method", "Dataset", "Branch", "Parameter"])
    summary = [{"dataset": "sceu_" + row["Dataset"], "method": kinetics_methods[row["Method"]],
                "branch": row["Branch"], "parameter": row["Parameter"], "score": numeric(row["Mean"]),
                "n_genes": int(row["N"]), "scope": "full", "comparison": "available",
                "source": KINETICS + "component_means.csv"} for row in components]
    read(KINETICS + "system_means.csv")

    # Independent supplementary analysis: restrict canonical rows to its frozen
    # gene universe and add only the registered TSvelo organoid alpha table.
    manifest_path = EXPANDED + "manifest.json"
    manifest = json.loads((root / manifest_path).read_text())
    register_file(manifest_path)
    if hashlib.sha256((root / KINETICS / "gene_metrics.csv").read_bytes()).hexdigest() != manifest["canonical_gene_metrics_sha256"]:
        raise ValueError("Canonical kinetics scores no longer match the expanded manifest")
    for relative, expected in manifest["outputs_sha256"].items():
        if hashlib.sha256((root / relative).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Expanded source hash mismatch: {relative}")
    gene_scope = read(EXPANDED + "gene_scope.csv")
    keyed(gene_scope, ["Dataset", "Gene"])
    universe = {ds: {r["Gene"] for r in gene_scope if r["Dataset"] == ds and r["Included"] == "True"}
                for ds in ["rpe1", "organoid"]}
    if len(universe["rpe1"]) != manifest["scope"]["rpe1"] or len(universe["organoid"]) != manifest["scope"]["organoid_each_branch"]:
        raise ValueError("Expanded reference gene scope disagrees with its manifest")
    extra = read(EXPANDED + "tsvelo_organoid_gene_metrics.csv")
    read(EXPANDED + "tsvelo_organoid_gene_mapping.csv")
    for row in extra:
        if (row["Method"], row["Dataset"], row["Parameter"]) != ("TSvelo", "organoid", "alpha") or row["Branch"] not in {"1", "2"}:
            raise ValueError("Expanded TSvelo scores must be organoid alpha, branches 1/2")
    keyed(full + extra, ["Method", "Dataset", "Branch", "Parameter", "Gene"])
    for row in full + extra:
        if row["Gene"] not in universe[row["Dataset"]]:
            continue
        is_extra = row["Method"] == "TSvelo" and row["Dataset"] == "organoid"
        kinetics.append({"dataset": "sceu_" + row["Dataset"], "method": kinetics_methods[row["Method"]],
                         "branch": row["Branch"], "parameter": row["Parameter"], "gene": row["Gene"],
                         "score": numeric(row["Spearman_Corr"]), "n_samples": int(row["Valid_Samples"]),
                         "status": row.get("Undefined_reason") or row.get("Status") or "",
                         "scope": "expanded", "source": EXPANDED + "tsvelo_organoid_gene_metrics.csv" if is_extra else KINETICS + "gene_metrics.csv"})
    keyed(kinetics, ["scope", "dataset", "branch", "parameter", "method", "gene"])
    coverage = read(EXPANDED + "coverage_summary.csv")
    keyed(coverage, ["Panel", "Method"])
    panels = {}
    for row in coverage:
        branch = str(int(float(row["Branch"]))) if row["Branch"] else ""
        binding = {"dataset": "sceu_" + row["Dataset"], "branch": branch, "parameter": row["Parameter"]}
        if row["Panel"] in panels and panels[row["Panel"]] != binding:
            raise ValueError("Expanded panel has inconsistent dataset/branch/parameter keys")
        panels[row["Panel"]] = binding
        summary.append({**binding, "method": kinetics_methods[row["Method"]], "scope": "expanded",
                        "comparison": "available", "score": numeric(row["Mean"]), "n_genes": int(row["Valid_N"]),
                        "candidate_n": int(row["Candidate_N"]), "present_n": int(row["Present_N"]),
                        "undefined_n": int(row["Undefined_N"]), "unavailable_n": int(row["Unavailable_N"]),
                        "source": EXPANDED + "coverage_summary.csv"})
    common = read(EXPANDED + "common_gene_summary.csv")
    keyed(common, ["Panel", "Comparison_scope", "Method"])
    membership = read(EXPANDED + "common_gene_membership.csv")
    keyed(membership, ["Panel", "Comparison_scope", "Gene"])
    common_genes = {}
    for row in membership:
        common_genes.setdefault((row["Panel"], row["Comparison_scope"]), []).append(row["Gene"])
    for row in common:
        members = common_genes[(row["Panel"], row["Comparison_scope"])]
        if len(members) != int(row["Common_genes_N"]) or len(members) != int(row["Method_valid_N"]):
            raise ValueError("Common-gene membership and summary counts disagree")
        summary.append({**panels[row["Panel"]], "method": kinetics_methods[row["Method"]],
                        "scope": "expanded", "comparison": row["Comparison_scope"], "score": numeric(row["Mean"]),
                        "n_genes": int(row["Method_valid_N"]), "source": EXPANDED + "common_gene_summary.csv"})
    common_membership = [{**panels[panel], "comparison": group, "genes": genes}
                         for (panel, group), genes in common_genes.items()]
    for relative in manifest["outputs_sha256"]:
        if relative.startswith("evaluation/kinetics_recovery/figures/"):
            register_file(relative)
    order_path = "evaluation/kinetics_recovery/plotting/figure6.py"
    order = next(ast.literal_eval(node.value) for node in ast.parse((root / order_path).read_text()).body
                 if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "PAPER_METHOD_ORDER" for t in node.targets))
    colors = next(ast.literal_eval(node.value) for node in ast.parse((root / order_path).read_text()).body
                  if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'TRAJECTORY_COLORS' for t in node.targets))
    kinetics_colors = {kinetics_methods[name]: color for name, color in colors.items()}
    register_file(order_path, download=False)
    kinetics_order = [kinetics_methods[name] for name in order]
    # Only public-facing fields from the artifact registry are exposed.
    artifact_rows = read("artifacts.tsv", download=False)
    artifacts = [{key: row[key] for key in ["artifact_id", "category", "sha256", "status", "uri"]}
                 for row in artifact_rows if row["status"] == "external"]
    release = {row['artifact_id']: row for row in read('data/release_map.tsv', download=False)}
    h5ad = [row for row in artifacts if release[row['artifact_id']]['member_path'].endswith('.h5ad')]
    for artifact in artifacts:
        if artifact["uri"] and not artifact["uri"].startswith("https://"):
            raise ValueError("Artifact downloads must use public HTTPS URIs")
    for entry in kinetics + summary:
        by_dataset[entry["dataset"]]
        if entry["method"] != "dynamo_labelled":
            by_method[entry["method"]]
    for relative in ["evaluation/cell/README.md", "evaluation/gene/README.md", "evaluation/kinetics_recovery/README.md", EXPANDED + "README.md"]:
        register_file(relative)
    content = json.loads((root / CONTENT).read_text())
    if content.get("schema_version") != 1:
        raise ValueError("Unsupported website content schema")
    if set(content["datasets"]) != set(by_dataset) or set(content["methods"]) != set(by_method):
        raise ValueError("Website content must cover the exact dataset and method registries")
    for dataset_id, scenario in content["simulations"].items():
        registered = by_dataset[dataset_id]
        if (scenario["cells"], scenario["genes"]) != (int(registered["n_obs"]), int(registered["n_vars"])):
            raise ValueError("Simulation content dimensions disagree with the registered input")
    register_file(CONTENT, download=False)
    register_file("site/content/README.md", download=False)
    cell_config_path = "evaluation/cell/config.json"
    cell_config = json.loads((root / cell_config_path).read_text())
    register_file(cell_config_path, download=False)
    ensemble_dataset_ids = {name: entry["dataset_id"] for name, entry in cell_config["datasets"].items()}
    ensemble_dataset_ids["dentate_lamanno"] = ensemble_dataset_ids["DentateGyrus_lamanno"]
    ensemble_scores = []
    for row in read("ensemble/results/scores.csv", download=False):
        ensemble_scores.append({"dataset": ensemble_dataset_ids[row["dataset"]], "variant": row["variant"],
                                "label": row["label"], "cbdir": numeric(row["cbdir"]), "batc": numeric(row["batc"])})
    scores_by_key = keyed(ensemble_scores, ["dataset", "variant"])
    ensemble_metadata_path = "ensemble/results/metadata.json"
    ensemble_metadata = json.loads((root / ensemble_metadata_path).read_text())
    register_file(ensemble_metadata_path, download=False)
    variants = ensemble_metadata["variants"]
    expected = {(ensemble_dataset_ids[name], variant) for name in ensemble_metadata["datasets"] for variant in variants}
    if set(scores_by_key) != expected or len(expected) != 28:
        raise ValueError("Ensemble scores must cover both variants on all 14 registered Cell datasets")
    ensemble_rankings = []
    for row in read("ensemble/results/rankings_vs_20.csv", download=False):
        dataset_id = ensemble_dataset_ids[row["dataset"]]
        if numeric(row["score"]) != scores_by_key[(dataset_id, row["variant"])][row["metric"]]:
            raise ValueError("Ensemble ranking score differs from its frozen score")
        ensemble_rankings.append({"dataset": dataset_id, "variant": row["variant"], "metric": row["metric"],
                                  "score": numeric(row["score"]), "rank": int(row["rank"]),
                                  "out_of": int(row["out_of"]), "tie_count": int(row["tie_count"])})
    if set(keyed(ensemble_rankings, ["dataset", "variant", "metric"])) != {(*key, metric) for key in expected for metric in ["cbdir", "batc"]}:
        raise ValueError("Ensemble rankings must cover both metrics for all frozen scores")
    ensemble = {"scores": ensemble_scores, "rankings": ensemble_rankings,
                "macro_average": ensemble_metadata["macro_average"], "variants": variants,
                "methods": ensemble_metadata["methods"]}
    # Frozen post-hoc stream panels are separate from the 20 Cell configurations.
    stream_artifacts = {r["path"]: r for r in artifact_rows if r["category"] == "ensemble-stream"}
    ensemble_streams = []
    for name in ensemble_metadata["datasets"]:
        for variant, folder in [(variants[0], "mean4"), (variants[1], "row_geometric_median4")]:
            relative = f"ensemble/results/figures/streams/{folder}/{name}.png"
            registered = stream_artifacts[relative]
            ensemble_streams.append({"dataset": ensemble_dataset_ids[name], "variant": variant,
                                     **media_asset(relative, registered["sha256"])})
    if set(keyed(ensemble_streams, ["dataset", "variant"])) != expected:
        raise ValueError("Ensemble images and scores must have identical keys")
    ensemble["streams"] = ensemble_streams

    commit = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    build_state = "working-tree-preview" if subprocess.check_output(["git","-C",str(root),"status","--porcelain"],text=True).strip() else "committed"
    downloads = {'tables': [], 'figures': [], 'metadata': [], 'audit': []}
    # Expose only publication result tables and summary figures as downloads.
    for item in sources:
        item['download'] = None
    table_downloads = [
        (CELL, 'Cell-level scores'), (GENE, 'Gene-level CBDir (reference 30-NN)'),
        (COSINE, 'Simulation CosSim (frozen table)'),
        (gene_summary_path, 'Gene-level dataset means'),
        (KINETICS + 'gene_metrics.csv', 'Canonical kinetic-parameter scores'),
        (KINETICS + 'component_means.csv', 'Canonical kinetic-parameter means'),
        (EXPANDED + 'coverage_summary.csv', 'Expanded reference: coverage and means'),
        (EXPANDED + 'common_gene_summary.csv', 'Expanded reference: common-gene means'),
        (EXPANDED + 'tsvelo_organoid_gene_metrics.csv', 'Expanded reference: TSvelo organoid alpha'),
    ]
    for relative, label in table_downloads:
        entry = next(s for s in sources if s['path'] == relative)
        entry['download'] = f'data/files/{relative}'
        downloads['tables'].append({'label': label, **entry})
    metadata_downloads = [
        ('data/input_release_manifest.tsv', 'Input release manifest',
         '28 source files with dataset roles, layer contracts, SHA-256 values, and pending publication status.'),
        ('data/kinetics_reference_gene_usage/README.md', 'Expanded reference-gene guide',
         'Scope and valid/undefined/unavailable status definitions for the 250/270 gene tables.'),
        ('data/kinetics_reference_gene_usage/manifest.json', 'Expanded reference-gene manifest',
         'Source commit, output hashes, coverage checks, and upstream reference provenance.'),
        ('data/kinetics_reference_gene_usage/rpe1_reference_gene_usage.csv', 'RPE1 reference-gene usage',
         '250 retained reference genes, Figure 6 display flags, and parameter coverage status.'),
        ('data/kinetics_reference_gene_usage/organoid_reference_gene_usage.csv', 'Organoid reference-gene usage',
         '270 retained genes shared by both branches, Figure 6 display flags, and method status.'),
    ]
    for relative, label, description in metadata_downloads:
        entry = register_file(relative)
        downloads['metadata'].append({'label': label, 'description': description, **entry})
    audit_downloads = [
        ("evaluation/cell/results/cell_metrics_supplementary.csv", "Cell supplementary metrics and evaluability counts", "Supplementary Cell metrics and valid-sample fields; it does not replace the five primary scores."),
        ("evaluation/cell/results/cell_metrics_edge_details.csv", "Cell edge-level audit", "Edge-level Cell evaluation details under the fixed reference contract."),
        ("evaluation/gene/results/gene_cbdir_edge_details_reference30_all.csv", "Gene edge-level audit", "Gene-level edge details underlying the reference-30 CBDir rows."),
        ("evaluation/gene/results/neighbor_audit_all.csv", "Gene neighbor audit", "Neighbor-set provenance and overlap checks for Gene-level evaluation."),
        ("ensemble/results/scores.csv", "Post-hoc ensemble scores", "Mean-4 and RowGeometricMedian-4 scores across the 14 real datasets."),
        ("ensemble/results/rankings_vs_20.csv", "Post-hoc ensemble rankings", "Independent per-dataset CBDir/BATC comparisons with the 20 registered configurations."),
        ("supplement/velovi/results/table_s6.csv", "veloVI coordinate-sensitivity summary", "A/B/C coordinate diagnostics and C − B differences across 14 datasets; primary Cell scoring is condition B."),
        ("evaluation/cell/results/cell_metrics_audit.csv", "Cell source and status audit", "Logical source identity, velocity layers, invalid-value handling and status for each method/dataset."),
        ("evaluation/cell/reference_neighbors/manifest.json", "Cell reference-neighbor manifest", "Hashes and construction settings for the frozen dataset-level neighbor assets."),
        ("evaluation/gene/results/provenance.json", "Gene result provenance", "Input and neighbor contracts for the frozen Gene evaluation."),
        ("ensemble/results/metadata.json", "Ensemble protocol and provenance", "Frozen input roles, transition settings and macro averages for both variants."),
        ("supplement/velovi/results/full_dataset_three_condition.csv", "veloVI coordinate diagnostic: full values", "42 dataset/condition rows; this post-hoc coordinate comparison does not replace Cell scores."),
        ("supplement/velovi/results/manifest.json", "veloVI coordinate diagnostic manifest", "Conditions, input hashes and output hashes for the coordinate-sensitivity analysis."),
    ]
    for relative, label, description in audit_downloads:
        entry = next((s for s in sources if s['path'] == relative), None)
        if entry is None:
            entry = register_file(relative)
        entry['download'] = f'data/files/{relative}'
        downloads['audit'].append({'label': label, 'description': description, **entry})
    figure_manifest_path = 'figures/fig7_composite_manifest.json'
    figure_manifest = json.loads((root / figure_manifest_path).read_text())
    register_file(figure_manifest_path, download=False)
    figure = figure_manifest['output']
    if hashlib.sha256((root / figure['path']).read_bytes()).hexdigest() != figure['sha256']:
        raise ValueError('Figure 7 does not match its composite manifest')
    entry = register_file(figure['path'])
    downloads['figures'].append({'label': 'Benchmark summary (Figure 7)', **entry})
    for relative, label in [
        ('rpe1_boxplot.png', 'Expanded reference: RPE1'),
        ('organoid_branch1_boxplot.png', 'Expanded reference: organoid branch 1'),
        ('organoid_branch2_boxplot.png', 'Expanded reference: organoid branch 2')]:
        matching = [s for s in sources if s['path'].endswith('/' + relative)]
        if len(matching) != 1:
            raise ValueError(f'Missing registered summary figure: {relative}')
        entry = matching[0]
        entry['download'] = f"data/files/{entry['path']}"
        downloads['figures'].append({'label': label, **entry})
    visuals = None
    directory = root / 'site/public/visuals'
    if (directory / 'index.json').exists():
        visual_index = validate_visuals(root, directory)
        visuals = {'index': 'visuals/index.json', 'materialized_at_commit': visual_index['source_commit'],
                   'validated_at_commit': commit, 'validation': 'source-and-display-sha256'}
    preview_path = 'evaluation/cell/streams/preview_manifest.json'
    preview = json.loads((root / preview_path).read_text())
    register_file(preview_path, download=False)
    if preview['schema_version'] != 1 or preview['inputs_sha256'] != hashlib.sha256((root / 'evaluation/cell/inputs.tsv').read_bytes()).hexdigest():
        raise ValueError('Cell preview inputs changed')
    for relative, field in [('evaluation/cell/render_stream_previews.py','renderer_sha256'), ('evaluation/cell/stream_preview_config.json','palette_sha256')]:
        if hashlib.sha256((root / relative).read_bytes()).hexdigest() != preview[field]:
            raise ValueError('Cell preview renderer or palette changed')
        register_file(relative, download=False)
    stream_keys = keyed(preview['entries'], ['dataset','method'])
    expected_keys = {(r['dataset'],r['method']) for r in cell}
    if set(stream_keys) != expected_keys:
        raise ValueError('Cell previews must cover the complete registered input set')
    streams = []
    stream_copies = []
    bindings = {(r['dataset_id'],r['config_id']):r for r in inputs}
    for key, entry in stream_keys.items():
        binding = bindings[key]
        if (entry['result_path'],entry['result_sha256']) != (binding['result_path'],binding['sha256']):
            raise ValueError(f'Cell preview uses a different method output: {key}')
        image = (root / entry['path']).resolve()
        if not image.is_relative_to((root / 'evaluation/cell/streams').resolve()):
            raise ValueError('Cell preview path escapes its source directory')
        if hashlib.sha256(image.read_bytes()).hexdigest() != entry['sha256']:
            raise ValueError(f'Cell preview image changed: {key}')
        target = f"previews/{entry['dataset']}/{entry['method']}.png"
        streams.append({'dataset':entry['dataset'],'method':entry['method'],'path':'data/'+target,'legend':entry['legend']})
        stream_copies.append((image,target))
    overview_path = 'site/public/overview.manifest.json'
    overview = json.loads((root / overview_path).read_text())
    for relative, checksum in [(overview['source'],overview['source_sha256']),('site/public/overview.png',overview['preview_sha256'])]:
        if hashlib.sha256((root / relative).read_bytes()).hexdigest() != checksum:
            raise ValueError('Overview figure changed')
        register_file(relative, download=False)
    register_file(overview_path,download=False)
    descriptions = {key: entry['context'].capitalize() + '.' for key, entry in content['datasets'].items()}
    labels = {
        CELL:('Cell-level method scores','Five metrics for 20 configurations across 14 real datasets.'),
        GENE:('Gene-level directional scores','Per-gene CBDir scores on the fixed reference neighborhood.'),
        COSINE:('Simulated gene-velocity agreement','CosSim values against known simulated ground-truth velocities.'),
        gene_summary_path:('Gene-score dataset summaries','Stored means across evaluated genes with finite CBDir scores.'),
        KINETICS+'gene_metrics.csv':('Main benchmark gene scores','Per-gene kinetic scores used by the canonical benchmark summaries.'),
        KINETICS+'component_means.csv':('Main benchmark parameter summaries','Finite-gene means by system, parameter, and organoid branch.'),
        EXPANDED+'coverage_summary.csv':('Extended-reference coverage and means','Candidate, present, valid, and missing counts with available-gene means.'),
        EXPANDED+'common_gene_summary.csv':('Fixed common-gene comparisons','Means for the frozen common-gene comparison groups.'),
        EXPANDED+'tsvelo_organoid_gene_metrics.csv':('TSvelo organoid transcription scores','Extended-reference alpha scores; undefined constant predictions remain missing.'),
    }
    for item in downloads['tables']:
        item['label'],item['description'] = labels[item['path']]
    for item in downloads['figures']:
        item['description'] = 'Published cross-task summary.' if item['path'].endswith('figSeven.pdf') else 'Distribution of finite kinetic scores across the extended reference-gene set.'
    # Diagnostic scores remain outside the primary Cell/Gene/Kinetics arrays.
    asset_path = "site/content/diagnostic_assets.json"
    asset_binding = json.loads((root / asset_path).read_text())
    register_file(asset_path, download=False)
    if asset_binding["schema_version"] != 1:
        raise ValueError("Unsupported diagnostic asset schema")
    diagnostic_images = {entry["id"]: media_asset(entry["source"], entry["sha256"])
                         for entry in asset_binding["assets"]}
    ae = "supplement/veloae/"
    projection_renderer = ae + "plot_shared_projection.py"
    if register_file(projection_renderer, download=False)["sha256"] != asset_binding["projection_renderer_sha256"]:
        raise ValueError("veloAE projection renderer changed")
    aliases = next(ast.literal_eval(n.value) for n in ast.parse((root / projection_renderer).read_text()).body
                   if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "ALIASES" for t in n.targets))
    projection_manifest = json.loads((root / (ae + "results/manifest.json")).read_text())
    register_file(ae + "results/manifest.json", download=False)
    if [r["dataset_id"] for r in projection_manifest["panels"]] != list(aliases.values()):
        raise ValueError("veloAE panel order differs from its renderer")
    if set(aliases.values()) != {r["dataset"] for r in cell}:
        raise ValueError("veloAE projections must cover the 14 registered Cell datasets")
    source_inputs = json.loads((root / (ae + "inputs/manifest.json")).read_text())
    register_file(ae + "inputs/manifest.json", download=False)
    for entry in source_inputs["files"]:
        if register_file(ae + entry["path"], download=False)["sha256"] != entry["sha256"]:
            raise ValueError("veloAE frozen diagnostic input changed")
    final_manifest = json.loads((root / (ae + "results/ppt_assets_final/manifest.json")).read_text())
    register_file(ae + "results/ppt_assets_final/manifest.json", download=False)
    for name in ["pca_summary.svg", "veloae_pca_velocity.png"]:
        expected_image = next(r["sha256"] for r in final_manifest["files"] if r["name"] == name)
        if register_file(ae + "results/ppt_assets_final/" + name, download=False)["sha256"] != expected_image:
            raise ValueError("veloAE final source image changed")
    pca_rows = read(ae + "inputs/pca_comparison.csv", download=False)
    if len(keyed(pca_rows, ["dataset", "method", "matrix_key"])) != 84:
        raise ValueError("veloAE PCA comparison requires 84 paired source rows")
    ae_panels = []
    for (stem, dataset_id), source in zip(aliases.items(), projection_manifest["panels"]):
        input_record = next(r for r in source_inputs["files"] if r["path"] == f"inputs/projections/{stem}.npz")
        if source["npz_sha256"] != input_record["sha256"]:
            raise ValueError("veloAE projection uses a different frozen input")
        ae_panels.append({"dataset": dataset_id, **diagnostic_images["ae-panel-"+dataset_id]})

    vi = "supplement/velovi/results/"
    vi_manifest = json.loads((root / (vi + "manifest.json")).read_text())
    register_file(vi + "manifest.json", download=False)
    for name, checksum in vi_manifest["outputs"].items():
        if register_file(vi + name, download=False)["sha256"] != checksum:
            raise ValueError("veloVI diagnostic table changed")
    vi_full = read(vi + "full_dataset_three_condition.csv", download=False)
    conditions = ["native_ms_native_velocity", "reference_ms_native_velocity", "reference_ms_affine_velocity"]
    if len(keyed(vi_full, ["dataset", "condition"])) != 42 or {r["condition"] for r in vi_full} != set(conditions):
        raise ValueError("veloVI requires 14 datasets and three distinct conditions")
    vi_table = read(vi + "table_s6.csv", download=False)
    vi_rows = [{"dataset": ensemble_dataset_ids[r["dataset"]],
                **{key: numeric(r[key]) for key in ["A_native_Ms_V", "B_reference_Ms_V", "C_reference_Ms_affine_V", "C_minus_B"]}}
               for r in vi_table if r["dataset"] != "Macro-average"]
    vi_macro = next(r for r in vi_table if r["dataset"] == "Macro-average")
    vi_macro = {key: numeric(vi_macro[key]) for key in ["A_native_Ms_V", "B_reference_Ms_V", "C_reference_Ms_affine_V", "C_minus_B"]}
    vi_full_rows = [{"dataset": ensemble_dataset_ids[r["dataset"]], "condition": r["condition"],
                     "transcosine": numeric(r["transcosine"]), "highdim_cbdir": numeric(r["highdim_cbdir"]),
                     "n_cells": int(r["n_cells"]), "n_common_genes": int(r["n_common_genes"])} for r in vi_full]
    if set(r["dataset"] for r in vi_rows) != set(aliases.values()):
        raise ValueError("veloVI scope must match the registered 14 Cell datasets")

    ti = "supplement/tivelo_gene/"
    assembly_path = "figures/assembly/suppfig13/assembly_manifest.json"
    assembly = json.loads((root / assembly_path).read_text())
    register_file(assembly_path, download=False)
    for entry in assembly["panels"]:
        if register_file(entry["source"], download=False)["sha256"] != entry["sha256"]:
            raise ValueError("TIvelo diagnostic source differs from its assembly")
    simulation_manifest = json.loads((root / (ti+"figures/simulation_source_manifest.json")).read_text())
    register_file(ti+"figures/simulation_source_manifest.json", download=False)
    for source, checksum in simulation_manifest["sources"].items():
        if register_file(source.removeprefix("repository:"), download=False)["sha256"] != checksum:
            raise ValueError("TIvelo simulation root-pair input changed")
    root_pairs = read(ti+"root_stage_intervention/root_stage_summary.csv", download=False)
    selected_pairs = [{"dataset": SIM_IDS[r["scenario"]], "condition": r["variant"], "n_genes": int(r["n_genes"]),
                       "cbdir": numeric(r["mean_cbdir_reference30"]), "cossim": numeric(r["mean_cossim"])}
                      for r in root_pairs if r["scenario"] in ["branching", "all"] and r["variant"] in ["correct_path_root", "oracle_early_root"]]
    if len(keyed(selected_pairs, ["dataset", "condition"])) != 4:
        raise ValueError("TIvelo displayed simulation pairs are incomplete")
    bone_rows = read(ti+"bonemarrow_root_control_20261001/paired_summary.csv", download=False)
    bone_pairs = [{"metric": r["metric"], "n_pairs": int(r["n_pairs"]),
                   **{key: numeric(r[key]) for key in ["original_mean", "hsc1_mean", "delta_mean"]}}
                  for r in bone_rows if r["metric"] in ["cbdir", "batc"]]
    if len(bone_pairs) != 2 or any(r["n_pairs"] != 3 for r in bone_pairs):
        raise ValueError("TIvelo bone-marrow comparison requires three paired optimizer seeds")
    bone_manifest_path = ti+"figures/source_manifest.json"
    bone_manifest = json.loads((root/bone_manifest_path).read_text())
    register_file(bone_manifest_path, download=False)
    if register_file(ti+"bonemarrow_root_control_20261001/paired_summary.csv", download=False)["sha256"] != bone_manifest["source_hashes"]["paired_summary.csv"]:
        raise ValueError("TIvelo bone-marrow summary differs from the displayed source")
    for path in [ae+"README.md", "supplement/velovi/README.md", ti+"README.md"]:
        register_file(path, download=False)
    analysis_downloads = [
        (ae+"inputs/pca_comparison.csv", "veloAE PC1 comparison"),
        (ti+"root_stage_intervention/root_stage_summary.csv", "TIvelo simulation root-pair summary"),
        (ti+"figures/simulation_plotted_pairs.csv", "TIvelo simulation paired gene scores"),
        (ti+"bonemarrow_root_control_20261001/paired_summary.csv", "TIvelo bone-marrow paired summary"),
        (ti+"figures/plotted_seed_pairs.csv", "TIvelo bone-marrow paired seed scores"),
    ]
    downloads["analyses"] = []
    for relative, label in analysis_downloads:
        entry = register_file(relative, download=False)
        entry["download"] = f"data/files/{relative}"
        downloads["analyses"].append({"label": label, **entry})
    diagnostics = {
        "veloae": {"method": "C11", "summary": diagnostic_images["ae-summary"],
                   "plate": diagnostic_images["ae-plate"], "panels": ae_panels,
                   "pca": pca_rows, "csv": "data/files/"+ae+"inputs/pca_comparison.csv"},
        "velovi": {"method": "C12", "rows": vi_rows, "full": vi_full_rows, "macro": vi_macro,
                   "conditions": conditions, "csv": "data/files/"+vi+"table_s6.csv",
                   "full_csv": "data/files/"+vi+"full_dataset_three_condition.csv"},
        "tivelo": {"method": "C08", "simulation_order": diagnostic_images["ti-simulation-order"],
                   "simulation_scores": diagnostic_images["ti-simulation-scores"],
                   "bone_image": diagnostic_images["ti-bonemarrow"], "simulation": selected_pairs,
                   "bonemarrow": bone_pairs, "displayed_seed": bone_manifest["displayed_seed"],
                   "simulation_csv": "data/files/"+ti+"figures/simulation_plotted_pairs.csv",
                   "bone_csv": "data/files/"+ti+"figures/plotted_seed_pairs.csv"},
    }
    payload = {"schema_version": 1, "source_commit": commit, "build_state": build_state, "methods": methods, "datasets": datasets,
               "content": content, "ensemble": ensemble, "diagnostics": diagnostics, "dataset_descriptions": descriptions, "streams": streams, "cell": cell, "gene": gene, "gene_summary": gene_summary, "visuals": visuals, 'downloads': downloads,
               "kinetics": kinetics, "kinetics_summary": summary,
               "kinetics_common_membership": common_membership, "kinetics_method_order": kinetics_order,
               "kinetics_colors": kinetics_colors,
               "kinetics_expanded_scope": {"sceu_rpe1": len(universe["rpe1"]), "sceu_organoid": len(universe["organoid"])},
               "artifacts": artifacts, 'h5ad': h5ad, "sources": sources}
    for source in sources:
        if hashlib.sha256((root / source["path"]).read_bytes()).hexdigest() != source["sha256"]:
            raise ValueError(f"Source changed during export: {source['path']}")
    output.mkdir(parents=True, exist_ok=True)
    for source in sources:
        if source["download"]:
            destination = output / "files" / source["path"]
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(root / source["path"], destination)
    for source, target in stream_copies + media_copies:
        destination = output / target
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source,destination)
    shutil.copyfile(root / "site/content/README.md", output / "content-guide.md")
    (output / "benchmark.json").write_text(json.dumps(payload, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    (output / "manifest.json").write_text(json.dumps({"schema_version": 1, "source_commit": commit, "build_state": build_state, "sources": sources, "cell_preview_manifest": preview_path, "visual_bundle": visuals}, indent=2), encoding="utf-8")
    return payload


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Generated static data directory")
    args = parser.parse_args()
    exported = export_data(output=args.output)
    print(f"Exported {len(exported['cell'])} cell rows, {len(exported['gene'])} gene rows and "
          f"{len(exported['kinetics'])} kinetic rows; source {exported['source_commit'][:7]}")
