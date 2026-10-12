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

    def register_file(relative, download=True):
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
    artifacts = [{key: row[key] for key in ["artifact_id", "category", "sha256", "status", "uri"]}
                 for row in read("artifacts.tsv", download=False) if row["status"] == "external"]
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
    commit = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    build_state = "working-tree-preview" if subprocess.check_output(["git","-C",str(root),"status","--porcelain"],text=True).strip() else "committed"
    downloads = {'tables': [], 'figures': [], 'metadata': []}
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
        validate_visuals(root, directory)
        visuals = {'index': 'visuals/index.json'}
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
    descriptions = {
        'pancreas':'Mouse pancreatic endocrinogenesis.',
        'gastrulation':'Mouse gastrulation erythroid lineage.',
        'forebrain':'Developing human forebrain glutamatergic lineage.',
        'dentategyrus_lamanno':'Mouse dentate gyrus neurogenesis.',
        'dentategyrus':'Mouse dentate gyrus neurogenesis.',
        'retina':'Mouse retinal development.',
        'mousebrain':'Embryonic E18 mouse brain.',
        'hematopoiesis':'Human hematopoiesis.',
        'bonemarrow':'Human bone marrow CD34+ differentiation.',
        'sceu_rpe1':'RPE1 cell cycle.',
        'sceu_organoid':'Mouse intestinal organoid differentiation.',
        'scnt_neuron':'Mouse cortical neuron activation.',
        'capo_pre':'Human dermal fibroblasts.',
        'rpe1_gen':'RPE1 cell cycle.',
        'sim_only_murk':'Normal + MURK simulation.',
        'sim_only_branching':'Normal + branching simulation.',
        'sim_hybrid':'Hybrid normal + MURK + branching simulation.',
    }
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
    payload = {"schema_version": 1, "source_commit": commit, "build_state": build_state, "methods": methods, "datasets": datasets,
               "dataset_descriptions": descriptions, "streams": streams, "cell": cell, "gene": gene, "gene_summary": gene_summary, "visuals": visuals, 'downloads': downloads,
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
    for source, target in stream_copies:
        destination = output / target
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source,destination)
    (output / "benchmark.json").write_text(json.dumps(payload, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    (output / "manifest.json").write_text(json.dumps({"schema_version": 1, "source_commit": commit, "build_state": build_state, "sources": sources, "cell_preview_manifest": preview_path}, indent=2), encoding="utf-8")
    return payload


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Generated static data directory")
    args = parser.parse_args()
    exported = export_data(output=args.output)
    print(f"Exported {len(exported['cell'])} cell rows, {len(exported['gene'])} gene rows and "
          f"{len(exported['kinetics'])} kinetic rows; source {exported['source_commit'][:7]}")
