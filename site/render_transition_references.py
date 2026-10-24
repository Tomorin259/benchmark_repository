#!/usr/bin/env python3
"""Draw the annotated evaluation transitions on the saved reference frames."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import runpy
import struct
import textwrap

ROOT = Path(__file__).resolve().parents[1]
SOURCES = (
    "data/datasets.tsv", "evaluation/cell/inputs.tsv",
    "evaluation/cell/config.json", "evaluation/cell/stream_preview_config.json",
    "evaluation/cell/simulation/config.json",
    "evaluation/cell/simulation/reference_manifest.json", "evaluation/gene/config.py",
    "evaluation/gene/plot_gene_assets.py", "evaluation/cell/reference_neighbors/forebrain.npz",
)
LABEL_OFFSETS = {
    "dentategyrus_lamanno": {"Nbl1": (-12, -16), "Nbl2": (-16, 18),
        "nIPC": (16, -3), "GlialProg": (-20, -24), "ImmAstro": (-5, 24), "OPC": (10, -16)},
    "forebrain": {"4": (-14, 0), "5": (-14, 0)},
    "retina": {"Neuroblast": (-12, 15), "PR": (0, 24), "AC/HC": (0, 24),
               "RGC": (0, 24), "Progenitor": (0, -24)},
    "sim_hybrid": {"BranchB_0": (0, 24), "BranchB_1": (0, 24),
                   "BranchB_5": (20, 0), "BranchA_8": (0, -24)},
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def reference_definitions(root):
    """Join the registered reference files to each task's actual edge set."""
    root = Path(root)
    cell = json.loads((root / "evaluation/cell/config.json").read_text())
    simulation = json.loads((root / "evaluation/cell/simulation/config.json").read_text())
    gene = runpy.run_path(str(root / "evaluation/gene/config.py"))
    palette = json.loads((root / "evaluation/cell/stream_preview_config.json").read_text())
    with (root / "data/datasets.tsv").open() as handle:
        datasets = {row["dataset_id"]: row["paper_name"] for row in csv.DictReader(handle, delimiter="\t")}
    with (root / "evaluation/cell/inputs.tsv").open() as handle:
        inputs = {row["evaluation_dataset"]: row for row in csv.DictReader(handle, delimiter="\t")
                  if row["operational_id"] == "scvelo_deterministic"}
    definitions = {}
    for name, settings in cell["datasets"].items():
        row = inputs[name]
        definitions[name] = {
            "dataset": row["dataset_id"], "name": datasets[row["dataset_id"]],
            "source": row["result_path"], "source_sha256": row["sha256"],
            "label": settings["label"], "label_config": cell["labels"],
            "evaluation_dataset": name, "frame": settings["embedding"],
            "palette": palette[name],
            "tasks": {"cell": [[str(a), str(b)] for a, b in settings["cluster_edges"]]},
        }
    for settings in gene["REAL_CONFIGS"].values():
        definitions[Path(settings["neighbor_asset"]).stem]["tasks"]["gene"] = [
            [str(a), str(b)] for a, b in settings["edges"]]
    references = json.loads((root / "evaluation/cell/simulation/reference_manifest.json").read_text())["scenarios"]
    for name, settings in simulation["datasets"].items():
        source = references[name]
        definitions[name] = {
            "dataset": settings["dataset_id"], "name": datasets[settings["dataset_id"]],
            "source": source["source"], "source_sha256": source["source_sha256"],
            "label": settings["label"], "label_config": simulation["labels"],
            "evaluation_dataset": name, "frame": "umap", "palette": None,
            "tasks": {"cell": settings["cluster_edges"],
                      "gene": [list(edge) for edge in gene["SIM_CONFIGS"][name]["edges"]]},
        }
    return list(definitions.values())


def read_reference(root, artifact_root, definition):
    import h5py
    import numpy as np

    path = Path(artifact_root) / definition["source"].split(":", 1)[1]
    if sha256(path) != definition["source_sha256"]:
        raise ValueError(f"Reference identity changed: {definition['dataset']}")
    with h5py.File(path, "r") as handle:
        obs = handle["obs"]
        names = obs[obs.attrs["_index"]][:].astype(str)
        label_config = definition["label_config"]
        override = label_config["dataset_specific_sources"].get(definition["evaluation_dataset"])
        key = override if override and override in obs else definition["label"]
        if key not in obs:
            key = next(candidate for candidate in label_config["fallback_sources"] if candidate in obs)
        column = obs[key]
        if isinstance(column, h5py.Group):
            codes = column["codes"][:]
            if (codes < 0).any():
                raise ValueError(f"Missing reference labels: {definition['dataset']}")
            labels = column["categories"][:].astype(str)[codes]
        else:
            labels = column[:].astype(str)
        xy = handle["obsm/X_umap"][:, :2]
    if definition["frame"] == "pca":
        with np.load(Path(root) / "evaluation/cell/reference_neighbors/forebrain.npz", allow_pickle=False) as saved:
            if not np.array_equal(names, saved["obs_names"].astype(str)):
                raise ValueError("Forebrain reference cell orders differ")
            xy = saved["reference_pca"][:, :2]
    if xy.shape != (len(names), 2) or not np.isfinite(xy).all():
        raise ValueError(f"Invalid reference coordinates: {definition['dataset']}")
    return xy, labels, names, key


def draw_reference(xy, labels, edges, palette, stem, dataset):
    """Place PAGA-style nodes at group medians; draw only configured arrows."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch
    import numpy as np

    if palette:
        colors = dict(zip(palette["categories"], palette["palette"]))
    else:
        progression = [plt.get_cmap("plasma")(value) for value in np.linspace(.06, .94, 10)]
        colors = {label: progression[int(label.rsplit("_", 1)[1])] for label in np.unique(labels)}
    nodes = list(dict.fromkeys(label for edge in edges for label in edge))
    if set(nodes) - set(labels) or set(labels) - set(colors):
        raise ValueError(f"Transition labels or colors are missing: {stem.name}")
    centers = {label: np.median(xy[labels == label], axis=0) for label in nodes}
    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    ax.set_rasterization_zorder(2)
    for label, color in colors.items():
        points = xy[labels == label]
        ax.scatter(points[:, 0], points[:, 1], s=7, alpha=.32,
                   color=color, linewidths=0, zorder=1)
    cyclic = len(nodes) == 3 and all(sum(edge[0] == node for edge in edges) == 1 for node in nodes)
    for source, target in edges:
        curve = .15 if cyclic else .08 if source.startswith("BranchA_") else -.08 if source.startswith("BranchB_") else .035
        arrow = FancyArrowPatch(centers[source], centers[target], arrowstyle="-|>",
                     mutation_scale=18, linewidth=2, color="#27323c", shrinkA=9, shrinkB=11,
                     connectionstyle=f"arc3,rad={curve}", zorder=3)
        ax.add_patch(arrow)
    for node in nodes:
        ax.scatter(*centers[node], s=80, color=colors[node], edgecolors="white", linewidths=1.3, zorder=4)
    spread = np.ptp(xy, axis=0)
    ax.set_xlim(xy[:, 0].min() - spread[0] * .12, xy[:, 0].max() + spread[0] * .12)
    ax.set_ylim(xy[:, 1].min() - spread[1] * .14, xy[:, 1].max() + spread[1] * .14)
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()
    fig.subplots_adjust(left=.015, right=.985, bottom=.015, top=.985)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    boxes = []
    samples = np.vstack([np.linspace(centers[a], centers[b], 25) for a, b in edges])
    samples = ax.transData.transform(samples)
    offsets = [(12, 15), (-12, 15), (12, -15), (-12, -15),
               (0, 24), (0, -24), (28, 0), (-28, 0), (24, 26), (-24, 26), (24, -26), (-24, -26)]
    for node in sorted(nodes, key=lambda n: -sum(n in edge for edge in edges)):
        display = node.replace("BranchA_", "A").replace("BranchB_", "B") if not palette else textwrap.fill(node, 16, break_long_words=False)
        candidates = []
        for offset in offsets:
            annotation = ax.annotate(display, centers[node], xytext=offset, textcoords="offset points",
                ha="left" if offset[0] > 0 else "right" if offset[0] < 0 else "center",
                va="center", fontsize=13, color="#27323c", zorder=5,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": .86, "pad": 1.2})
            box = annotation.get_window_extent(renderer).expanded(1.06, 1.12)
            overlap = sum(max(0, min(box.x1, b.x1) - max(box.x0, b.x0)) *
                          max(0, min(box.y1, b.y1) - max(box.y0, b.y0)) for b in boxes)
            outside = max(0, -box.x0) + max(0, box.x1 - fig.bbox.width) + max(0, -box.y0) + max(0, box.y1 - fig.bbox.height)
            cost = 100 * overlap + 1000 * outside + 15 * box.count_contains(samples) + abs(offset[0]) + abs(offset[1])
            candidates.append((cost, offset, box))
            annotation.remove()
        _, offset, box = min(candidates, key=lambda item: item[0])
        if node in LABEL_OFFSETS.get(dataset, {}):
            offset = LABEL_OFFSETS[dataset][node]
        elif not palette and any(label.startswith("BranchA_") for label in nodes):
            stage = int(node.rsplit("_", 1)[1])
            offset = ((-24, -18) if node.startswith("BranchA_") else (24, -18)) if stage == 0 else \
                     ((-24, 0) if node.startswith("BranchA_") else (24, 0)) if stage == 1 else \
                     (-14, 10) if node.startswith("BranchA_") else (0, 24)
        ax.annotate(display, centers[node], xytext=offset, textcoords="offset points",
            ha="left" if offset[0] > 0 else "right" if offset[0] < 0 else "center",
            va="center", fontsize=13, color="#27323c", zorder=5,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": .86, "pad": 1.2},
            arrowprops={"arrowstyle": "-", "color": "#87939b", "linewidth": .6,
                        "shrinkA": 2, "shrinkB": 7})
        boxes.append(box)
    for extension in ("png", "svg"):
        metadata = {"Date": None} if extension == "svg" else {}
        fig.savefig(stem.with_suffix("." + extension), dpi=160, facecolor="white", metadata=metadata,
                    bbox_inches="tight", pad_inches=.06)
    plt.close(fig)
    return {node: center.tolist() for node, center in centers.items()}


def render_references(root, artifact_root, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    root, output = Path(root), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none", "svg.hashsalt": "evaluation-transitions"})
    entries = []
    for definition in reference_definitions(root):
        xy, labels, names, label_key = read_reference(root, artifact_root, definition)
        groups = {}
        for task, edges in definition["tasks"].items():
            key = tuple(sorted(map(tuple, edges)))
            groups.setdefault(key, {"tasks": [], "edges": edges})["tasks"].append(task)
        for group in groups.values():
            suffix = "" if len(groups) == 1 else "_" + group["tasks"][0]
            stem = output / (definition["dataset"] + suffix)
            centers = draw_reference(xy, labels, group["edges"], definition["palette"], stem, definition["dataset"])
            files = {ext: {"path": "transitions/" + stem.with_suffix("." + ext).name,
                          "sha256": sha256(stem.with_suffix("." + ext)),
                          "bytes": stem.with_suffix("." + ext).stat().st_size} for ext in ("png", "svg")}
            with stem.with_suffix(".png").open("rb") as image:
                width, height = struct.unpack(">II", image.read(24)[16:24])
            files["png"].update(width=width, height=height)
            entries.append({"dataset": definition["dataset"], "tasks": group["tasks"],
                "name": definition["name"], "frame": definition["frame"], "edges": group["edges"],
                "nodes": centers, "files": files,
                "reference": {"source": definition["source"], "sha256": definition["source_sha256"],
                    "label_key": label_key, "n_cells": len(names),
                    "coordinate_sha256": hashlib.sha256(np.ascontiguousarray(xy).tobytes()).hexdigest(),
                    "cell_order_sha256": hashlib.sha256("\n".join(names).encode()).hexdigest()}})
            print(f"Rendered {definition['dataset']}: {' / '.join(group['tasks'])}, {len(group['edges'])} edges", flush=True)
    sources = [{"path": rel, "sha256": sha256(root / rel)} for rel in SOURCES]
    sources.append({"path": "site/render_transition_references.py", "sha256": sha256(__file__)})
    manifest = {"schema_version": 1, "node_position": "group-wise coordinate median",
                "sources": sources, "entries": entries}
    (output / "index.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True, help="The local velocity artifact directory")
    parser.add_argument("--output", type=Path, required=True, help="Directory for PNG, SVG and index.json")
    parser.add_argument("--repository-root", type=Path, default=ROOT)
    args = parser.parse_args()
    render_references(args.repository_root, args.artifact_root, args.output)
