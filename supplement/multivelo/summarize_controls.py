"""Rebuild full-precision diagnostic summaries using only the standard library."""
from pathlib import Path
import argparse
import csv
import statistics

ROOT = Path(__file__).resolve().parent
SUBSETS = [("branching", "branching", 20), ("all", "branching", 20), ("murk", "MURK", 10), ("all", "MURK", 20)]
CONDITIONS = [
    ("scVelo_dyn", "Original", "saved"),
    ("scVelo_dyn", "Induction only", "force_induction_fit_tau_smoothed"),
    ("MultiVelo (RNA_only)", "Default", "refit_default"),
    ("MultiVelo (RNA_only)", "Gate open", "on_allow_repression"),
    ("MultiVelo (RNA_only)", "No output smoothing", "refit_no_output_smoothing"),
    ("MultiVelo (RNA_only)", "Complete trajectory", "refit_complete"),
]


def read(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def write(path, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(result_root=ROOT / "results"):
    scv = read(result_root / "scvelo_controls.csv")
    mv = read(result_root / "multivelo_refits.csv")
    summaries = []
    for dataset, category, n in SUBSETS:
        for method, condition, source in CONDITIONS:
            if source == "on_allow_repression" and category == "MURK":
                continue
            frame = scv if method == "scVelo_dyn" else mv
            key = "variant" if method == "scVelo_dyn" else "condition"
            rows = [r for r in frame if r["dataset"] == dataset and r[key] == source and ("MURK" if r["gene"].startswith("OneSidedBoost") else "branching") == category]
            if len(rows) != n or len({r["gene"] for r in rows}) != n:
                raise ValueError(f"Unexpected gene coverage: {dataset}, {category}, {source}")
            row = dict(dataset=dataset, gene_class=category, method=method, condition=condition, n_genes=n)
            for metric in ["cbdir", "cossim"]:
                values = [float(r[metric]) for r in rows]
                row[f"mean_{metric}"] = statistics.mean(values)
                row[f"sd_{metric}"] = statistics.stdev(values)
            summaries.append(row)
    effects = []
    for dataset in ["branching", "all"]:
        base = {r["gene"]: r for r in mv if r["dataset"] == dataset and r["condition"] == "refit_default" and not r["gene"].startswith("OneSidedBoost")}
        changed = {r["gene"]: r for r in mv if r["dataset"] == dataset and r["condition"] == "on_allow_repression"}
        assert base.keys() == changed.keys() and len(base) == 20
        row = dict(dataset=dataset, n_gene_pairs=20)
        for metric in ["cbdir", "cossim", "fit_loss"]:
            row[f"n_lower_{metric}"] = sum(float(changed[g][metric]) < float(base[g][metric]) for g in base)
        effects.append(row)
    return summaries, effects


def manuscript_table(module_root=ROOT):
    """Select the four branching controls for S9; retain the full result CSVs."""
    rows = read(module_root / "results/condition_summary.csv")
    selected = {(r["dataset"], r["method"], r["condition"]): r for r in rows if r["gene_class"] == "branching"}
    latex = [r"\begin{table}[p]", r"\centering", r"\footnotesize", r"\setlength{\tabcolsep}{3pt}", r"\renewcommand{\arraystretch}{1.15}", r"\begin{tabular}{@{}llrrrr@{}}", r"\toprule",
             r" & & \multicolumn{2}{c}{\shortstack{normal + branching\\branching ($n=20$)}} & \multicolumn{2}{c}{\shortstack{hybrid\\branching ($n=20$)}} \\",
             r"\cmidrule(lr){3-4}\cmidrule(lr){5-6}", r"Method & Condition & $\mathrm{CBDir}_{\mathrm{gene}}$ & CosSim & $\mathrm{CBDir}_{\mathrm{gene}}$ & CosSim \\", r"\midrule"]
    for method, condition, _ in CONDITIONS[:4]:
        values = []
        for dataset in ["branching", "all"]:
            row = selected[dataset, method, condition]
            if int(row["n_genes"]) != 20:
                raise ValueError("S9 requires 20 evaluated branching genes in each subset")
            values.extend(f"{float(row[key]):.4f}" for key in ["mean_cbdir", "mean_cossim"])
        label = "MultiVelo (RNA only)" if method.startswith("MultiVelo") else method
        latex.append(label.replace("_", r"\_") + " & " + condition + " & " + " & ".join(values) + r" \\")
    caption = (module_root / "CAPTION.md").read_text().split("## Supplementary Table S9\n\n", 1)[1].strip()
    title, body = caption.split("\n\n", 1)
    latex += [r"\bottomrule", r"\end{tabular}", r"\caption{\textbf{\normalsize " + title.replace("scVelo_dyn", r"scVelo\_dyn") + r"}\label{tab:suppMultiveloPhaseControls}", r"\footnotesize\normalfont", body + "}", r"\end{table}"]
    return "\n".join(latex) + "\n"


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--result-root", type=Path, default=ROOT / "results")
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--table-only", action="store_true", help="Export compact manuscript S9 without regenerating scientific CSVs")
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    if args.table_only:
        (args.output_root / "tableS9.tex").write_text(manuscript_table(args.result_root.parent))
        print("Wrote four-row manuscript Table S9")
        return
    summaries, effects = summarize(args.result_root)
    write(args.output_root / "condition_summary.csv", summaries)
    write(args.output_root / "candidate_gate_effects.csv", effects)
    print("Wrote", len(summaries), "condition summaries and", len(effects), "paired-effect summaries")


if __name__ == "__main__":
    main()
