"""Assemble the accepted S15 A/B materials without rescaling source panels.

The author retains control of the final PPTX layout and PDF export. This entry
point reuses the vector-merge and raster-placement operations from that assembly.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil

from PIL import Image
from pypdf import PdfReader, PdfWriter, Transformation
from pypdf._page import PageObject

REPO = Path(__file__).resolve().parents[3]
MODULE = REPO / "supplement/multivelo"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise ValueError("Use an empty output directory")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    material_manifest = json.loads((MODULE / "materials_manifest.json").read_text())
    placements = [
        ("phase_truth", 0, 0),
        ("phase_multivelo_default", 174.6, 0),
        ("phase_multivelo_gate_open", 349.2, 0),
        ("progression_key", 431.262, 177.12),
    ]
    inputs = {}
    for stem in [row[0] for row in placements] + ["paired_controls"]:
        for extension in ["pdf", "png"]:
            relative = f"figures/{stem}.{extension}"
            source = MODULE / relative
            assert sha(source) == material_manifest["outputs"][relative]
            inputs[str(source.relative_to(REPO))] = sha(source)

    width, height, dpi = 518.4, 203.04, 600
    page = PageObject.create_blank_page(width=width, height=height)
    canvas = Image.new("RGB", (round(width * dpi / 72), round(height * dpi / 72)), "white")
    for stem, x, y in placements:
        source = PdfReader(MODULE / "figures" / f"{stem}.pdf").pages[0]
        page.merge_transformed_page(source, Transformation().translate(x, height - y - float(source.mediabox.height)))
        with Image.open(MODULE / "figures" / f"{stem}.png") as image:
            assert image.width == round(float(source.mediabox.width) * dpi / 72)
            assert image.height == round(float(source.mediabox.height) * dpi / 72)
            canvas.paste(image.convert("RGB"), (round(x * dpi / 72), round(y * dpi / 72)))
    writer = PdfWriter()
    writer.add_page(page)
    writer.write(args.output_dir / "A_phase_portraits.pdf")
    canvas.save(args.output_dir / "A_phase_portraits.png", dpi=(dpi, dpi))
    for extension in ["pdf", "png"]:
        shutil.copyfile(MODULE / "figures" / f"paired_controls.{extension}", args.output_dir / f"B_paired_controls.{extension}")
    record = dict(
        figure="S15", inputs=inputs, placements_A=placements,
        A_size_points=[width, height], panel_scale=1.0,
        transform="translation only; native vector panels and original raster samples retained",
        panel_letters=False, primary_results_modified=False,
        script_sha256=sha(Path(__file__)),
        outputs={p.name: sha(p) for p in args.output_dir.iterdir() if p.is_file()},
        final_export="Author-manual PPTX layout and PDF export; see suppfig15_composite_manifest.json",
    )
    (args.output_dir / "assembly_manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    print("Assembled S15 A/B materials:", args.output_dir)


if __name__ == "__main__":
    main()
