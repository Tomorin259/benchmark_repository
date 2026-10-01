"""Assemble the frozen eight simulation panels with one progression colorbar."""
import argparse
import copy
import hashlib
import io
import json
from pathlib import Path

from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import RectangleObject
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

PT = 72
WIDTH, HEIGHT = 7.5 * PT, 4.28 * PT
LEFT, PANEL, GAP = .23, 1.70, .08
ROW_TOPS = [.27, 2.02]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--materials-root', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--font-file', type=Path, required=True,
                        help='Arial TTF used by the approved original assembly')
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pdfmetrics.registerFont(TTFont('Arial', str(args.font_file)))
    writer = PdfWriter()
    page = writer.add_blank_page(WIDTH, HEIGHT)
    sources = []

    def place(filename, x, top, width, height, clip=None):
        path = args.materials_root / filename
        source = copy.deepcopy(PdfReader(path).pages[0])
        box = source.mediabox if clip is None else RectangleObject(clip)
        if clip is not None:
            source.cropbox = RectangleObject(clip)
        scale = min(width / float(box.width), height / float(box.height))
        tx = x + (width - float(box.width) * scale) / 2 - float(box.left) * scale
        ty = HEIGHT - top - height + (height - float(box.height) * scale) / 2 - float(box.bottom) * scale
        page.merge_transformed_page(source, Transformation().scale(scale).translate(tx, ty))
        sources.append({'source': filename, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                        'target_points': [x, top, width, height], 'clip_points': clip})

    for row, scenario in enumerate(['branching', 'hybrid']):
        for col, kind in enumerate(['states', 'known', 'main_dpt', 'branch_dpt']):
            place(f'ab_assets/individual/{scenario}_{kind}.pdf',
                  (LEFT + col * (PANEL + GAP)) * PT,
                  ROW_TOPS[row] * PT, PANEL * PT, PANEL * PT)
    place('ab_assets/order_comparison_material.pdf',
          WIDTH - .18 * PT - 193, 3.91 * PT, 193, 21, [149, 13, 342, 34])
    stream = io.BytesIO()
    c = canvas.Canvas(stream, pagesize=(WIDTH, HEIGHT))
    c.setFont('Arial', 8)
    for col, title in enumerate(['Root selection', 'Known order', 'Main DPT', 'Branch DPT']):
        c.drawCentredString((LEFT + col * (PANEL + GAP) + PANEL / 2) * PT,
                           HEIGHT - .18 * PT, title)
    c.save()
    page.merge_page(PdfReader(stream).pages[0])
    writer.add_metadata({'/Title': 'TIvelo root and ordering materials',
                         '/Subject': 'Frozen panels; one normalized progression scale'})
    with (args.output_dir / 'ab_right_legend_material.pdf').open('wb') as handle:
        writer.write(handle)
    (args.output_dir / 'right_legend_source_manifest.json').write_text(json.dumps({
        'canvas_inches': [7.5, 4.28], 'sources': sources,
        'font_sha256': hashlib.sha256(args.font_file.read_bytes()).hexdigest(),
        'data_changes': False, 'panel_letters': False,
        'normalization': 'Known developmental bins / 9; DPT unchanged, normalized per path',
        'colorbar': 'One unlabeled 0--1 scale at bottom right'
    }, indent=2) + '\n')


if __name__ == '__main__':
    main()
