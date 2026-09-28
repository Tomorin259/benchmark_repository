"""Assemble A/B/C; retain native SVG panels and C's 300-dpi source detail."""
import argparse
import hashlib
import json
import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
from matplotlib import font_manager
from plot_shared_projection import ORDER


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, required=True)
    root = parser.parse_args().input_dir
    ns = 'http://www.w3.org/2000/svg'
    ET.register_namespace('', ns)
    width, top_height, c_top, height = 180, 180*420/1060, 78, 259
    svg = ET.Element(f'{{{ns}}}svg', width=f'{width}mm', height=f'{height}mm',
                     viewBox=f'0 0 {width} {height}')
    ET.SubElement(svg, f'{{{ns}}}rect', width='180', height=str(height), fill='white')
    placements = [(root/'pca_summary.svg', 0, 3, width, top_height)]
    for i, stem in enumerate(ORDER):
        placements.append((root/f'panel_{i+1:02d}_{stem}.svg',
                           (i % 4)*45, c_top+(i//4)*43, 45, 43))
    placements.extend([(root/'shared_key.svg', 135, c_top+129, 45, 43),
                       (root/'orientation.svg', 0, c_top+172, 45, 9)])
    for path, x, y, w, h in placements:
        panel = ET.parse(path).getroot()
        prefix = path.stem+'_'
        for element in panel.iter():
            for key, value in list(element.attrib.items()):
                if key == 'id':
                    element.set(key, prefix+value)
                elif 'url(#' in value:
                    element.set(key, value.replace('url(#', 'url(#'+prefix))
                elif key.endswith('href') and value.startswith('#'):
                    element.set(key, '#'+prefix+value[1:])
        panel.set('x', str(x)); panel.set('y', str(y))
        panel.set('width', str(w)); panel.set('height', str(h))
        svg.append(panel)
    labels = [('A', .8, 5), ('B', 90.8, 5), ('C', .8, 77)]
    for label, x, y in labels:
        text = ET.SubElement(svg, f'{{{ns}}}text', x=str(x), y=str(y),
                             style='font-family:DejaVu Sans,sans-serif;font-size:3.2px;font-weight:bold')
        text.text = label
    svg_path = root/'veloae_supplement.svg'
    ET.ElementTree(svg).write(svg_path, encoding='utf-8', xml_declaration=True)
    manifest = {'size_mm': [width, height], 'png_dpi': 600,
                'scatter_source_dpi': 300, 'sources': [
                    {'name': p.name, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest(),
                     'placement_mm': [x, y, w, h]} for p, x, y, w, h in placements]}
    (root/'assembly_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    pixels_per_mm = 600/25.4
    canvas = Image.new('RGB', (round(width*pixels_per_mm), round(height*pixels_per_mm)), 'white')
    with Image.open(root/'pca_summary.png') as top:
        canvas.paste(top.convert('RGB'), (0, round(3*pixels_per_mm)))
    with Image.open(root/'veloae_pca_velocity_4x4.png') as plate:
        target = (round(180*pixels_per_mm), round(181*pixels_per_mm))
        canvas.paste(plate.convert('RGB').resize(target, Image.Resampling.LANCZOS),
                     (0, round(c_top*pixels_per_mm)))
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.truetype(font_manager.findfont('DejaVu Sans'), round(3.2*pixels_per_mm))
    for label, x, y in labels:
        draw.text((round(x*pixels_per_mm), round((y-3.2)*pixels_per_mm)), label, font=font, fill='#262626')
    canvas.save(root/'veloae_supplement.png', dpi=(600, 600))
    if shutil.which('rsvg-convert'):
        subprocess.run(['rsvg-convert', '-f', 'pdf', '-o', str(root/'veloae_supplement.pdf'),
                        str(svg_path)], check=True)


if __name__ == '__main__':
    main()
