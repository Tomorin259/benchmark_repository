"""Export the two native assets for manual A4 layout in PowerPoint."""
import argparse
import hashlib
import json
import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

from plot_shared_projection import ORDER, COLUMNS, ROWS, PANEL_WIDTH_MM, PANEL_HEIGHT_MM, FOOTER_MM


def prefix_svg_ids(root, prefix):
    for element in root.iter():
        for key, value in list(element.attrib.items()):
            if key == 'id':
                element.set(key, prefix+value)
            elif 'url(#' in value:
                element.set(key, value.replace('url(#', 'url(#'+prefix))
            elif key.endswith('href') and value.startswith('#'):
                element.set(key, '#'+prefix+value[1:])


def write_manifest(output):
    files = sorted(path for path in output.iterdir()
                   if path.is_file() and path.name != 'manifest.json')
    manifest = {'purpose': 'two native assets for manual A4 PowerPoint layout',
                'assets': ['pca_summary', 'veloae_pca_velocity'],
                'projection_layout': {'columns': COLUMNS, 'rows': ROWS},
                'files': [{'name': path.name,
                           'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                          for path in files]}
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--refresh-manifest', action='store_true',
                        help='Only refresh checksums after external SVG-to-PDF conversion')
    args = parser.parse_args()
    source, output = args.input_dir, args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    if args.refresh_manifest:
        write_manifest(output)
        return
    for suffix in ('png', 'pdf', 'svg'):
        shutil.copy2(source/f'pca_summary.{suffix}', output/f'pca_summary.{suffix}')

    ns = 'http://www.w3.org/2000/svg'
    ET.register_namespace('', ns)
    width_mm = COLUMNS*PANEL_WIDTH_MM
    height_mm = ROWS*PANEL_HEIGHT_MM+FOOTER_MM
    root = ET.Element(f'{{{ns}}}svg', width=f'{width_mm}mm', height=f'{height_mm}mm',
                      viewBox=f'0 0 {width_mm} {height_mm}')
    ET.SubElement(root, f'{{{ns}}}rect', width=str(width_mm), height=str(height_mm), fill='white')
    placements = []
    for i, stem in enumerate(ORDER):
        placements.append((source/f'panel_{i+1:02d}_{stem}.svg',
                           (i % COLUMNS)*PANEL_WIDTH_MM, (i//COLUMNS)*PANEL_HEIGHT_MM,
                           PANEL_WIDTH_MM, PANEL_HEIGHT_MM))
    placements.extend([(source/'shared_key.svg', (COLUMNS-1)*PANEL_WIDTH_MM,
                        (ROWS-1)*PANEL_HEIGHT_MM, PANEL_WIDTH_MM, PANEL_HEIGHT_MM),
                       (source/'orientation.svg', 0, ROWS*PANEL_HEIGHT_MM, PANEL_WIDTH_MM, FOOTER_MM)])
    for path, x, y, width, height in placements:
        panel = ET.parse(path).getroot()
        prefix_svg_ids(panel, path.stem+'_')
        panel.set('x', str(x)); panel.set('y', str(y))
        panel.set('width', str(width)); panel.set('height', str(height))
        root.append(panel)
    projection_svg = output/'veloae_pca_velocity.svg'
    ET.ElementTree(root).write(projection_svg, encoding='utf-8', xml_declaration=True)
    shutil.copy2(source/'veloae_pca_velocity_5x3.png',
                 output/'veloae_pca_velocity.png')
    if shutil.which('rsvg-convert'):
        subprocess.run(['rsvg-convert', '-f', 'pdf', '-o',
                        str(output/'veloae_pca_velocity.pdf'), str(projection_svg)], check=True)

    write_manifest(output)


if __name__ == '__main__':
    main()
