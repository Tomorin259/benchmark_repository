"""Verify the website projection against its scientific source tables."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import shutil

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("site_export", ROOT / "site/export_data.py")
EXPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EXPORT)


class SiteExportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.output = Path(cls.temp.name)
        cls.payload = EXPORT.export_data(ROOT, cls.output)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_every_cell_value_matches_frozen_input_identity(self):
        with (ROOT / EXPORT.CELL).open() as handle:
            raw = {(r["dataset"], r["method"]): r for r in csv.DictReader(handle)}
        with (ROOT / "evaluation/cell/inputs.tsv").open() as handle:
            inputs = list(csv.DictReader(handle, delimiter="\t"))
        actual = EXPORT.keyed(self.payload["cell"], ["dataset", "method"])
        self.assertEqual(len(actual), len(inputs))
        for binding in inputs:
            exported = actual[(binding["dataset_id"], binding["config_id"])]
            source = raw[(binding["evaluation_dataset"], binding["operational_id"])]
            for metric in ["consistency", "iccoh", "cbdir", "transcosine", "batc"]:
                self.assertEqual(exported[metric], float(source[metric]))

    def test_gene_missing_values_and_simulation_cosine_are_preserved(self):
        actual = EXPORT.keyed(self.payload["gene"], ["dataset", "method", "gene"])
        with (ROOT / EXPORT.GENE).open() as handle:
            for source in csv.DictReader(handle):
                ds = EXPORT.SIM_IDS.get(source["dataset"], source["dataset"].lower())
                row = actual[(ds, EXPORT.GENE_IDS[source["method"]], source["gene"])]
                self.assertEqual(row["cbdir"], EXPORT.numeric(source["new_cbdir"]))
        self.assertIsNone(actual[("retina", "C18", "Anks1b")]["cbdir"])
        with (ROOT / EXPORT.COSINE).open() as handle:
            for source in csv.DictReader(handle):
                row = actual[(EXPORT.SIM_IDS[source["Dataset"]], EXPORT.GENE_IDS[source["Method"]], source["Gene"])]
                self.assertEqual(row["cossim"], float(source["CosSim"]))

    def test_selected_gene_scopes_and_labelled_control_stay_distinct(self):
        selected = [r for r in self.payload["kinetics"] if r["scope"] == "selected"]
        for dataset, branch, n_genes in [("sceu_rpe1", "", 19), ("sceu_organoid", "1", 6), ("sceu_organoid", "2", 8)]:
            genes = {r["gene"] for r in selected if r["dataset"] == dataset and r["branch"] == branch}
            self.assertEqual(len(genes), n_genes)
        self.assertTrue(any(r["method"] == "dynamo_labelled" for r in self.payload["kinetics"]))
        self.assertFalse(any(r["method"] == "C07" for r in self.payload["kinetics"]))
        with (ROOT / (EXPORT.KINETICS + "component_means.csv")).open() as handle:
            source = list(csv.DictReader(handle))
        for original, exported in zip(source, self.payload["kinetics_summary"]):
            self.assertEqual(exported["score"], float(original["Mean"]))
            self.assertEqual(exported["n_genes"], int(original["N"]))
            self.assertEqual(exported["branch"], original["Branch"])

    def test_original_download_bytes_and_source_hashes_are_unchanged(self):
        for source in self.payload["sources"]:
            original = (ROOT / source["path"]).read_bytes()
            self.assertEqual(hashlib.sha256(original).hexdigest(), source["sha256"])
            if source["download"]:
                self.assertEqual(original, (self.output / source["download"].removeprefix("data/")).read_bytes())
        serialized = (self.output / "benchmark.json").read_text()
        json.loads(serialized, parse_constant=lambda value: self.fail(f"Non-JSON number: {value}"))
        for private_prefix in ["/Users/", "/mnt/", "DB_PASSWORD", "JWT_SECRET"]:
            self.assertNotIn(private_prefix, serialized)
        with (ROOT / "artifacts.tsv").open() as handle:
            external = [r for r in csv.DictReader(handle, delimiter="\t") if r["status"] == "external"]
        self.assertEqual(len(self.payload["artifacts"]), len(external))
        self.assertEqual([r["uri"] for r in self.payload["artifacts"]], [r["uri"] for r in external])

    def test_all_kinetic_values_preserve_their_own_source_scope(self):
        actual = EXPORT.keyed(self.payload["kinetics"], ["scope", "dataset", "branch", "parameter", "method", "gene"])
        method_ids = {"DeepVelo": "C16", "DeepKINET": "C13", "TSvelo": "C14",
                      "cellDancer": "C15", "Dynamo(label data)": "dynamo_labelled", "Dynamo": "dynamo_labelled"}
        with (ROOT / (EXPORT.KINETICS + "gene_metrics.csv")).open() as handle:
            for source in csv.DictReader(handle):
                key = ("full", "sceu_" + source["Dataset"], source["Branch"], source["Parameter"], method_ids[source["Method"]], source["Gene"])
                self.assertEqual(actual[key]["score"], EXPORT.numeric(source["Spearman_Corr"]))
                self.assertEqual(actual[key]["n_samples"], int(source["Valid_Samples"]))
        for parameter in ["alpha", "gamma"]:
            with (ROOT / (EXPORT.PLOT + f"rpe1_{parameter}_heatmap.csv")).open() as handle:
                for source in csv.DictReader(handle):
                    for gene, score in source.items():
                        if gene != "Method":
                            key = ("selected", "sceu_rpe1", "", parameter, method_ids[source["Method"]], gene)
                            self.assertEqual(actual[key]["score"], EXPORT.numeric(score))
        for branch in ["1", "2"]:
            with (ROOT / (EXPORT.PLOT + f"organoid_branch{branch}_gene_metrics.csv")).open() as handle:
                for source in csv.DictReader(handle):
                    key = ("selected", "sceu_organoid", branch, source["Parameter"].lower(),
                           method_ids[source["Method"].removesuffix(f"_branch{branch}")], source["Gene"])
                    self.assertEqual(actual[key]["score"], EXPORT.numeric(source["Correlation"]))

    def test_duplicate_and_scientific_output_paths_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            EXPORT.keyed([{"id": "a"}, {"id": "a"}], ["id"])
        with self.assertRaises(ValueError):
            EXPORT.export_data(ROOT, ROOT / "evaluation/cell/results")

    def test_expanded_scores_are_separate_and_tsvelo_keeps_undefined_values(self):
        expanded = [r for r in self.payload["kinetics"] if r["scope"] == "expanded"]
        canonical = [r for r in self.payload["kinetics"] if r["scope"] == "full"]
        self.assertEqual(len(expanded), 4748)
        self.assertEqual(sum(r["score"] is not None for r in expanded), 4567)
        self.assertEqual(len(canonical), 4736)
        self.assertFalse(any(r["method"] == "C14" and r["dataset"] == "sceu_organoid" for r in canonical))
        self.assertEqual(self.payload["kinetics_expanded_scope"], {"sceu_rpe1": 250, "sceu_organoid": 270})
        tsvelo = [r for r in expanded if r["method"] == "C14" and r["dataset"] == "sceu_organoid"]
        for branch in ["1", "2"]:
            rows = [r for r in tsvelo if r["branch"] == branch]
            self.assertEqual(len(rows), 37)
            self.assertEqual(sum(r["score"] is not None for r in rows), 20)
            self.assertEqual(sum(r["status"] == "constant_prediction" for r in rows), 17)
            self.assertTrue(all(r["parameter"] == "alpha" for r in rows))
        actual = EXPORT.keyed(expanded, ["dataset", "branch", "parameter", "method", "gene"])
        with (ROOT / (EXPORT.EXPANDED + "tsvelo_organoid_gene_metrics.csv")).open() as handle:
            for row in csv.DictReader(handle):
                record = actual[("sceu_organoid", row["Branch"], row["Parameter"], "C14", row["Gene"])]
                self.assertEqual(record["score"], EXPORT.numeric(row["Spearman_Corr"]))
                self.assertEqual(record["n_samples"], int(row["Valid_Samples"]))

    def test_expanded_coverage_and_common_views_match_committed_tables(self):
        ids = {"DeepVelo": "C16", "DeepKINET": "C13", "TSvelo": "C14", "cellDancer": "C15", "Dynamo(label data)": "dynamo_labelled"}
        actual = EXPORT.keyed([r for r in self.payload["kinetics_summary"] if r["scope"] == "expanded"],
                              ["dataset", "branch", "parameter", "comparison", "method"])
        with (ROOT / (EXPORT.EXPANDED + "coverage_summary.csv")).open() as handle:
            coverage = list(csv.DictReader(handle))
        panels = {}
        for row in coverage:
            branch = str(int(float(row["Branch"]))) if row["Branch"] else ""
            key = ("sceu_" + row["Dataset"], branch, row["Parameter"])
            panels[row["Panel"]] = key
            record = actual[key + ("available", ids[row["Method"]])]
            self.assertEqual(record["score"], EXPORT.numeric(row["Mean"]))
            for field, source in [("candidate_n", "Candidate_N"), ("present_n", "Present_N"),
                                  ("n_genes", "Valid_N"), ("undefined_n", "Undefined_N"), ("unavailable_n", "Unavailable_N")]:
                self.assertEqual(record[field], int(row[source]))
        with (ROOT / (EXPORT.EXPANDED + "common_gene_summary.csv")).open() as handle:
            common = list(csv.DictReader(handle))
        self.assertEqual(len(actual), len(coverage) + len(common))
        for row in common:
            record = actual[panels[row["Panel"]] + (row["Comparison_scope"], ids[row["Method"]])]
            self.assertEqual(record["score"], EXPORT.numeric(row["Mean"]))
            self.assertEqual(record["n_genes"], int(row["Common_genes_N"]))
        members = {}
        with (ROOT / (EXPORT.EXPANDED + "common_gene_membership.csv")).open() as handle:
            for row in csv.DictReader(handle):
                members.setdefault(panels[row["Panel"]] + (row["Comparison_scope"],), []).append(row["Gene"])
        exported = {(r["dataset"], r["branch"], r["parameter"], r["comparison"]): r["genes"]
                    for r in self.payload["kinetics_common_membership"]}
        self.assertEqual(exported, members)

    def test_cell_previews_match_the_exact_registered_method_outputs(self):
        with (ROOT / 'evaluation/cell/inputs.tsv').open() as handle:
            bindings = {(r['dataset_id'],r['config_id']):r for r in csv.DictReader(handle,delimiter='\t')}
        index = json.loads((ROOT / 'evaluation/cell/streams/preview_manifest.json').read_text())
        registered = EXPORT.keyed(index['entries'],['dataset','method'])
        displayed = EXPORT.keyed(self.payload['streams'],['dataset','method'])
        self.assertEqual(set(bindings),set(registered))
        self.assertEqual(set(bindings),set(displayed))
        self.assertEqual(len(displayed),280)
        self.assertEqual(sum(int(r['n_common'])<int(r['n_reference']) for r in bindings.values()),4)
        for key, entry in registered.items():
            self.assertEqual(entry['result_path'],bindings[key]['result_path'])
            self.assertEqual(entry['result_sha256'],bindings[key]['sha256'])
            self.assertEqual((ROOT/entry['path']).read_bytes(),(self.output/displayed[key]['path'].removeprefix('data/')).read_bytes())
            self.assertEqual(entry['n_common'],int(bindings[key]['n_common']))
        forebrain = [e for e in registered.values() if e['dataset']=='forebrain']
        self.assertEqual(len(forebrain),20)
        self.assertTrue(all(e['projection']=='fixed_reference30_projection_display_only' for e in forebrain))

    def test_changed_canonical_source_blocks_expanded_export(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = Path(directory) / "repository"
            for source in self.payload["sources"]:
                destination = fixture / source["path"]
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(ROOT / source["path"], destination)
            canonical = fixture / (EXPORT.KINETICS + "gene_metrics.csv")
            canonical.write_bytes(canonical.read_bytes() + b"\n")
            output = Path(directory) / "output"
            with self.assertRaisesRegex(ValueError, "Canonical kinetics scores"):
                EXPORT.export_data(fixture, output)
            self.assertFalse(output.exists())

    def test_gene_summary_uses_committed_finite_gene_means(self):
        actual = EXPORT.keyed(self.payload['gene_summary'], ['dataset', 'method'])
        path = ROOT / 'evaluation/gene/results/gene_cbdir_summary_by_dataset_method.csv'
        with path.open() as handle:
            for source in csv.DictReader(handle):
                dataset = EXPORT.SIM_IDS.get(source['dataset'], source['dataset'].lower())
                row = actual[(dataset, EXPORT.GENE_IDS[source['method']])]
                self.assertEqual(row['cbdir'], float(source['new_mean']))
                self.assertEqual(row['n_genes'], int(source['n_genes']))

    def test_downloads_only_expose_curated_tables_and_summary_figures(self):
        tables = self.payload['downloads']['tables']
        figures = self.payload['downloads']['figures']
        metadata = self.payload['downloads']['metadata']
        self.assertEqual(len(tables), 9)
        self.assertEqual(len(figures), 4)
        self.assertEqual(len(metadata), 5)
        exposed = {r['path'] for r in tables + figures + metadata}
        self.assertEqual(exposed, {s['path'] for s in self.payload['sources'] if s['download']})
        self.assertTrue(all(r['path'].endswith('.csv') for r in tables))
        self.assertTrue(all(r['path'].endswith(('.png', '.pdf')) for r in figures))
        self.assertEqual({r['path'] for r in metadata}, {
            'data/input_release_manifest.tsv',
            'data/kinetics_reference_gene_usage/README.md',
            'data/kinetics_reference_gene_usage/manifest.json',
            'data/kinetics_reference_gene_usage/rpe1_reference_gene_usage.csv',
            'data/kinetics_reference_gene_usage/organoid_reference_gene_usage.csv',
        })
        self.assertFalse(any('/trajectories/' in p or p.endswith(('.svg', '.zip')) for p in exposed))
        with (ROOT / 'data/release_map.tsv').open() as handle:
            h5ad_ids = {r['artifact_id'] for r in csv.DictReader(handle, delimiter='\t') if r['member_path'].endswith('.h5ad')}
        self.assertEqual({r['artifact_id'] for r in self.payload['h5ad']}, h5ad_ids)

    def test_visual_source_or_asset_drift_blocks_publication(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            visuals = root / 'site/public/visuals'
            visuals.mkdir(parents=True)
            (root / 'source.csv').write_text('frozen\n')
            (visuals / 'curve.json').write_text('{"y":[1,null,3]}')
            index = {'schema_version': 1,
                     'sources': [{'path': 'source.csv', 'sha256': hashlib.sha256((root / 'source.csv').read_bytes()).hexdigest()}],
                     'files': [{'path': 'curve.json', 'sha256': hashlib.sha256((visuals / 'curve.json').read_bytes()).hexdigest()}]}
            (visuals / 'index.json').write_text(json.dumps(index))
            EXPORT.validate_visuals(root, visuals)
            (root / 'source.csv').write_text('changed\n')
            with self.assertRaisesRegex(ValueError, 'Visual source changed'):
                EXPORT.validate_visuals(root, visuals)
            (root / 'source.csv').write_text('frozen\n')
            (visuals / 'curve.json').write_text('{"y":[1,0,3]}')
            with self.assertRaisesRegex(ValueError, 'Visual asset changed'):
                EXPORT.validate_visuals(root, visuals)

    @unittest.skipUnless((ROOT / 'site/public/visuals/index.json').exists(), 'Static visual bundle not materialized')
    def test_materialized_visuals_match_paper_display_data_and_scopes(self):
        directory = ROOT / 'site/public/visuals'
        index = EXPORT.validate_visuals(ROOT, directory)
        phase_keys = {(p['dataset'], p['method'], p['gene']) for p in index['phases']}
        self.assertEqual(phase_keys, {(r['dataset'], r['method'], r['gene']) for r in self.payload['gene']})
        self.assertEqual(len({p['gene'] for p in index['phases'] if p['dataset'] == 'sim_hybrid'}), 40)
        curves = EXPORT.keyed(index['curves'], ['dataset', 'branch', 'parameter', 'gene'])
        methods = {'scEU-seq': 'reference', 'DeepVelo': 'C16', 'DeepKINET': 'C13',
                   'TSvelo': 'C14', 'cellDancer': 'C15', 'Dynamo(label data)': 'dynamo_labelled'}
        for parameter in ['alpha', 'gamma']:
            loaded = {}
            with (ROOT / (EXPORT.PLOT + f'rpe1_{parameter}_trajectories.csv')).open() as handle:
                for row in csv.DictReader(handle):
                    if row['Gene'] not in loaded:
                        item = curves[('sceu_rpe1', '', parameter, row['Gene'])]
                        document = json.loads((ROOT / 'site/public' / item['path']).read_text())
                        loaded[row['Gene']] = {s['method']: s for s in document['series']}
                    series = loaded[row['Gene']][methods[row['Source']]]
                    position = float(row['Relative_Position'])
                    sample = round(position * (len(series['x']) - 1))
                    self.assertAlmostEqual(series['x'][sample], position, places=14)
                    expected = EXPORT.numeric(row['Display_Z'])
                    if expected is None:
                        self.assertIsNone(series['y'][sample])
                    else:
                        self.assertAlmostEqual(series['y'][sample], expected, places=12)
        for branch in ['1', '2']:
            tsvelo = []
            for key, item in curves.items():
                if key[:3] != ('sceu_organoid', branch, 'alpha'):
                    continue
                document = json.loads((ROOT / 'site/public' / item['path']).read_text())
                self.assertEqual(document['transform'], 'Saved profiles: prediction baseline ratio; reference Table S4 values')
                self.assertTrue(all(s['scale'] == ('reference' if s['method'] == 'reference' else 'prediction') for s in document['series']))
                tsvelo.extend(s for s in document['series'] if s['method'] == 'C14')
            self.assertEqual(len(tsvelo), 37)
            self.assertTrue(all(len(s['y']) == 200 for s in tsvelo))
            self.assertEqual(sum(len({v for v in s['y'] if v is not None}) == 1 for s in tsvelo), 17)


if __name__ == "__main__":
    unittest.main()
