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
        self.assertTrue(all(e['projection']=='fixed_reference30_pca30_projection_pc1_pc2' for e in forebrain))
        self.assertTrue(all(e['frame']['projection_dimensions']==30 and
                            e['frame']['display_components']==[1,2] for e in forebrain))

    def test_previous_figure7_export_is_hidden_while_author_export_is_pending(self):
        manifest = json.loads((ROOT / 'figures/fig7_composite_manifest.json').read_text())
        if manifest.get('pending_pdf_export'):
            self.assertNotIn(manifest['output']['path'], {r['path'] for r in self.payload['downloads']['figures']})

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
        self.assertEqual(len(tables), 17)
        figure_manifest = json.loads((ROOT / 'figures/fig7_composite_manifest.json').read_text())
        self.assertEqual(len(figures), 3 if figure_manifest.get('pending_pdf_export') else 4)
        self.assertEqual(len(metadata), 11)
        audit = self.payload['downloads']['audit']
        analyses = self.payload['downloads']['analyses']
        self.assertEqual(len(analyses), 5)
        exposed = {r['path'] for r in tables + figures + metadata + audit + analyses}
        self.assertEqual(exposed, {s['path'] for s in self.payload['sources'] if s['download']})
        self.assertTrue(all(r['path'].endswith('.csv') for r in tables))
        self.assertTrue(all(r['path'].endswith(('.png', '.pdf')) for r in figures))
        self.assertEqual({r['path'] for r in metadata}, {
            'data/RESULT_DOWNLOADS.md',
            'data/result_package_manifest.tsv',
            'data/result_member_manifest.tsv',
            'data/input_release_manifest.tsv',
            'data/input_file_manifest.tsv',
            'data/upstream_references.tsv',
            'data/UPSTREAM_REFERENCES.md',
            'data/kinetics_reference_gene_usage/README.md',
            'data/kinetics_reference_gene_usage/manifest.json',
            'data/kinetics_reference_gene_usage/rpe1_reference_gene_usage.csv',
            'data/kinetics_reference_gene_usage/organoid_reference_gene_usage.csv',
        })
        self.assertEqual(len(audit), 13)
        self.assertEqual(len(self.payload['ensemble']['scores']), 28)
        self.assertEqual(self.payload['content']['schema_version'], 1)
        self.assertEqual(len(self.payload['content']['datasets']), 17)
        self.assertFalse(any('/trajectories/' in p or p.endswith(('.svg', '.zip')) for p in exposed))
        with (ROOT / 'data/release_map.tsv').open() as handle:
            h5ad_ids = {r['artifact_id'] for r in csv.DictReader(handle, delimiter='\t') if r['member_path'].endswith('.h5ad')}
        self.assertEqual({r['artifact_id'] for r in self.payload['h5ad']}, h5ad_ids)

    def test_individual_inputs_and_upstream_sources_are_separate(self):
        files = self.payload['input_files']
        self.assertEqual(len(files), 25)
        self.assertEqual(len({f['artifact_id'] for f in files}), 25)
        self.assertEqual({ext: sum(f['format'] == ext for f in files) for ext in ['H5AD', 'LOOM', 'CSV']},
                         {'H5AD': 22, 'LOOM': 1, 'CSV': 2})
        references = self.payload['upstream_references']
        self.assertEqual({r['table'] for r in references}, {'S1', 'S3', 'S4'})
        self.assertFalse({r['artifact_id'] for r in references} & {f['artifact_id'] for f in files})
        self.assertTrue(all(f['url'].endswith(f['path'] + '?download=true') for f in files))
        self.assertFalse(any('.tar.gz' in f['url'] for f in files))

    def test_result_packages_preserve_archive_and_member_identities(self):
        packages = self.payload['result_packages']
        members = self.payload['result_members']
        self.assertEqual(len(packages), 339)
        self.assertEqual(len(members), 388)
        self.assertEqual(sum(row['bytes'] for row in packages), 149576197521)
        self.assertEqual(sum(row['member_bytes'] for row in members), 252167762795)
        self.assertEqual(sum(row['format'] == 'tar.gz' for row in packages), 325)
        self.assertEqual(sum(row['format'] == 'npz' for row in packages), 14)
        self.assertEqual({row['release_version'] for row in packages}, {'v1.0.0'})
        base = [row for row in packages if '/diagnostics/multivelo/' not in row['repository_path']]
        self.assertEqual(len(base), 338)
        self.assertEqual(sum(row['bytes'] for row in base), 149558222960)
        supplement = [row for row in packages if '/diagnostics/multivelo/' in row['repository_path']]
        self.assertEqual(len(supplement), 1)
        self.assertEqual(supplement[0]['sha256'], '34ea57fb9f11c7b7dfa60bae78ebccfe2f2eaf680fe3c8fe9b92abfe3b962fcf')
        self.assertIn('/resolve/v1.0.0/v1.0.0/diagnostics/multivelo/', supplement[0]['url'])
        indexed = {row['package_id']: row for row in packages}
        for package in packages:
            selected = [row for row in members if row['package_id'] == package['package_id']]
            self.assertEqual(len(selected), package['scientific_files'])
            self.assertEqual(sum(row['member_bytes'] for row in selected), package['source_bytes'])
        registry = {row['artifact_id']: row for row in self.payload['artifacts']}
        for member in members:
            self.assertIn(member['package_id'], indexed)
            self.assertNotIn('..', Path(member['archive_member']).parts)
            for identifier in member['registered_artifact_ids'].split(';') if member['registered_artifact_ids'] else []:
                self.assertEqual(registry[identifier]['sha256'], member['member_sha256'])

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

    def test_content_scope_matches_registry_and_frozen_method_coverage(self):
        content = self.payload['content']
        with (ROOT / 'data/datasets.tsv').open() as handle:
            datasets = {r['dataset_id']: r for r in csv.DictReader(handle, delimiter='\t')}
        with (ROOT / 'methods/configurations.tsv').open() as handle:
            methods = {r['config_id'] for r in csv.DictReader(handle, delimiter='\t')}
        self.assertEqual(set(content['datasets']), set(datasets))
        self.assertEqual(set(content['methods']), methods)
        scenarios = json.loads((ROOT / 'evaluation/gene/simulation/config.json').read_text())['scenarios']
        for source_name, dataset_id in EXPORT.SIM_IDS.items():
            scenario = content['simulations'][dataset_id]
            generator = scenarios['hybrid' if source_name == 'all' else source_name]
            self.assertEqual([scenario['cells'], scenario['genes']], generator['expected_shape'])
            self.assertEqual(scenario['normal_genes'] + scenario['murk_genes'] + scenario['branching_genes'], scenario['genes'])
        scored_gene_methods = {r['method'] for r in self.payload['gene']}
        marked_components = {key for key, entry in content['methods'].items() if entry['gene_components'] == 'spliced + unspliced velocity components'}
        self.assertEqual(marked_components, scored_gene_methods)
        self.assertEqual(len(scored_gene_methods), 9)
        self.assertEqual(content['methods']['C13']['kinetic']['parameters'], ['beta', 'gamma'])
        self.assertEqual(content['methods']['C14']['kinetic']['parameters'], ['alpha'])
        self.assertNotIn('dynamo_labelled', methods)
        self.assertNotIn('C07', {r['method'] for r in self.payload['kinetics']})

    def test_ensemble_projection_preserves_scores_ranks_and_macro_averages(self):
        ensemble = self.payload['ensemble']
        aliases = {'capofibro': 'capo_pre', 'dentate_lamanno': 'dentategyrus_lamanno',
                   'neuron': 'scnt_neuron', 'organoid': 'sceu_organoid', 'rpe1': 'sceu_rpe1'}
        actual = EXPORT.keyed(ensemble['scores'], ['dataset', 'variant'])
        with (ROOT / 'ensemble/results/scores.csv').open() as handle:
            for row in csv.DictReader(handle):
                key = (aliases.get(row['dataset'], row['dataset']), row['variant'])
                for metric in ['cbdir', 'batc']:
                    self.assertEqual(actual[key][metric], float(row[metric]))
        self.assertEqual(len(actual), 28)
        ranks = EXPORT.keyed(ensemble['rankings'], ['dataset', 'variant', 'metric'])
        with (ROOT / 'ensemble/results/rankings_vs_20.csv').open() as handle:
            for row in csv.DictReader(handle):
                key = (aliases.get(row['dataset'], row['dataset']), row['variant'], row['metric'])
                self.assertEqual(ranks[key]['score'], float(row['score']))
                self.assertEqual(ranks[key]['rank'], int(row['rank']))
                self.assertEqual(ranks[key]['out_of'], 21)
                self.assertEqual(ranks[key]['tie_count'], int(row['tie_count']))
        self.assertEqual(len(ranks), 56)
        metadata = json.loads((ROOT / 'ensemble/results/metadata.json').read_text())
        self.assertEqual(ensemble['macro_average'], metadata['macro_average'])
        for variant in ensemble['variants']:
            rows = [r for r in ensemble['scores'] if r['variant'] == variant]
            for metric in ['cbdir', 'batc']:
                self.assertAlmostEqual(sum(r[metric] for r in rows) / len(rows), ensemble['macro_average'][variant][metric], places=14)
        self.assertTrue(all(r['method'] in {m['config_id'] for m in self.payload['methods']} for r in self.payload['cell']))

    def test_ensemble_images_cover_and_preserve_all_frozen_results(self):
        scores = EXPORT.keyed(self.payload['ensemble']['scores'], ['dataset', 'variant'])
        images = EXPORT.keyed(self.payload['ensemble']['streams'], ['dataset', 'variant'])
        self.assertEqual(set(images), set(scores))
        self.assertEqual(len(images), 28)
        with (ROOT / 'artifacts.tsv').open() as handle:
            registered = {r['path']: r['sha256'] for r in csv.DictReader(handle, delimiter='\t') if r['category'] == 'ensemble-stream'}
        for image in images.values():
            self.assertEqual(image['sha256'], registered[image['source']])
            original = (ROOT / image['source']).read_bytes()
            self.assertEqual(original, (self.output / image['path'].removeprefix('data/')).read_bytes())

    def test_veloae_panels_follow_registered_dataset_and_frozen_input_scope(self):
        ae = self.payload['diagnostics']['veloae']
        panel_manifest = json.loads((ROOT / 'supplement/veloae/results/manifest.json').read_text())
        self.assertEqual([r['dataset'] for r in ae['panels']], [r['dataset_id'] for r in panel_manifest['panels']])
        self.assertEqual(len(ae['panels']), 14)
        with (ROOT / 'supplement/veloae/inputs/pca_comparison.csv').open() as handle:
            self.assertEqual(ae['pca'], list(csv.DictReader(handle)))
        self.assertEqual(len(ae['pca']), 84)
        self.assertTrue(ae['summary']['source'].endswith('pca_summary.svg'))
        for image in ae['panels'] + [ae['summary'], ae['plate']]:
            self.assertEqual((ROOT / image['source']).read_bytes(), (self.output / image['path'].removeprefix('data/')).read_bytes())

    def test_velovi_diagnostic_conditions_and_archived_macro_are_preserved(self):
        vi = self.payload['diagnostics']['velovi']
        self.assertEqual(len(vi['full']), 42)
        self.assertEqual(len(vi['rows']), 14)
        config = json.loads((ROOT / 'evaluation/cell/config.json').read_text())['datasets']
        ids = {key: row['dataset_id'] for key, row in config.items()}
        actual = EXPORT.keyed(vi['full'], ['dataset', 'condition'])
        with (ROOT / 'supplement/velovi/results/full_dataset_three_condition.csv').open() as handle:
            for row in csv.DictReader(handle):
                exported = actual[(ids[row['dataset']], row['condition'])]
                for metric in ['transcosine', 'highdim_cbdir']:
                    self.assertEqual(exported[metric], float(row[metric]))
        with (ROOT / 'supplement/velovi/results/table_s6.csv').open() as handle:
            macro = next(r for r in csv.DictReader(handle) if r['dataset'] == 'Macro-average')
        self.assertEqual(vi['macro'], {k: float(v) for k, v in macro.items() if k != 'dataset'})
        self.assertEqual(set(vi['conditions']), {'native_ms_native_velocity','reference_ms_native_velocity','reference_ms_affine_velocity'})

    def test_tivelo_diagnostic_refits_keep_scope_and_source_values(self):
        ti = self.payload['diagnostics']['tivelo']
        self.assertEqual(len(ti['simulation']), 4)
        with (ROOT / 'supplement/tivelo_gene/root_stage_intervention/root_stage_summary.csv').open() as handle:
            source = {(EXPORT.SIM_IDS[r['scenario']], r['variant']): r for r in csv.DictReader(handle)}
        for row in ti['simulation']:
            original = source[(row['dataset'], row['condition'])]
            self.assertEqual(row['n_genes'], 20 if row['dataset'] == 'sim_only_branching' else 40)
            self.assertEqual(row['cbdir'], float(original['mean_cbdir_reference30']))
            self.assertEqual(row['cossim'], float(original['mean_cossim']))
        with (ROOT / 'supplement/tivelo_gene/bonemarrow_root_control_20261001/paired_summary.csv').open() as handle:
            source = {r['metric']: r for r in csv.DictReader(handle)}
        for row in ti['bonemarrow']:
            self.assertEqual(row['n_pairs'], 3)
            for field in ['original_mean', 'hsc1_mean', 'delta_mean']:
                self.assertEqual(row[field], float(source[row['metric']][field]))
        for field in ['simulation_order','simulation_scores','bone_image']:
            image = ti[field]
            self.assertEqual((ROOT / image['source']).read_bytes(), (self.output / image['path'].removeprefix('data/')).read_bytes())
        self.assertEqual(len(self.payload['cell']), 280)

    def test_visual_bundle_generation_and_validation_commits_are_distinct(self):
        historical = json.loads((ROOT / 'site/public/visuals/index.json').read_text())['source_commit']
        visual = self.payload['visuals']
        self.assertEqual(visual['materialized_at_commit'], historical)
        self.assertEqual(visual['validated_at_commit'], self.payload['source_commit'])
        manifest = json.loads((self.output / 'manifest.json').read_text())
        self.assertEqual(manifest['visual_bundle'], visual)

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


    def test_cell_overall_statistics_are_unshifted_committed_values(self):
        path = ROOT / EXPORT.SUMMARY / 'fig7_cell_metric_summary.csv'
        actual = EXPORT.keyed(self.payload['cell_summary'], ['method', 'metric'])
        self.assertEqual(len(actual), 100)
        with path.open() as handle:
            for source in csv.DictReader(handle):
                row = actual[(source['config_id'], source['metric'])]
                self.assertEqual(row['n_datasets'], 14)
                for key in ['mean', 'median', 'sample_sd', 'q1', 'q3', 'iqr']:
                    self.assertEqual(row[key], float(source[key]))
                self.assertNotIn('mean_plus_1', row)
        summary = self.payload['overall_summary']
        self.assertEqual(len(summary), 34)
        self.assertTrue(any(r['method']=='dynamo_labelled' for r in summary))
        self.assertFalse(any(r['method']=='C07' for r in summary))

    def test_simulation_cell_stays_separate_and_preserves_frozen_scores(self):
        self.assertEqual(len(self.payload['cell']), 280)
        self.assertFalse(any(r['dataset'].startswith('sim_') for r in self.payload['cell']))
        actual = EXPORT.keyed(self.payload['cell_simulation'], ['dataset','method'])
        self.assertEqual(len(actual), 27)
        with (ROOT / EXPORT.SIM_CELL / 'inputs.tsv').open() as handle:
            bindings = EXPORT.keyed(list(csv.DictReader(handle, delimiter='\t')), ['evaluation_dataset','operational_id'])
        with (ROOT / EXPORT.SIM_CELL / 'results/cell_metrics.csv').open() as handle:
            for source in csv.DictReader(handle):
                binding = bindings[(source['dataset'],source['method'])]
                row = actual[(binding['dataset_id'],binding['config_id'])]
                for key in ['consistency','iccoh','cbdir','transcosine','batc']:
                    self.assertEqual(row[key], float(source[key]))
        self.assertEqual(len(self.payload['streams']), 280)

    def test_cross_level_rows_preserve_pairing_gene_sd_and_hybrid_scope(self):
        actual = EXPORT.keyed(self.payload['cross_level'], ['case','method'])
        self.assertEqual(len(actual), 36)
        categories = {r['category']:r['id'] for r in self.payload['cross_cases']}
        with (ROOT / EXPORT.SIM_CELL / 'inputs.tsv').open() as handle:
            bindings = EXPORT.keyed(list(csv.DictReader(handle, delimiter='\t')), ['evaluation_dataset','operational_id'])
        with (ROOT / EXPORT.SIM_CELL / 'results/cross_level_summary.csv').open() as handle:
            for source in csv.DictReader(handle):
                method = bindings[(source['dataset'],source['operational_id'])]['config_id']
                row = actual[(categories[source['category']],method)]
                for key in ['cell_cbdir','cell_batc','gene_cbdir_mean','gene_cbdir_sd','gene_cossim_mean','gene_cossim_sd']:
                    self.assertEqual(row[key],float(source[key]))
                self.assertEqual(row['n_genes'],int(source['n_genes']))
        for method in {r['method'] for r in self.payload['cross_level']}:
            for key in ['cell_cbdir','cell_batc']:
                self.assertEqual(actual[('hybrid-murk',method)][key],actual[('hybrid-branching',method)][key])

    def test_real_cases_keep_registered_edge_images_and_gene_scope(self):
        cases = self.payload['real_cases']
        self.assertEqual(len(cases), 3)
        self.assertEqual({(c['dataset'],c['gene'],tuple(c['transition'])) for c in cases}, {
            ('capo_pre','FAM111A',('G1','S')),
            ('dentategyrus_lamanno','Ncald',('ImmGranule1','ImmGranule2')),
            ('gastrulation','Fech',('Erythroid2','Erythroid3')),
        })
        with (ROOT / 'figures/assembly/assembly_manifest.tsv').open() as handle:
            assembly = {r['source_path']:r for r in csv.DictReader(handle,delimiter='\t') if r['figure']=='fig3'}
        for case in cases:
            image = case['image']
            self.assertEqual(image['sha256'],assembly[image['source']]['sha256'])
            self.assertEqual((self.output / image['path'].removeprefix('data/')).read_bytes(),(ROOT / image['source']).read_bytes())
        self.assertEqual(len(self.payload['cross_level']),36)
        self.assertEqual(len(self.payload['cross_cases']),4)

    def test_kinetic_fields_use_their_own_images_and_labelled_identity(self):
        rows = self.payload['kinetics_streams']
        self.assertEqual(len(rows),10)
        self.assertEqual({(r['dataset'],r['method']) for r in rows},{
            (dataset,method) for dataset in ['sceu_rpe1','sceu_organoid']
            for method in ['C13','C14','C15','C16','dynamo_labelled']})
        with (ROOT / 'figures/assembly/assembly_manifest.tsv').open() as handle:
            assembly = {r['source_path']:r for r in csv.DictReader(handle,delimiter='\t') if r['figure']=='fig6' and r['panel'] in ['A','E']}
        for row in rows:
            image = row['image']
            self.assertEqual(image['sha256'],assembly[image['source']]['sha256'])
            self.assertEqual((self.output / image['path'].removeprefix('data/')).read_bytes(),(ROOT / image['source']).read_bytes())
            self.assertNotIn('cbdir',row)
        self.assertFalse(any(r['method']=='C07' for r in rows))

    def test_strategy_labels_and_macro_f1_preserve_committed_class_sets(self):
        actual = self.payload['kinetics_strategy']
        self.assertEqual(actual['n_genes'],250)
        self.assertEqual(len(actual['genes']),250)
        genes = {r['gene']:r for r in actual['genes']}
        for method,stem,name in [('C15','celldancer','cellDancer'),('C16','deepvelo','DeepVelo')]:
            with (ROOT / EXPORT.PLOT / ('rpe1_'+stem+'_kinetic_class_membership.csv')).open() as handle:
                members = list(csv.DictReader(handle))
            reference = {r['Gene']:r['Class'] for r in members if r['scEU']=='True'}
            predicted = {r['Gene']:r['Class'] for r in members if r[name]=='True'}
            self.assertEqual(set(reference),set(genes))
            for gene in genes:
                self.assertEqual(genes[gene]['reference'],reference[gene])
                self.assertEqual(genes[gene]['predictions'][method],predicted[gene])
            with (ROOT / EXPORT.PLOT / ('rpe1_'+stem+'_kinetic_classes.csv')).open() as handle:
                counts = list(csv.DictReader(handle))
            expected = sum(2*int(r['intersection'])/(2*int(r['intersection'])+int(r['scEU_only'])+int(r[name+'_only'])) for r in counts)/3
            projected = next(r for r in actual['methods'] if r['method']==method)
            self.assertAlmostEqual(projected['macro_f1'],expected,places=14)
            self.assertEqual({r['class'] for r in projected['counts']},set(actual['classes']))

if __name__ == "__main__":
    unittest.main()
