import json
import sys
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SIM = ROOT / "evaluation/gene/simulation"
sys.path.insert(0, str(ROOT))
from evaluation.gene.simulation.generate import _bin_labels, _one_sided, _truth_layers


class SimulationGenerationTests(unittest.TestCase):
    def test_config_and_output_contract(self):
        config = json.loads((SIM / "config.json").read_text())
        self.assertEqual(config["upstream"]["version"], "1.1.7")
        self.assertEqual(config["upstream"]["numpy_version"], "1.20.3")
        self.assertEqual(config["upstream"]["scipy_version"], "1.7.2")
        self.assertEqual(config["upstream"]["pandas_version"], "1.3.4")
        self.assertEqual(config["upstream"]["anndata_version"], "0.8.0")
        self.assertEqual(len(config["upstream"]["audited_sources"]), 3)
        self.assertTrue(all(len(source["sha256"]) == 64 for source in config["upstream"]["audited_sources"]))
        self.assertEqual(set(config["scenarios"]), {"murk", "branching", "hybrid"})
        self.assertEqual(config["scenarios"]["murk"]["output_name"], "sim_only_murk.h5ad")
        self.assertEqual(config["scenarios"]["branching"]["output_name"], "sim_only_branching.h5ad")
        self.assertEqual(config["scenarios"]["hybrid"]["output_name"], "sim_hybrid.h5ad")
        self.assertEqual(config["scenarios"]["branching"]["expected_shape"], [2000, 620])
        self.assertEqual(config["scenarios"]["hybrid"]["expected_shape"], [2000, 720])
        self.assertEqual(config["gene_parameters"]["boost_all"]["alpha2"], 4)
        self.assertEqual(config["gene_parameters"]["boost_all"]["beta2"], 2)
        self.assertEqual(config["gene_parameters"]["boost_all"]["gamma2"], 1)
        self.assertEqual(config["validation"]["seed"], 123)
        self.assertEqual(set(config["validation"]["output_sha256"]), {"murk", "branching", "hybrid"})
        self.assertTrue(all(len(value) == 64 for value in config["validation"]["output_sha256"].values()))

    def test_bins_are_ordered_and_match_branch_contract(self):
        labels, categories = _bin_labels(range(2000), [1000, 1000], "all")
        self.assertEqual(labels[:101], ["BranchA_0"] * 100 + ["BranchA_1"])
        self.assertEqual(labels[1000], "BranchB_0")
        self.assertTrue(categories == ["BranchA_%d" % i for i in range(10)] + ["BranchB_%d" % i for i in range(10)])
        labels, categories = _bin_labels(range(1000), [10, 1000], "branch_b")
        self.assertEqual(labels[-1], "BranchB_9")
        self.assertEqual(len(categories), 10)

    def test_truth_equations(self):
        velocity_s, velocity_u = _truth_layers(2.0, 3.0, 5.0, 1.0, .5)
        self.assertEqual(velocity_s, 2.0)
        self.assertEqual(velocity_u, 2.0)

    def test_random_alpha_is_drawn_before_expression_noise(self):
        import numpy as np
        import pandas as pd

        events = []

        class FakeSimulation:
            @staticmethod
            def simulate(kinetic_type, **kwargs):
                events.append((kinetic_type, kwargs))
                size = kwargs["path1_cell_number"] + kwargs["path2_cell_number"]
                return pd.DataFrame({"cellID": range(size)})

        def uniform(low, high):
            events.append(("uniform", (low, high)))
            return (low + high) / 2

        config = json.loads((SIM / "config.json").read_text())
        parameters = config["gene_parameters"]
        scenario = config["scenarios"]["hybrid"]
        with patch.object(np.random, "uniform", side_effect=uniform):
            frame = _one_sided(
                FakeSimulation,
                scenario,
                {
                    "quiet": parameters["one_sided_quiet"],
                    "boost": parameters["one_sided_boost"],
                },
                scenario["one_sided_alpha_ranges"],
            )

        self.assertEqual([event[0] for event in events], ["uniform", "uniform", "mono", "tran_boost"])
        self.assertEqual(events[-1][1]["alpha1"], 3.0)
        self.assertEqual(events[-1][1]["alpha2"], 10.0)
        self.assertEqual(frame["cellID"].tolist(), list(range(2000)))

    def test_generation_entrypoint_compiles_without_running_scientific_dependencies(self):
        self.assertTrue(__import__("py_compile").compile(str(SIM / "generate.py"), doraise=True))

    def test_generation_entrypoint_has_no_server_paths(self):
        text = (SIM / "generate.py").read_text()
        self.assertNotIn("/mnt/", text)
        self.assertNotIn("/Users/", text)
