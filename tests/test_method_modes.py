"""Protect role-dependent training, velocity scaling, and output isolation."""
import argparse
import ast
import contextlib
import io
import importlib.machinery
import sys
import json
import tempfile
import unittest
from datetime import datetime, timezone
from inspect import signature
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import anndata as ad
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def load_workflow(name, function_names, **namespace):
    # Method environments are separate. Exercise the direct script with a fixed
    # model boundary rather than importing or training upstream GPU packages.
    path = ROOT / "methods" / f"{name}.py"
    tree = ast.parse(path.read_text())
    tree.body = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in function_names
    ]
    exec(compile(tree, str(path), "exec"), namespace)
    return SimpleNamespace(**namespace)


class ArrayTensor:
    def __init__(self, values):
        self.values = values

    def detach(self):
        return self

    def cpu(self):
        return self

    def numpy(self):
        return self.values


class PackageImportTests(unittest.TestCase):
    def test_direct_launchers_resolve_installed_packages_instead_of_sibling_scripts(self):
        stdlib = {"argparse", "datetime", "importlib", "inspect", "json", "os",
                  "pathlib", "random", "sys", "time"}
        methods = ROOT / "methods"
        checked = 0
        with tempfile.TemporaryDirectory() as temp:
            package_root = Path(temp)
            for path in methods.glob("*.py"):
                tree = ast.parse(path.read_text())
                imports = [
                    name.split(".")[0]
                    for node in tree.body
                    for name in (
                        [alias.name for alias in node.names] if isinstance(node, ast.Import)
                        else [node.module] if isinstance(node, ast.ImportFrom) else []
                    )
                ]
                shadowed = {name for name in imports if (methods / f"{name}.py").exists()}
                if not shadowed:
                    continue
                for name in shadowed:
                    (package_root / f"{name}.py").touch()
                first_external = next(
                    i for i, node in enumerate(tree.body)
                    if isinstance(node, ast.Import) and any(
                        alias.name.split(".")[0] not in stdlib for alias in node.names
                    ) or isinstance(node, ast.ImportFrom) and node.module.split(".")[0] not in stdlib
                )
                prefix = ast.Module(body=tree.body[:first_external], type_ignores=[])
                with self.subTest(launcher=path.name), patch.object(
                    sys, "path", [str(methods), str(package_root), *sys.path]
                ):
                    exec(compile(prefix, str(path), "exec"), {"__file__": str(path)})
                    for name in shadowed:
                        spec = importlib.machinery.PathFinder.find_spec(name, sys.path)
                        self.assertEqual(spec.origin, str(package_root / f"{name}.py"))
                checked += 1
        self.assertEqual(checked, 15)


class VeloVIRoleTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.adata = ad.AnnData(np.ones((3, 2)))
        for key in ("spliced", "unspliced", "Ms", "Mu"):
            self.adata.layers[key] = np.ones((3, 2))
        self.adata.obs["clusters"] = ["a", "a", "b"]
        self.adata.obsm["X_umap"] = np.zeros((3, 2))
        self.latent = np.array([[1., 1.], [2., 3.], [1., 4.]])
        self.spliced = np.array([[10., 20.], [30., 40.], [50., 60.]])
        owner = self

        class Model:
            setup_anndata = Mock()

            def __init__(self, adata):
                self.train = Mock(spec=self.train)
                self.get_latent_time = Mock(return_value=owner.latent)
                self.get_velocity = Mock(side_effect=lambda **kw: (
                    owner.spliced + 100 if kw.get("velo_mode") == "unspliced" else owner.spliced
                ))
                self.get_rates = Mock(return_value={key: np.array([20., 30.]) for key in ("alpha", "beta", "gamma")})
                self.module = SimpleNamespace(switch_time_unconstr=np.zeros(2))
                owner.model = self

            def train(self, max_epochs=500, batch_size=256):
                pass

        self.cuda = Mock(return_value=True)

        def graph(adata, **kwargs):
            adata.uns["velocity_graph"] = np.eye(3)

        def stream(adata, **kwargs):
            adata.obsm["velocity_umap"] = np.zeros((3, 2))

        self.workflow = load_workflow(
            "velovi", {"run_dataset", "main"},
            argparse=argparse, signature=signature, Path=Path, np=np, VELOVI=Model,
            velovi=SimpleNamespace(__version__="0.3.1"),
            torch=SimpleNamespace(
                cuda=SimpleNamespace(is_available=self.cuda),
                nn=SimpleNamespace(functional=SimpleNamespace(
                    softplus=lambda values: ArrayTensor(np.logaddexp(0, values)),
                )),
            ),
            sc=SimpleNamespace(read=Mock(return_value=self.adata)),
            scv=SimpleNamespace(tl=SimpleNamespace(velocity_graph=graph),
                                pl=SimpleNamespace(velocity_embedding_stream=stream),
                                settings=SimpleNamespace()),
            plt=SimpleNamespace(close=Mock()), preprocess_data=lambda adata: adata,
            DATASETS={"Pancreas": {"input": "Pancreas/Pancreas.h5ad", "label_key": "clusters",
                                  "embedding_key": None, "preprocessed": True, "expected_shape": (3, 2)}},
            DATASETS_TO_RUN=["Pancreas"], DATA_ROOT=self.root,
            GENE_MAX_EPOCHS=1000, GENE_GPU_BATCH_SIZE=64, GENE_CPU_BATCH_SIZE=32,
            N_POSTERIOR_SAMPLES=25, N_JOBS=8,
        )

    def run_role(self, role):
        with contextlib.redirect_stdout(io.StringIO()):
            self.workflow.run_dataset("Pancreas", role, self.root, self.root / role)
        return ad.read_h5ad(self.root / role / "Pancreas/veloVI_Pancreas_all.h5ad")

    def test_cell_preserves_default_training_and_per_gene_scaling(self):
        result = self.run_role("cell")
        self.model.train.assert_called_once_with()
        self.model.get_latent_time.assert_called_once_with(n_samples=25)
        self.model.get_velocity.assert_called_once_with(n_samples=25, velo_statistic="mean")
        scaling = np.array([10., 5.])
        np.testing.assert_array_equal(result.layers["velocity"], self.spliced / scaling)
        np.testing.assert_array_equal(result.layers["fit_t"], self.latent * scaling)
        np.testing.assert_array_equal(result.var["fit_alpha"], np.array([20., 30.]) / scaling)
        self.assertNotIn("velocity_u", result.layers)
        self.assertEqual(result.uns["benchmark_run"]["max_epochs"], 500)
        self.assertEqual(result.uns["benchmark_run"]["batch_size"], 256)

    def test_gene_gpu_preserves_unscaled_velocities_and_training(self):
        result = self.run_role("gene")
        self.model.train.assert_called_once_with(max_epochs=1000, batch_size=64)
        self.model.get_latent_time.assert_not_called()
        np.testing.assert_array_equal(result.layers["velocity"], self.spliced)
        np.testing.assert_array_equal(result.layers["velocity_u"], self.spliced + 100)
        self.assertEqual(result.uns["benchmark_run"]["evaluation_role"], "gene")
        self.assertEqual(result.uns["benchmark_run"]["batch_size"], 64)

    def test_gene_cpu_uses_historical_batch_size(self):
        self.cuda.return_value = False
        self.run_role("gene")
        self.model.train.assert_called_once_with(max_epochs=1000, batch_size=32)

    def test_nonpositive_cell_time_is_rejected(self):
        self.latent[:, 0] = 0
        with self.assertRaisesRegex(ValueError, "finite and positive"):
            self.run_role("cell")

    def test_unknown_role_is_rejected_before_input_read(self):
        with self.assertRaisesRegex(ValueError, "unknown veloVI evaluation role"):
            self.run_role("kinetics")
        self.workflow.sc.read.assert_not_called()

    def test_existing_output_is_rejected_before_training(self):
        output = self.root / "existing"
        output.mkdir()
        sentinel = output / "keep"
        sentinel.write_text("existing result")
        with patch("sys.argv", ["velovi.py", "--evaluation-role", "gene", "--output-root", str(output)]):
            with self.assertRaisesRegex(RuntimeError, "output root is not empty"):
                self.workflow.main()
        self.assertEqual(sentinel.read_text(), "existing result")
        self.workflow.sc.read.assert_not_called()


class DeepVeloRoleTests(unittest.TestCase):
    def test_role_controls_training_head_and_recorded_configuration(self):
        for role, predict_u in (("cell", False), ("gene", True), ("kinetics", True)):
            with self.subTest(role=role), tempfile.TemporaryDirectory() as temp:
                output = Path(temp) / role
                run = Mock()
                workflow = load_workflow(
                    "deepvelo", {"main"}, Path=Path, json=json, datetime=datetime, timezone=timezone,
                    parse_args=lambda: SimpleNamespace(
                        datasets=["Pancreas"], evaluation_role=role, output_root=output,
                        reset_seed_per_dataset=False,
                    ),
                    DATASETS_TO_RUN=["Pancreas"], DATASETS={"Pancreas": {}},
                    scv=SimpleNamespace(__version__="0.2.5", settings=SimpleNamespace(), set_figure_params=Mock()),
                    seed_everything=Mock(), train=Mock(), run_dataset=run, SEED=123, N_JOBS=8,
                )
                with contextlib.redirect_stdout(io.StringIO()):
                    workflow.main()
                self.assertEqual(run.call_args.args[2:4], (predict_u, role))
                metadata = json.loads((output / "run_config.json").read_text())
                self.assertEqual(metadata["pred_unspliced"], predict_u)
                self.assertEqual(metadata["evaluation_role"], role)


if __name__ == "__main__":
    unittest.main()
