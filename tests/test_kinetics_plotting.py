import json
import unittest
from pathlib import Path


PLOTTING_ROOT = (
    Path(__file__).resolve().parents[1]
    / "evaluation"
    / "kinetics_recovery"
    / "plotting"
)


class KineticsPlottingContractTests(unittest.TestCase):
    def test_stream_config_has_explicit_method_keys(self):
        config = json.loads((PLOTTING_ROOT / "stream_config.json").read_text())
        specs = config["organoid"]
        self.assertEqual([spec["method"] for spec in specs], ["cellDancer", "DeepKINET", "DeepVelo"])
        self.assertEqual(
            [spec["output_stem"] for spec in specs],
            ["celldancer", "DeepKINET", "deepvelo"],
        )
        self.assertEqual(
            [spec["velocity_layer"] for spec in specs],
            ["velocity", "DeepKINET_velocity", "velocity"],
        )
        manual_rpe1_specs = config["rpe1_manual"]
        self.assertEqual(
            [spec["method"] for spec in manual_rpe1_specs],
            ["Dynamo", "cellDancer", "DeepVelo", "TSvelo"],
        )
        rpe1_specs = config["rpe1_deepkinet"]
        self.assertEqual([spec["method"] for spec in rpe1_specs], ["DeepKINET"])
        self.assertEqual([spec["output_stem"] for spec in rpe1_specs], ["deepkinet"])
        self.assertEqual(
            [spec["velocity_layer"] for spec in rpe1_specs],
            ["DeepKINET_velocity"],
        )

    def test_stream_paths_are_relative(self):
        config = json.loads((PLOTTING_ROOT / "stream_config.json").read_text())
        for panel in config.values():
            for spec in panel:
                path = Path(spec.get("result") or spec.get("asset"))
                self.assertFalse(path.is_absolute())
                self.assertNotIn("..", path.parts)

    def test_kinetic_class_panels_use_mapped_trajectories(self):
        config = json.loads(
            (PLOTTING_ROOT / "figure5_config.json").read_text()
        )
        venn = config["inputs"]["rpe1"]["cellDancer_venn"]
        self.assertEqual(
            set(venn),
            {"reference_alpha", "reference_gamma", "prediction_alpha", "prediction_gamma"},
        )
        self.assertTrue(
            all("static" not in Path(path).parts for path in venn.values())
        )


if __name__ == "__main__":
    unittest.main()
