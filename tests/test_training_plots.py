"""Tests for AUPRC history plots generated from Lightning CSV logs."""

from __future__ import annotations

import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

import pandas as pd

from scripts.training_pipeline.execute_pipeline import save_auprc_plots


class SaveAuprcPlotsTest(unittest.TestCase):
    """Verify sparse logger rows produce the appropriate run-level PNGs."""

    def test_sparse_interaction_and_contact_metrics_write_both_plots(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary_path = Path(temporary_directory)
            metrics_path = temporary_path / "metrics.csv"
            output_dir = temporary_path / "run"

            # CSVLogger stores different epoch-level metrics on sparse rows. The
            # duplicate epoch 1 also verifies that plotting tolerates repeated
            # rows from logging/resume behavior and keeps the latest value.
            pd.DataFrame(
                {
                    "epoch": [1, 0, 1, 2],
                    "val_interaction_auprc": [0.70, 0.60, 0.75, None],
                    "test_interaction_auprc": [None, None, None, 0.80],
                    "val_contact_auprc_macro": [0.30, 0.20, 0.35, None],
                    "test_contact_auprc_macro": [None, None, None, 0.40],
                    "val_contact_auprc_micro": [0.10, 0.05, 0.12, None],
                    "test_contact_auprc_micro": [None, None, None, 0.15],
                }
            ).to_csv(metrics_path, index=False)

            save_auprc_plots(metrics_path, output_dir)

            interaction_plot = output_dir / "interaction_auprc_history.png"
            contact_plot = output_dir / "contact_auprc_history.png"
            self.assertTrue(interaction_plot.is_file())
            self.assertTrue(contact_plot.is_file())
            self.assertGreater(interaction_plot.stat().st_size, 0)
            self.assertGreater(contact_plot.stat().st_size, 0)

    def test_interaction_only_run_skips_contact_plot(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary_path = Path(temporary_directory)
            metrics_path = temporary_path / "metrics.csv"
            output_dir = temporary_path / "run"
            pd.DataFrame(
                {
                    "epoch": [0, 1, 2],
                    "val_interaction_auprc": [0.55, 0.65, None],
                    "test_interaction_auprc": [None, None, 0.70],
                }
            ).to_csv(metrics_path, index=False)

            output = StringIO()
            with redirect_stdout(output):
                save_auprc_plots(metrics_path, output_dir)

            self.assertTrue((output_dir / "interaction_auprc_history.png").is_file())
            self.assertFalse((output_dir / "contact_auprc_history.png").exists())
            self.assertIn("Skipping contact-map auprc by epoch plot", output.getvalue())

    def test_test_only_metric_ignores_malformed_validation_rows(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary_path = Path(temporary_directory)
            metrics_path = temporary_path / "metrics.csv"
            output_dir = temporary_path / "run"

            # Invalid epochs/values must not abort post-training artifact
            # generation. With no valid validation history, the final test
            # result is rendered as a categorical point rather than on an epoch
            # axis or as a fabricated validation observation.
            pd.DataFrame(
                {
                    "epoch": ["malformed", 2],
                    "step": ["malformed", 100],
                    "val_interaction_auprc": ["not-a-number", 1.2],
                    "test_interaction_auprc": [None, 0.7],
                }
            ).to_csv(metrics_path, index=False)

            save_auprc_plots(metrics_path, output_dir)

            interaction_plot = output_dir / "interaction_auprc_history.png"
            self.assertTrue(interaction_plot.is_file())
            self.assertGreater(interaction_plot.stat().st_size, 0)
            self.assertFalse((output_dir / "contact_auprc_history.png").exists())


if __name__ == "__main__":
    unittest.main()
