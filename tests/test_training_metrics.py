"""Tests for streaming interaction/contact average-precision metrics."""

from __future__ import annotations

import unittest

import torch

from scripts.training_pipeline.metrics import (
    BinnedBinaryAveragePrecision,
    MacroContactAveragePrecision,
)


class BinnedBinaryAveragePrecisionTest(unittest.TestCase):
    """Verify histogram AP values, streaming behavior, and validation."""

    def test_streaming_updates_match_expected_average_precision(self) -> None:
        # Every score occupies a distinct bin, so binned AP is exact here. In
        # descending order the positive labels occur at ranks one and three:
        # AP = (1/1 + 2/3) / 2 = 5/6.
        predictions = torch.tensor([0.9, 0.8, 0.7, 0.1])
        targets = torch.tensor([1, 0, 1, 0])

        one_update = BinnedBinaryAveragePrecision(num_bins=16)
        one_update.update(predictions, targets)

        streamed = BinnedBinaryAveragePrecision(num_bins=16)
        streamed.update(predictions[:2], targets[:2])
        streamed.update(predictions[2:], targets[2:])

        expected = torch.tensor(5.0 / 6.0, dtype=torch.float64)
        torch.testing.assert_close(one_update.compute(), expected)
        torch.testing.assert_close(streamed.compute(), expected)
        torch.testing.assert_close(streamed.compute(), one_update.compute())

    def test_same_bin_predictions_are_treated_as_tied(self) -> None:
        # With two bins, both high probabilities share one threshold. The first
        # threshold therefore contains one positive and one negative, producing
        # precision 1/2 at full recall and AP 1/2.
        metric = BinnedBinaryAveragePrecision(num_bins=2)
        metric.update(
            torch.tensor([0.90, 0.75, 0.10]),
            torch.tensor([1, 0, 0]),
        )

        torch.testing.assert_close(metric.compute(), torch.tensor(0.5, dtype=torch.float64))

    def test_all_negative_and_empty_inputs_return_zero(self) -> None:
        metric = BinnedBinaryAveragePrecision(num_bins=8)
        metric.update(torch.tensor([], dtype=torch.float32), torch.tensor([], dtype=torch.long))
        metric.update(torch.tensor([0.8, 0.2]), torch.tensor([0, 0]))

        torch.testing.assert_close(metric.compute(), torch.tensor(0.0, dtype=torch.float64))

    def test_reset_clears_both_histograms(self) -> None:
        metric = BinnedBinaryAveragePrecision(num_bins=8)
        metric.update(torch.tensor([0.9, 0.1]), torch.tensor([1, 0]))
        self.assertGreater(metric.compute().item(), 0.0)

        metric.reset()

        self.assertEqual(metric.positive_histogram.sum().item(), 0)
        self.assertEqual(metric.negative_histogram.sum().item(), 0)

    def test_rejects_logits_and_unknown_targets(self) -> None:
        metric = BinnedBinaryAveragePrecision(num_bins=8)
        with self.assertRaisesRegex(ValueError, "probabilities"):
            metric.update(torch.tensor([1.1]), torch.tensor([1]))
        with self.assertRaisesRegex(ValueError, "binary"):
            metric.update(torch.tensor([0.5]), torch.tensor([-1]))


class MacroContactAveragePrecisionTest(unittest.TestCase):
    """Verify per-map exact AP, masks, inclusion, and macro weighting."""

    def test_macro_ap_masks_unknowns_excludes_maps_and_counts_all_negative(self) -> None:
        # Map 0 has positive labels at ranks one and four after unknown cells are
        # removed: AP = (1 + 2/4) / 2 = 0.75. Map 1 is included but all-negative,
        # so it contributes zero. Map 2 is excluded despite its positive label.
        predictions = torch.tensor(
            [
                [[0.90, 0.80, 0.10], [0.70, 0.20, 0.99]],
                [[0.90, 0.20, 0.10], [0.80, 0.30, 0.05]],
                [[0.99, 0.01, 0.01], [0.01, 0.01, 0.01]],
            ]
        )
        targets = torch.tensor(
            [
                [[1, 0, 1], [0, -1, -1]],
                [[0, 0, 0], [0, 0, 0]],
                [[1, 0, 0], [0, 0, 0]],
            ]
        )
        known_mask = targets >= 0
        include_map = torch.tensor([True, True, False])

        metric = MacroContactAveragePrecision()
        metric.update(predictions, targets, known_mask, include_map)

        self.assertEqual(metric.evaluated_map_count.item(), 2)
        torch.testing.assert_close(
            metric.average_precision_sum,
            torch.tensor(0.75, dtype=torch.float64),
        )
        torch.testing.assert_close(metric.compute(), torch.tensor(0.375, dtype=torch.float64))

    def test_streaming_updates_weight_maps_not_batches(self) -> None:
        # The first update contains one perfect map (AP=1). The second contains
        # two all-negative maps (AP=0 each). A mean of batch means would be 0.5,
        # whereas the correct map-macro result is 1/3.
        metric = MacroContactAveragePrecision()
        metric.update(
            torch.tensor([[[0.9, 0.1]]]),
            torch.tensor([[[1, 0]]]),
            torch.tensor([[[True, True]]]),
            torch.tensor([True]),
        )
        metric.update(
            torch.tensor([[[0.8, 0.2]], [[0.7, 0.1]]]),
            torch.tensor([[[0, 0]], [[0, 0]]]),
            torch.tensor([[[True, True]], [[True, True]]]),
            torch.tensor([True, True]),
        )

        torch.testing.assert_close(
            metric.compute(),
            torch.tensor(1.0 / 3.0, dtype=torch.float64),
        )

    def test_included_map_without_known_cells_is_skipped(self) -> None:
        metric = MacroContactAveragePrecision()
        metric.update(
            torch.tensor([[[0.9, 0.1]]]),
            torch.tensor([[[-1, -1]]]),
            torch.tensor([[[False, False]]]),
            torch.tensor([True]),
        )

        self.assertEqual(metric.evaluated_map_count.item(), 0)
        self.assertEqual(metric.compute().item(), 0.0)

    def test_reset_clears_sum_and_count(self) -> None:
        metric = MacroContactAveragePrecision()
        metric.update(
            torch.tensor([[[0.9, 0.1]]]),
            torch.tensor([[[1, 0]]]),
            torch.tensor([[[True, True]]]),
            torch.tensor([True]),
        )

        metric.reset()

        self.assertEqual(metric.average_precision_sum.item(), 0.0)
        self.assertEqual(metric.evaluated_map_count.item(), 0)

    def test_validates_contact_shapes_and_masks(self) -> None:
        metric = MacroContactAveragePrecision()
        predictions = torch.tensor([[[0.5]]])
        targets = torch.tensor([[[1]]])

        with self.assertRaisesRegex(ValueError, "known_mask must be a boolean"):
            metric.update(
                predictions,
                targets,
                torch.ones_like(targets),
                torch.tensor([True]),
            )
        with self.assertRaisesRegex(ValueError, "include_map must be shaped"):
            metric.update(
                predictions,
                targets,
                torch.ones_like(targets, dtype=torch.bool),
                torch.tensor([[True]]),
            )


if __name__ == "__main__":
    unittest.main()
