"""Memory-efficient average-precision metrics for interaction training.

Two aggregation levels are needed by the joint interaction/contact pipeline:

* :class:`BinnedBinaryAveragePrecision` pools a potentially very large set of
  binary predictions into a fixed-size score histogram. It is intended for
  cell-micro contact-map evaluation, where retaining every residue-pair score
  until the end of an epoch would consume hundreds of megabytes.
* :class:`MacroContactAveragePrecision` computes exact average precision for
  each selected structural contact map and then gives every evaluated map equal
  weight. Unknown or padded residue pairs are omitted through an explicit mask.

Both metrics store only additive states. TorchMetrics can therefore sum their
states across distributed workers before ``compute`` without gathering raw
predictions or making results depend on batch boundaries.
"""

from __future__ import annotations

import torch
from torchmetrics import Metric
from torchmetrics.functional.classification import binary_average_precision


class BinnedBinaryAveragePrecision(Metric):
    """Approximate binary average precision with a streaming score histogram.

    ``predictions`` must be probabilities in ``[0, 1]`` and ``targets`` must be
    binary labels with the same shape. Inputs may have any number of dimensions;
    all values are pooled into one micro-aggregated precision-recall curve.

    The probability interval is divided into ``num_bins`` equal-width bins.
    ``update`` assigns every observation to one class-specific bin using one
    ``torch.bincount`` call, giving O(N) update time and O(num_bins) persistent
    memory. ``compute`` traverses bins from high to low probability and returns
    the step-wise precision-recall integral, commonly called average precision
    (AP) or AUPRC. Predictions within one bin are treated as tied, so the result
    approaches exact AP as the number of bins increases.

    Unknown labels must be removed by the caller before ``update``. This keeps
    the metric generic and prevents contact-map padding conventions such as
    ``-1`` from being confused with binary targets.

    Args:
        num_bins: Number of equally spaced probability bins. The default 4096
            keeps the metric state small while resolving probabilities to about
            2.4e-4. At least two bins are required.
        **kwargs: Additional keyword arguments accepted by
            :class:`torchmetrics.Metric`.
    """

    is_differentiable = False
    higher_is_better = True
    full_state_update = False

    def __init__(self, num_bins: int = 4096, **kwargs: object) -> None:
        super().__init__(**kwargs)
        if isinstance(num_bins, bool) or not isinstance(num_bins, int) or num_bins < 2:
            raise ValueError("num_bins must be an integer greater than or equal to 2")

        self.num_bins = num_bins

        # Positive and negative counts are separate because precision requires
        # cumulative true- and false-positive counts at every score threshold.
        # Integer states preserve exact observation counts even for very large
        # contact-map validation sets. DDP workers sum both histograms.
        self.add_state(
            "positive_histogram",
            default=torch.zeros(num_bins, dtype=torch.long),
            dist_reduce_fx="sum",
        )
        self.add_state(
            "negative_histogram",
            default=torch.zeros(num_bins, dtype=torch.long),
            dist_reduce_fx="sum",
        )

    def update(self, predictions: torch.Tensor, targets: torch.Tensor) -> None:
        """Add arbitrary-shaped probability predictions and binary targets."""
        if predictions.shape != targets.shape:
            raise ValueError(
                "predictions and targets must have identical shapes, got "
                f"{tuple(predictions.shape)} and {tuple(targets.shape)}"
            )

        flat_predictions = predictions.detach().reshape(-1)
        flat_targets = targets.detach().reshape(-1)
        if flat_predictions.numel() == 0:
            return

        # The bin mapping assumes calibrated probability bounds. Rejecting
        # logits and non-finite values is safer than silently clamping them into
        # endpoint bins, which would produce a plausible but incorrect PR curve.
        if not flat_predictions.is_floating_point():
            raise ValueError("predictions must be floating-point probabilities in [0, 1]")
        if bool((~torch.isfinite(flat_predictions)).any()):
            raise ValueError("predictions must contain only finite values")
        if bool(((flat_predictions < 0) | (flat_predictions > 1)).any()):
            raise ValueError("predictions must be probabilities in [0, 1]")

        binary_targets = (flat_targets == 0) | (flat_targets == 1)
        if not bool(binary_targets.all()):
            raise ValueError("targets must contain only binary values 0 and 1")
        flat_targets = flat_targets.to(torch.long)

        # Multiplication maps p in [i/B, (i+1)/B) to bin i. Probability 1.0
        # initially maps to B and is explicitly placed in the final bin.
        bin_indices = (flat_predictions * self.num_bins).to(torch.long)
        bin_indices.clamp_(max=self.num_bins - 1)

        # Encode the target class as an offset so a single bincount pass builds
        # both histograms: [negative bins | positive bins]. This avoids an
        # O(num_bins * N) comparison against every possible threshold.
        class_bin_indices = bin_indices + flat_targets * self.num_bins
        class_bin_counts = torch.bincount(
            class_bin_indices,
            minlength=2 * self.num_bins,
        )
        self.negative_histogram += class_bin_counts[: self.num_bins]
        self.positive_histogram += class_bin_counts[self.num_bins :]

    def compute(self) -> torch.Tensor:
        """Return histogram-approximated micro average precision."""
        # A descending sweep adds every observation in a score bin at once.
        # Consequently, observations quantized to the same bin are handled as
        # tied predictions at one precision-recall operating point.
        positive_descending = self.positive_histogram.flip(0).to(torch.float64)
        negative_descending = self.negative_histogram.flip(0).to(torch.float64)
        total_positive = positive_descending.sum()

        cumulative_true_positive = positive_descending.cumsum(dim=0)
        cumulative_false_positive = negative_descending.cumsum(dim=0)
        predicted_positive = cumulative_true_positive + cumulative_false_positive
        precision = torch.where(
            predicted_positive > 0,
            cumulative_true_positive / predicted_positive.clamp_min(1.0),
            torch.zeros_like(predicted_positive),
        )

        # AP is sum(delta_recall * precision). A bin adds
        # positive_count / total_positive recall. clamp_min makes the no-positive
        # case numerically safe; torch.where defines its result explicitly as 0.
        recall_increment = positive_descending / total_positive.clamp_min(1.0)
        average_precision = (precision * recall_increment).sum()
        return torch.where(
            total_positive > 0,
            average_precision,
            torch.zeros_like(average_precision),
        )


class MacroContactAveragePrecision(Metric):
    """Compute exact AP per selected contact map and average maps equally.

    Args passed to :meth:`update` have the following conventions:

    * ``predictions``: contact probabilities shaped ``[B, L1, L2]`` in
      ``[0, 1]``;
    * ``targets``: contact labels with the same shape, normally ``1`` for a
      contact, ``0`` for a known non-contact, and ``-1`` for unknown/padding;
    * ``known_mask``: Boolean tensor of the same shape. Only ``True`` cells enter
      a map's precision-recall curve, so target values at other cells are ignored;
    * ``include_map``: Boolean tensor shaped ``[B]`` selecting real structural
      maps. It lets callers exclude optional synthetic all-zero maps used only
      for negative-example training.

    AP is evaluated exactly with TorchMetrics' functional
    ``binary_average_precision`` independently for every selected map containing
    at least one known cell. An included map with known cells but no positive
    contacts contributes AP=0 and remains in the macro denominator. An included
    map with no known cells is skipped because it has no evaluable labels.

    Only an AP sum and evaluated-map count persist across steps. Both use sum
    reduction under DDP, so ``compute`` returns a true map-level macro average
    without retaining any residue-pair predictions after ``update`` returns.
    """

    is_differentiable = False
    higher_is_better = True
    full_state_update = False

    def __init__(self, **kwargs: object) -> None:
        super().__init__(**kwargs)
        self.add_state(
            "average_precision_sum",
            default=torch.tensor(0.0, dtype=torch.float64),
            dist_reduce_fx="sum",
        )
        self.add_state(
            "evaluated_map_count",
            default=torch.tensor(0, dtype=torch.long),
            dist_reduce_fx="sum",
        )

    def update(
        self,
        predictions: torch.Tensor,
        targets: torch.Tensor,
        known_mask: torch.Tensor,
        include_map: torch.Tensor,
    ) -> None:
        """Add one batch of padded contact maps to the macro AP state."""
        if predictions.ndim != 3:
            raise ValueError(
                "predictions, targets, and known_mask must be shaped [B, L1, L2], "
                f"got predictions with shape {tuple(predictions.shape)}"
            )
        if predictions.shape != targets.shape or predictions.shape != known_mask.shape:
            raise ValueError(
                "predictions, targets, and known_mask must have identical shapes, got "
                f"{tuple(predictions.shape)}, {tuple(targets.shape)}, and {tuple(known_mask.shape)}"
            )
        if include_map.ndim != 1 or include_map.shape[0] != predictions.shape[0]:
            raise ValueError(
                "include_map must be shaped [B], got "
                f"{tuple(include_map.shape)} for batch size {predictions.shape[0]}"
            )
        if known_mask.dtype != torch.bool:
            raise ValueError("known_mask must be a boolean tensor")
        if include_map.dtype != torch.bool:
            raise ValueError("include_map must be a boolean tensor")
        if not predictions.is_floating_point():
            raise ValueError("predictions must be floating-point probabilities in [0, 1]")

        predictions = predictions.detach()
        targets = targets.detach()
        known_mask = known_mask.detach()
        include_map = include_map.detach()

        # Validate only cells that can enter the metric. Padding and excluded
        # maps may legitimately contain sentinel targets or placeholder scores.
        evaluation_mask = known_mask & include_map[:, None, None]
        evaluated_predictions = predictions[evaluation_mask]
        evaluated_targets = targets[evaluation_mask]
        if evaluated_predictions.numel() > 0:
            if bool((~torch.isfinite(evaluated_predictions)).any()):
                raise ValueError("known predictions in included maps must be finite")
            if bool(((evaluated_predictions < 0) | (evaluated_predictions > 1)).any()):
                raise ValueError("known predictions in included maps must be probabilities in [0, 1]")
            binary_targets = (evaluated_targets == 0) | (evaluated_targets == 1)
            if not bool(binary_targets.all()):
                raise ValueError("known targets in included maps must contain only 0 and 1")

        batch_ap_sum = self.average_precision_sum.new_zeros(())
        batch_map_count = self.evaluated_map_count.new_zeros(())

        for map_index in include_map.nonzero(as_tuple=False).flatten():
            map_known_mask = known_mask[map_index]
            if not bool(map_known_mask.any()):
                # No threshold can be evaluated without at least one known
                # target, so this map must not enter the macro denominator.
                continue

            map_predictions = predictions[map_index][map_known_mask]
            map_targets = targets[map_index][map_known_mask].to(torch.long)
            batch_map_count += 1

            if bool((map_targets == 1).any()):
                map_average_precision = binary_average_precision(
                    map_predictions,
                    map_targets,
                )
                batch_ap_sum += map_average_precision.to(batch_ap_sum.dtype)
            # A known, included all-negative map is a legitimate difficult case.
            # Its AP is explicitly zero, while its count added above ensures it
            # lowers the macro average instead of disappearing from evaluation.

        self.average_precision_sum += batch_ap_sum
        self.evaluated_map_count += batch_map_count

    def compute(self) -> torch.Tensor:
        """Return exact map-macro average precision, or zero with no maps."""
        map_count = self.evaluated_map_count.to(self.average_precision_sum.dtype)
        macro_average_precision = self.average_precision_sum / map_count.clamp_min(1.0)
        return torch.where(
            self.evaluated_map_count > 0,
            macro_average_precision,
            torch.zeros_like(macro_average_precision),
        )
