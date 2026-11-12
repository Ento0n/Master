"""Exact contact-map average precision with equal weight for each map.

Each selected map contributes one metric value after unknown and padded cells
are removed. Only metric sums and map counts persist across updates;
distributed workers sum these states without retaining residue-pair predictions.
"""

from __future__ import annotations

import torch
from torchmetrics import Metric
from torchmetrics.functional.classification import binary_average_precision


class ContactMetricMean(Metric):
    """Accumulate batch means using the number of evaluated contact-maps.

    ``batch_mean`` is a scalar mean of per-map scores and ``map_count`` is the
    number of maps included in that mean. Batches/ranks without known cells use
    zero for both inputs, which registers the same epoch reductions everywhere
    without adding weight. An epoch with no evaluable maps returns zero.
    """

    is_differentiable = False
    higher_is_better = True
    full_state_update = False

    def __init__(self, **kwargs: object) -> None:
        super().__init__(**kwargs)
        self.add_state("score_sum", default=torch.tensor(0.0, dtype=torch.float64), dist_reduce_fx="sum")
        self.add_state("map_count", default=torch.tensor(0, dtype=torch.long), dist_reduce_fx="sum")

    def update(self, batch_mean: torch.Tensor, map_count: int) -> None:
        """Add one batch without retaining its computation graph or predictions."""
        self.score_sum += batch_mean.detach().to(self.score_sum.dtype) * map_count
        self.map_count += map_count

    def compute(self) -> torch.Tensor:
        """Return the mean across all evaluated maps in all batches/workers."""
        return self.score_sum / self.map_count.clamp_min(1)


class ContactAveragePrecision(Metric):
    """Compute exact AP per selected contact-map and average maps equally.

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
    contacts contributes AP=0 and remains in the map count. An included
    map with no known cells is skipped because it has no evaluable labels.

    Only an AP sum and evaluated-map count persist across steps. Both use sum
    reduction under DDP, so ``compute`` returns a mean across all evaluated maps
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
        """Add one batch of padded contact-maps to the AP state."""
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
                # target, so this map must not enter the map count.
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
            # lowers the mean instead of disappearing from evaluation.

        self.average_precision_sum += batch_ap_sum
        self.evaluated_map_count += batch_map_count

    def compute(self) -> torch.Tensor:
        """Return exact average precision across maps, or zero with no maps."""
        map_count = self.evaluated_map_count.to(self.average_precision_sum.dtype)
        average_precision = self.average_precision_sum / map_count.clamp_min(1.0)
        return torch.where(
            self.evaluated_map_count > 0,
            average_precision,
            torch.zeros_like(average_precision),
        )
