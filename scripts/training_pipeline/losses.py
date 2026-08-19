"""Reusable contact-map losses shared by the training pipeline and tests."""

from __future__ import annotations

import torch
import torch.nn.functional as F


def balanced_binary_cross_entropy_per_map(
    logits: torch.Tensor,
    targets: torch.Tensor,
    known_mask: torch.Tensor,
    omega: float,
) -> torch.Tensor:
    """Return class-balanced BCE averaged equally across supervised maps.

    Tensors share shape ``[B, ...]``: either contact maps ``[B,L1,L2]`` or
    chain-wise interface labels ``[B,L]``. ``known_mask`` selects supervised
    binary targets. Within each example, positive and negative mean BCE receive
    weights ``omega`` and ``1 - omega``. If one class is absent, the active
    class is renormalized rather than shrinking the example's contribution.
    """
    if logits.shape != targets.shape or logits.shape != known_mask.shape:
        raise ValueError("logits, targets, and known_mask must have matching shapes")

    per_map_losses = []
    for logits_i, targets_i, known_i in zip(logits, targets, known_mask):
        known_logits = logits_i[known_i]
        known_targets = targets_i[known_i]
        if known_targets.numel() == 0:
            continue

        positive = known_targets == 1
        negative = known_targets == 0
        weighted_loss = logits_i.new_zeros(())
        active_weight = 0.0

        if positive.any():
            weighted_loss = weighted_loss + omega * F.binary_cross_entropy_with_logits(
                known_logits[positive],
                torch.ones_like(known_logits[positive]),
            )
            active_weight += omega

        if negative.any():
            weighted_loss = weighted_loss + (1.0 - omega) * F.binary_cross_entropy_with_logits(
                known_logits[negative],
                torch.zeros_like(known_logits[negative]),
            )
            active_weight += 1.0 - omega

        if active_weight > 0:
            per_map_losses.append(weighted_loss / active_weight)

    if per_map_losses:
        return torch.stack(per_map_losses).mean()
    return logits.sum() * 0.0
