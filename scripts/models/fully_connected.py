"""Small fully connected model for one-hot encoded protein interaction pairs."""

from __future__ import annotations

import torch
from torch import nn


class FullyConnectedInteractionModel(nn.Module):
    """Three Linear layers that classify whether two proteins interact.

    The input is expected to be one flat vector containing sequence 1 followed
    by sequence 2. The final layer returns one raw logit per protein pair.
    """

    def __init__(self, input_dim: int, hidden_dim: int = 512, dropout: float = 0.2) -> None:
        super().__init__()

        # Keep the first model intentionally simple: Linear -> ReLU -> Linear -> ReLU -> Linear.
        self.network = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, 1),
        )

    def forward(self, features: torch.Tensor) -> torch.Tensor:
        """Return one interaction logit for each row in the batch."""
        return self.network(features).squeeze(dim=-1)
