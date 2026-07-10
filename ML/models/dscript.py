"""D-SCRIPT style model for protein interaction prediction.

The model works on precomputed per-residue embeddings instead of raw sequences.
For one protein pair the important shapes are:

    embeddings_1: [batch, residues_1, embedding_dim]
    embeddings_2: [batch, residues_2, embedding_dim]
    contact map:  [batch, residues_1, residues_2]

The contact map contains one score for every residue-residue pair across the two
proteins. The public ``forward`` method returns an interaction logit because the
training code uses ``binary_cross_entropy_with_logits`` for numerical stability.
Use ``contact_map`` or ``map_predict`` when probabilities in [0, 1] are needed.
"""

from __future__ import annotations

import torch
from torch import nn


class LogisticActivation(nn.Module):
    """Generalized sigmoid used by D-SCRIPT for the final interaction score."""

    def __init__(self, midpoint: float = 0.5, slope: float = 20.0, trainable_slope: bool = False) -> None:
        super().__init__()
        # The final D-SCRIPT score is first pooled into a value around [0, 1].
        # This activation turns that pooled value into a sharp interaction
        # probability centered around ``midpoint``. A high slope makes scores
        # above the midpoint quickly approach 1 and scores below it approach 0.
        self.midpoint = float(midpoint)
        self.slope = nn.Parameter(torch.tensor(float(slope)))
        self.slope.requires_grad = trainable_slope

    def logits(self, values: torch.Tensor) -> torch.Tensor:
        """Return pre-sigmoid logits for BCEWithLogits-style training."""
        # BCEWithLogitsLoss expects raw logits, so this method applies only the
        # affine part of the generalized sigmoid. ``forward`` applies sigmoid.
        return torch.clamp(self.slope, min=0.0) * (values - self.midpoint)

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        """Return probabilities in [0, 1]."""
        return torch.sigmoid(self.logits(values))

    def clip_parameters(self) -> None:
        """Keep the sigmoid slope non-negative when it is trainable."""
        with torch.no_grad():
            self.slope.clamp_(min=0.0)


class DScriptProjection(nn.Module):
    """Project per-residue input embeddings to D-SCRIPT's compact residue space."""

    def __init__(self, input_dim: int, projection_dim: int = 100, dropout: float = 0.5) -> None:
        super().__init__()
        # ESM-style embeddings are high dimensional. D-SCRIPT first learns a
        # compact residue representation so the pairwise contact module is much
        # cheaper and can focus on interaction-specific features.
        self.network = nn.Sequential(
            nn.Linear(input_dim, projection_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
        )

    def forward(self, embeddings: torch.Tensor) -> torch.Tensor:
        """Project embeddings shaped [batch, residues, input_dim]."""
        return self.network(embeddings)


class DScriptContactModule(nn.Module):
    """Predict inter-protein residue contact logits from projected embeddings."""

    def __init__(self, projection_dim: int = 100, hidden_dim: int = 50, contact_width: int = 7) -> None:
        super().__init__()
        if contact_width <= 0 or contact_width % 2 == 0:
            raise ValueError("contact_width must be a positive odd integer")

        # For every residue pair (i, j), D-SCRIPT builds two feature vectors:
        # abs(z_i - z_j) captures similarity/distance between residue features,
        # and z_i * z_j captures feature agreement. Concatenating both gives
        # 2 * projection_dim channels for each residue-residue pair.
        #
        # Conv2d sees the two residue axes as spatial dimensions:
        # [batch, channels, residues_1, residues_2].
        self.broadcast_conv = nn.Conv2d(2 * projection_dim, hidden_dim, kernel_size=1)
        self.broadcast_activation = nn.ReLU()
        self.broadcast_norm = nn.BatchNorm2d(hidden_dim)

        # The 1x1 convolution above mixes features at each residue pair
        # independently. This wider convolution then lets nearby residue pairs
        # influence one another, matching the local-contact idea from D-SCRIPT.
        self.contact_conv = nn.Conv2d(hidden_dim, 1, kernel_size=contact_width, padding=contact_width // 2)
        self.contact_norm = nn.BatchNorm2d(1)
        self.clip_parameters()

    def clip_parameters(self) -> None:
        """Keep the local contact filter invariant to residue-axis transposition."""
        with torch.no_grad():
            # Symmetrizing the spatial kernel makes the local filter treat the
            # two residue axes consistently. This is the constrained-filter idea
            # used in D-SCRIPT's contact prediction module.
            weights = self.contact_conv.weight
            self.contact_conv.weight.copy_(0.5 * (weights + weights.transpose(2, 3)))

    def forward(self, projected_1: torch.Tensor, projected_2: torch.Tensor) -> torch.Tensor:
        """Return contact logits shaped [batch, residues_1, residues_2]."""
        # Move the residue representation dimension into the Conv2d channel
        # position. After this, z1 is [batch, projection_dim, residues_1] and
        # z2 is [batch, projection_dim, residues_2].
        z1 = projected_1.transpose(1, 2)
        z2 = projected_2.transpose(1, 2)

        # Broadcasting creates one feature vector for every residue pair:
        # residue_diff/product are [batch, projection_dim, residues_1, residues_2].
        residue_diff = torch.abs(z1.unsqueeze(3) - z2.unsqueeze(2))
        residue_product = z1.unsqueeze(3) * z2.unsqueeze(2)
        pair_features = torch.cat((residue_diff, residue_product), dim=1)

        # First map pair features to hidden contact features, then predict one
        # raw contact logit per residue pair. Sigmoid is intentionally not used
        # here so the training loss can use BCEWithLogits directly.
        hidden = self.broadcast_conv(pair_features)
        hidden = self.broadcast_activation(hidden)
        hidden = self.broadcast_norm(hidden)

        contact_logits = self.contact_conv(hidden)
        contact_logits = self.contact_norm(contact_logits)
        return contact_logits.squeeze(dim=1)


class DScriptInteractionModel(nn.Module):
    """D-SCRIPT interaction model over two per-residue embedding tensors.

    The default ``forward`` output is a raw interaction logit to match the
    current training code's BCEWithLogits convention. Use ``contact_map`` or
    ``map_predict`` when the predicted residue contact map is needed.
    """

    def __init__(
        self,
        embedding_dim: int,
        projection_dim: int = 100,
        contact_hidden_dim: int = 50,
        contact_width: int = 7,
        projection_dropout: float = 0.5,
        use_weight_matrix: bool = True,
        theta_init: float = 1.0,
        lambda_init: float = 0.0,
        final_midpoint: float = 0.5,
        final_slope: float = 20.0,
        trainable_final_slope: bool = False,
    ) -> None:
        super().__init__()
        # The model has three conceptual parts:
        # 1. project high-dimensional residue embeddings,
        # 2. predict a residue-residue contact map,
        # 3. pool the contact map into one protein-pair interaction score.
        self.projection = DScriptProjection(embedding_dim, projection_dim, projection_dropout)
        self.contact_module = DScriptContactModule(projection_dim, contact_hidden_dim, contact_width)
        self.use_weight_matrix = use_weight_matrix

        if self.use_weight_matrix:
            # The optional D-SCRIPT weight matrix down/up-weights contact scores
            # based on position in the two proteins. Its size is created
            # dynamically in _weight_matrix, so it always matches the current
            # [residues_1, residues_2] contact map.
            self.theta = nn.Parameter(torch.tensor(float(theta_init)))
            self.lambda_ = nn.Parameter(torch.tensor(float(lambda_init)))

        # This activation converts the pooled contact evidence into a final
        # interaction logit/probability.
        self.final_activation = LogisticActivation(final_midpoint, final_slope, trainable_final_slope)
        self.clip_parameters()

    def clip_parameters(self) -> None:
        """Clamp constrained parameters to the ranges used by D-SCRIPT."""
        self.contact_module.clip_parameters()
        self.final_activation.clip_parameters()
        with torch.no_grad():
            if self.use_weight_matrix:
                self.theta.clamp_(min=0.0, max=1.0)
                self.lambda_.clamp_(min=0.0)

    def contact_logits(
        self,
        embeddings_1: torch.Tensor,
        embeddings_2: torch.Tensor,
        mask_1: torch.Tensor | None = None,
        mask_2: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """Return pre-sigmoid residue contact logits."""
        self._validate_inputs(embeddings_1, embeddings_2)

        # The same projection network is applied to both proteins so both
        # residue sets live in the same learned D-SCRIPT feature space.
        projected_1 = self.projection(embeddings_1)
        projected_2 = self.projection(embeddings_2)
        logits = self.contact_module(projected_1, projected_2)

        # Padding residues should never contribute to contact loss or
        # interaction pooling. Masked logits are set to zero; downstream code
        # always carries the mask along when the actual values matter.
        pair_mask = self._pair_mask(logits, mask_1, mask_2)
        return logits.masked_fill(~pair_mask, 0.0)

    def contact_map(
        self,
        embeddings_1: torch.Tensor,
        embeddings_2: torch.Tensor,
        mask_1: torch.Tensor | None = None,
        mask_2: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """Return the predicted residue contact map with values in [0, 1]."""
        logits = self.contact_logits(embeddings_1, embeddings_2, mask_1, mask_2)
        # Sigmoid turns contact logits into probabilities. These are the values
        # to threshold at e.g. 0.5 when reporting contact-map metrics.
        contacts = torch.sigmoid(logits)
        pair_mask = self._pair_mask(contacts, mask_1, mask_2)
        return contacts.masked_fill(~pair_mask, 0.0)

    def map_predict(
        self,
        embeddings_1: torch.Tensor,
        embeddings_2: torch.Tensor,
        mask_1: torch.Tensor | None = None,
        mask_2: torch.Tensor | None = None,
        interaction_pair_mask: torch.Tensor | None = None,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """Return ``(contact_map, interaction_probability)``."""
        logits, contacts = self.forward(
            embeddings_1,
            embeddings_2,
            mask_1,
            mask_2,
            interaction_pair_mask=interaction_pair_mask,
            return_contact_map=True,
        )
        return contacts, torch.sigmoid(logits)

    def predict_probability(
        self,
        embeddings_1: torch.Tensor,
        embeddings_2: torch.Tensor,
        mask_1: torch.Tensor | None = None,
        mask_2: torch.Tensor | None = None,
        interaction_pair_mask: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """Return interaction probabilities in [0, 1]."""
        return torch.sigmoid(
            self.forward(
                embeddings_1,
                embeddings_2,
                mask_1,
                mask_2,
                interaction_pair_mask=interaction_pair_mask,
            )
        )

    def forward(
        self,
        embeddings_1: torch.Tensor,
        embeddings_2: torch.Tensor,
        mask_1: torch.Tensor | None = None,
        mask_2: torch.Tensor | None = None,
        return_contact_map: bool = False,
        return_contact_logits: bool = False,
        interaction_pair_mask: torch.Tensor | None = None,
    ) -> torch.Tensor | tuple[torch.Tensor, torch.Tensor] | tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Return interaction logits, optionally with contact probabilities/logits."""
        # Contact logits are needed for BCEWithLogits contact-map supervision.
        # Contact probabilities are needed for D-SCRIPT-style pooling and for
        # human-interpretable/thresholded metrics.
        contact_logits = self.contact_logits(embeddings_1, embeddings_2, mask_1, mask_2)
        contacts = torch.sigmoid(contact_logits)
        pair_mask = self._pair_mask(contacts, mask_1, mask_2)
        contacts = contacts.masked_fill(~pair_mask, 0.0)

        # Pool only residue pairs allowed for interaction prediction. The base
        # pair_mask removes batch padding; interaction_pair_mask can further
        # remove unresolved/unknown structure cells from true contact maps.
        pooling_mask = self._pooling_pair_mask(pair_mask, interaction_pair_mask)
        logits = self._interaction_logits(contacts, pooling_mask)

        if return_contact_map and return_contact_logits:
            return logits, contacts, contact_logits
        if return_contact_logits:
            return logits, contact_logits
        if return_contact_map:
            return logits, contacts
        return logits

    def _interaction_logits(self, contacts: torch.Tensor, pair_mask: torch.Tensor) -> torch.Tensor:
        """Pool contact probabilities into one raw interaction logit per pair."""
        scores = contacts
        if self.use_weight_matrix:
            # The positional weight matrix has shape [residues_1, residues_2].
            # It is unsqueezed so the same matrix broadcasts across the batch.
            scores = scores * self._weight_matrix(contacts).unsqueeze(0)
        scores = scores.masked_fill(~pair_mask, 0.0)

        # D-SCRIPT does not average all contacts directly. It first computes the
        # mean score over valid residue pairs and then focuses on contacts above
        # that mean, so weak background contact probabilities contribute less.
        valid = pair_mask.to(dtype=scores.dtype)
        valid_count = valid.sum(dim=(1, 2)).clamp_min(1.0)
        mean = (scores * valid).sum(dim=(1, 2)) / valid_count

        # ``high_contacts`` keeps only contact evidence above the per-example
        # mean. The +1 in the denominator is the same stabilizing convention used
        # by D-SCRIPT to avoid division by zero when no score exceeds the mean.
        high_contacts = torch.relu(scores - mean[:, None, None]) * valid
        high_contact_count = (high_contacts > 0).to(dtype=scores.dtype).sum(dim=(1, 2))
        raw_probability = high_contacts.sum(dim=(1, 2)) / (high_contact_count + 1.0)
        return self.final_activation.logits(raw_probability)

    def _weight_matrix(self, contacts: torch.Tensor) -> torch.Tensor:
        """Create the D-SCRIPT positional weight matrix for the current map size."""
        # The map size comes from the actual contact tensor, so no fixed maximum
        # sequence length is baked into the model. This is why variable-length
        # batches can still use the correct [residues_1, residues_2] weights.
        residues_1, residues_2 = contacts.shape[1:]
        device = contacts.device
        dtype = contacts.dtype

        # Positions are 1-based to follow the D-SCRIPT formula.
        positions_1 = torch.arange(1, residues_1 + 1, device=device, dtype=dtype)
        positions_2 = torch.arange(1, residues_2 + 1, device=device, dtype=dtype)

        # The Gaussian-like weights are centered in each sequence. lambda_
        # controls how strongly position matters, and theta mixes the learned
        # positional weighting with a uniform all-ones matrix.
        center_1 = (residues_1 + 1) / 2.0
        center_2 = (residues_2 + 1) / 2.0
        scale_1 = -((residues_1 + 1) / 2.0)
        scale_2 = -((residues_2 + 1) / 2.0)

        lambda_ = torch.clamp(self.lambda_, min=0.0)
        weights_1 = torch.exp(lambda_ * -torch.square((positions_1 - center_1) / scale_1))
        weights_2 = torch.exp(lambda_ * -torch.square((positions_2 - center_2) / scale_2))

        theta = torch.clamp(self.theta, min=0.0, max=1.0)
        weights = weights_1[:, None] * weights_2[None, :]
        return (1.0 - theta) * weights + theta

    @staticmethod
    def _validate_inputs(embeddings_1: torch.Tensor, embeddings_2: torch.Tensor) -> None:
        """Check that the two protein tensors can form a batched pair map."""
        if embeddings_1.ndim != 3 or embeddings_2.ndim != 3:
            raise ValueError("D-SCRIPT expects tensors shaped [batch, residues, embedding_dim]")
        if embeddings_1.shape[0] != embeddings_2.shape[0]:
            raise ValueError("Both embedding tensors must have the same batch size")
        if embeddings_1.shape[2] != embeddings_2.shape[2]:
            raise ValueError("Both embedding tensors must have the same embedding dimension")

    @staticmethod
    def _pair_mask(
        contacts: torch.Tensor,
        mask_1: torch.Tensor | None,
        mask_2: torch.Tensor | None,
    ) -> torch.Tensor:
        """Return a [batch, residues_1, residues_2] mask for valid residue pairs."""
        batch_size, residues_1, residues_2 = contacts.shape
        device = contacts.device

        # When no masks are passed, every residue in the tensor is treated as
        # real. The training pipeline passes masks because it pads variable-
        # length proteins to the largest protein length in the batch.
        if mask_1 is None:
            mask_1 = torch.ones((batch_size, residues_1), dtype=torch.bool, device=device)
        else:
            mask_1 = mask_1.to(device=device, dtype=torch.bool)

        if mask_2 is None:
            mask_2 = torch.ones((batch_size, residues_2), dtype=torch.bool, device=device)
        else:
            mask_2 = mask_2.to(device=device, dtype=torch.bool)

        if mask_1.shape != (batch_size, residues_1):
            raise ValueError("mask_1 must be shaped [batch, residues_1]")
        if mask_2.shape != (batch_size, residues_2):
            raise ValueError("mask_2 must be shaped [batch, residues_2]")

        # Combine the two 1D residue masks into a 2D residue-pair mask. A pair is
        # valid only if both the residue from protein 1 and the residue from
        # protein 2 are real, non-padding residues.
        return mask_1[:, :, None] & mask_2[:, None, :]

    @staticmethod
    def _pooling_pair_mask(
        pair_mask: torch.Tensor,
        interaction_pair_mask: torch.Tensor | None,
    ) -> torch.Tensor:
        """Return the residue-pair mask used for interaction-score pooling."""
        if interaction_pair_mask is None:
            return pair_mask

        # interaction_pair_mask has the same [batch, residues_1, residues_2]
        # layout as the contact map. True means the residue pair is allowed to
        # influence the interaction score. False removes padding or unresolved
        # structure cells such as contact-map targets encoded as -1.
        interaction_pair_mask = interaction_pair_mask.to(device=pair_mask.device, dtype=torch.bool)
        if interaction_pair_mask.shape != pair_mask.shape:
            raise ValueError("interaction_pair_mask must be shaped [batch, residues_1, residues_2]")
        return pair_mask & interaction_pair_mask
