"""Transformer-query contact-map model with D-SCRIPT-compatible outputs.

Data flow for one protein pair:

    cached residue embeddings
        -> shared D-SCRIPT projection
        -> learned queries attend to both projected proteins
        -> each query predicts one mask on each protein and a patch strength
        -> matched query patches form residue-pair query evidence
        -> gate compatibility probabilities (default), or add both logit evidences
        -> the existing D-SCRIPT normal or max interaction head pools that map

The transformer never receives structural annotation masks. ``mask_1`` and
``mask_2`` identify real sequence residues versus batch padding only. Unknown
contact-map cells (target value -1) remain a loss/metric concern in the training
pipeline and must not influence attention or interaction pooling.

``query_evidence_mode="additive"`` restores QP-Additive: query evidence plus
compatibility logits, including the contact bias once. Its training wrapper uses
only the main contact/interaction losses; gated mode also has branch supervision.
"""

from __future__ import annotations

import math

import torch
import torch.nn.functional as F
from torch import nn

from .dscript import DScriptContactModule, DScriptInteractionModel


class QueryPatchContactModule(nn.Module):
    """Construct contact logits from learned, partner-conditioned patch queries.

    Inputs are projected residue tensors ``[B, L1, D]`` and ``[B, L2, D]``.
    The output is a single dense contact-logit tensor ``[B, L1, L2]``. The
    implementation contracts query masks with ``torch.bmm`` and therefore does
    not materialize a large ``[B, K, L1, L2]`` tensor.
    """

    def __init__(
        self,
        projection_dim: int = 100,
        num_queries: int = 32,
        num_heads: int = 4,
        num_layers: int = 1,
        dropout: float = 0.1,
        contact_bias_init: float = -6.0,
        compatibility_rank: int = 32,
        compatibility_scale: float = 2,
        query_gate_bias_init: float = 0.0,
        compatibility_head: str = "bilinear",
        contact_hidden_dim: int = 50,
        contact_width: int = 7,
        query_evidence_mode: str = "gated",
    ) -> None:
        super().__init__()
        if projection_dim <= 0:
            raise ValueError("projection_dim must be positive")
        if num_queries <= 0:
            raise ValueError("num_queries must be positive")
        if num_heads <= 0 or projection_dim % num_heads != 0:
            raise ValueError("num_heads must be positive and divide projection_dim")
        if num_layers <= 0:
            raise ValueError("num_layers must be positive")
        if not 0.0 <= dropout < 1.0:
            raise ValueError("dropout must be in [0, 1)")
        if compatibility_head not in {"bilinear", "dscript"}:
            raise ValueError(
                "compatibility_head must be either 'bilinear' or 'dscript'"
            )
        if query_evidence_mode not in {"gated", "additive"}:
            raise ValueError("query_evidence_mode must be either 'gated' or 'additive'")

        self.projection_dim = projection_dim
        self.num_queries = num_queries
        self.query_evidence_mode = query_evidence_mode

        # The K vectors are model parameters shared by all protein pairs. They
        # start as distinct random slots and become generic interface-search
        # prompts through gradients from contact and interaction objectives.
        self.queries = nn.Parameter(torch.empty(num_queries, projection_dim))
        nn.init.normal_(self.queries, mean=0.0, std=0.02)

        # One decoder is enough for the initial model. Its query self-attention
        # lets slots coordinate, while cross-attention lets every slot summarize
        # both proteins. No chain-role embedding is added, so exchanging protein
        # A and B only reorders the memory and preserves swap equivariance.
        decoder_layer = nn.TransformerDecoderLayer(
            d_model=projection_dim,
            nhead=num_heads,
            dim_feedforward=4 * projection_dim,
            dropout=dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.decoder = nn.TransformerDecoder(
            decoder_layer,
            num_layers=num_layers,
            norm=nn.LayerNorm(projection_dim),
        )
        self.memory_norm = nn.LayerNorm(projection_dim)

        # Query and residue states play different roles in mask prediction, so
        # they receive separate D->D transformations. The residue transformation
        # is shared across both proteins. Identity initialization starts the head
        # as a normalized dot product while still allowing a learned bilinear
        # compatibility metric during training.
        self.query_mask_norm = nn.LayerNorm(projection_dim)
        self.residue_mask_norm = nn.LayerNorm(projection_dim)
        self.query_mask_projection = nn.Linear(projection_dim, projection_dim, bias=False)
        self.residue_mask_projection = nn.Linear(projection_dim, projection_dim, bias=False)
        nn.init.eye_(self.query_mask_projection.weight)
        nn.init.eye_(self.residue_mask_projection.weight)

        # Each query contributes a non-negative patch. A negative initial
        # strength keeps those patches modest. The sparse background contact
        # bias belongs to compatibility logits and is used once in either mode.
        self.patch_strength = nn.Linear(projection_dim, 1)
        nn.init.constant_(self.patch_strength.bias, -2.0)
        self.contact_bias = nn.Parameter(torch.tensor(float(contact_bias_init)))

        # Only gated mode learns a gate scale/bias. Additive mode uses raw query
        # evidence, as in QP-Additive, and has no inactive gate parameters in its
        # optimizer or state dict. Gated checkpoint parameter names stay intact.
        if query_evidence_mode == "gated":
            initial_raw_scale = math.log(math.expm1(1.0))
            self.query_gate_raw_scale = nn.Parameter(torch.tensor(initial_raw_scale))
            self.query_gate_bias = nn.Parameter(torch.tensor(float(query_gate_bias_init)))
        else:
            self.register_parameter("query_gate_raw_scale", None)
            self.register_parameter("query_gate_bias", None)

        # Construct only the selected compatibility head. This keeps old
        # bilinear checkpoints unchanged and avoids adding inactive parameters
        # to the optimizer when the local D-SCRIPT head is selected.
        self.compatibility_head = compatibility_head
        self.compatibility_module: DScriptContactModule | None = None
        if compatibility_head == "bilinear":
            self.compatibility_rank = compatibility_rank
            self.compatibility_norm = nn.LayerNorm(projection_dim)
            self.compatibility_projection = nn.Linear(
                projection_dim,
                compatibility_rank,
                bias=False,
            )

            # Alternating signs let the low-rank metric learn both matching and
            # opposing residue features while retaining the old initialization.
            metric = torch.ones(compatibility_rank)
            metric[1::2] = -1.0
            self.compatibility_metric = nn.Parameter(metric)
            self.register_buffer(
                "compatibility_scale",
                torch.tensor(float(compatibility_scale)),
            )
        else:
            # This head builds |z1-z2| and z1*z2 pair features, then refines
            # them with a local convolution over the two residue axes.
            self.compatibility_module = DScriptContactModule(
                projection_dim=projection_dim,
                hidden_dim=contact_hidden_dim,
                contact_width=contact_width,
            )

    def forward(
        self,
        projected_1: torch.Tensor,
        projected_2: torch.Tensor,
        mask_1: torch.Tensor | None = None,
        mask_2: torch.Tensor | None = None,
        return_patch_details: bool = False,
    ) -> (
        torch.Tensor
        | tuple[
            torch.Tensor,
            torch.Tensor,
            torch.Tensor,
            torch.Tensor,
            torch.Tensor,
            torch.Tensor | None,
            torch.Tensor,
            torch.Tensor,
            torch.Tensor,
        ]
    ):
        """Return contact logits, optionally with both masks and strengths.

        ``mask_1`` and ``mask_2`` are boolean sequence-padding masks with shapes
        ``[B, L1]`` and ``[B, L2]``. ``True`` means a real residue. When patch
        details are requested, the additional outputs are the two query masks
        ``[B,K,L1]``/``[B,K,L2]``, strengths ``[B,K]``, query evidence
        ``[B,L1,L2]``, gate logits ``[B,L1,L2]``, compatibility logits
        ``[B,L1,L2]``, and chain-wise interface logits ``[B,L1]``/``[B,L2]``.
        The gate-logit entry is ``None`` in additive mode because no gate exists.
        """
        if projected_1.ndim != 3 or projected_2.ndim != 3:
            raise ValueError("projected residue tensors must be shaped [batch, residues, projection_dim]")
        if projected_1.shape[0] != projected_2.shape[0]:
            raise ValueError("projected residue tensors must have the same batch size")
        if projected_1.shape[2] != self.projection_dim or projected_2.shape[2] != self.projection_dim:
            raise ValueError("projected residue tensors do not match projection_dim")

        batch_size, residues_1, _ = projected_1.shape
        residues_2 = projected_2.shape[1]
        device = projected_1.device

        # Standalone, unpadded calls may omit masks. Batches from the training
        # collator provide them explicitly so attention never reads zero padding.
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
        if not mask_1.any(dim=1).all() or not mask_2.any(dim=1).all():
            raise ValueError("every protein must contain at least one unmasked residue")

        # Concatenation is only a transformer-memory layout; the two proteins
        # remain separate for mask decoding below. Cross-attention is invariant
        # to memory ordering in evaluation mode, which supports A/B swap symmetry.
        memory = self.memory_norm(torch.cat((projected_1, projected_2), dim=1))
        memory_mask = torch.cat((mask_1, mask_2), dim=1)
        queries = self.queries.unsqueeze(0).expand(batch_size, -1, -1)
        query_states = self.decoder(
            tgt=queries,
            memory=memory,
            memory_key_padding_mask=~memory_mask,
        )
        # query_states: [B, K, D], conditioned jointly on both proteins.

        # A shared dot-product mask head compares every joint query state with
        # every residue in each protein. The sigmoid outputs are soft patch
        # memberships rather than normalized attention probabilities.
        query_features = self.query_mask_projection(self.query_mask_norm(query_states))
        residue_features_1 = self.residue_mask_projection(self.residue_mask_norm(projected_1))
        residue_features_2 = self.residue_mask_projection(self.residue_mask_norm(projected_2))
        scale = math.sqrt(self.projection_dim)
        mask_logits_1 = torch.bmm(query_features, residue_features_1.transpose(1, 2)) / scale
        mask_logits_2 = torch.bmm(query_features, residue_features_2.transpose(1, 2)) / scale
        # shape: [B, K, L]

        # A residue belongs to the predicted interface when at least one query
        # selects it. Max aggregation lets queries specialize: positive residues
        # need one matching slot, while a negative residue forces every slot low.
        interface_logits_1 = mask_logits_1.amax(dim=1).masked_fill(~mask_1, -20.0)
        interface_logits_2 = mask_logits_2.amax(dim=1).masked_fill(~mask_2, -20.0)

        patch_masks_1 = torch.sigmoid(mask_logits_1).masked_fill(~mask_1[:, None, :], 0.0)
        patch_masks_2 = torch.sigmoid(mask_logits_2).masked_fill(~mask_2[:, None, :], 0.0)
        patch_strengths = F.softplus(self.patch_strength(query_states).squeeze(-1))

        # A^T diag(strength) B composes K rank-one patches into one [B,L1,L2]
        # evidence tensor. Scaling keeps its initial magnitude comparable
        # when num_queries changes.
        weighted_masks_1 = patch_masks_1 * patch_strengths[:, :, None]
        query_evidence = torch.bmm(weighted_masks_1.transpose(1, 2), patch_masks_2)
        query_evidence = query_evidence / math.sqrt(self.num_queries)
        # Function: A^T * diag(strength) * B / sqrt(K)

        pair_mask = mask_1[:, :, None] & mask_2[:, None, :]

        query_evidence = query_evidence.masked_fill(~pair_mask, 0.0)

        # Both heads return one compatibility logit per residue pair, including
        # the background bias. Query evidence is combined with it afterward.
        if self.compatibility_head == "bilinear":
            normalized_1 = self.compatibility_norm(projected_1)
            normalized_2 = self.compatibility_norm(projected_2)
            # shape: [B, Lx, D]

            compatibility_1 = F.normalize(
                self.compatibility_projection(normalized_1),
                p=2,
                dim=-1,
            )
            compatibility_2 = F.normalize(
                self.compatibility_projection(normalized_2),
                p=2,
                dim=-1,
            )
            # shape: [B, Lx, compatibility_rank]

            compatibility_metric = torch.tanh(self.compatibility_metric)
            compatibility_evidence = torch.bmm(
                compatibility_1 * compatibility_metric,
                compatibility_2.transpose(1, 2),
            )
            compatibility_evidence = compatibility_evidence.masked_fill(
                ~pair_mask,
                0.0,
            )
            compatibility_contact_logits = (
                self.compatibility_scale * compatibility_evidence
                + self.contact_bias
            )
        else:
            if self.compatibility_module is None:
                raise RuntimeError("D-SCRIPT compatibility module is not initialized")
            compatibility_contact_logits = (
                self.compatibility_module(projected_1, projected_2)
                + self.contact_bias
            )

        # =====================================================================
        # Evidence combination: historical additive logits or gated probability
        # =====================================================================
        if self.query_evidence_mode == "additive":
            # QP-Additive: E_query + alpha * E_compatibility + contact_bias.
            # Compatibility logits already contain the scale and bias. Query
            # evidence receives no learned gate scale, bias, or extra sigmoid.
            contact_logits = query_evidence + compatibility_contact_logits
            query_gate_logits = None
        else:
            # Query evidence is non-negative; a positive learned scale preserves
            # the interpretation that stronger matched patches open the gate.
            gate_scale = F.softplus(self.query_gate_raw_scale)
            query_gate_logits = gate_scale * query_evidence + self.query_gate_bias

            # p(contact) = sigmoid(gate_logits) * sigmoid(compatibility_logits).
            # Log-add-exp converts the product back to stable BCE logits and
            # ensures contact probability cannot exceed gate probability.
            contact_logits = -torch.logaddexp(
                -query_gate_logits,
                torch.logaddexp(
                    -compatibility_contact_logits,
                    -query_gate_logits - compatibility_contact_logits,
                ),
            )

        if return_patch_details:
            return (
                contact_logits,
                patch_masks_1,
                patch_masks_2,
                patch_strengths,
                query_evidence,
                query_gate_logits,
                compatibility_contact_logits,
                interface_logits_1,
                interface_logits_2,
            )
        return contact_logits

    def clip_parameters(self) -> None:
        """Apply the selected compatibility head's parameter constraints."""
        if self.compatibility_module is not None:
            self.compatibility_module.clip_parameters()


class QueryPatchInteractionModel(DScriptInteractionModel):
    """D-SCRIPT projection and interaction pooling with a query-patch map head."""

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
        interaction_module_type: str = "normal",
        gamma_init: float = 0.0,
        max_pooling: bool = True,
        num_queries: int = 32,
        query_heads: int = 4,
        query_layers: int = 1,
        query_dropout: float = 0.1,
        contact_bias_init: float = -6.0,
        compatibility_rank: int = 32,
        compatibility_scale: float = 2,
        query_gate_bias_init: float = 0.0,
        compatibility_head: str = "bilinear",
        query_evidence_mode: str = "gated",
    ) -> None:
        # Constructing the common parent first retains its projection, public
        # prediction API, constrained interaction parameters, and both pooling
        # heads. Only the pair-feature/CNN contact module is replaced afterward.
        super().__init__(
            embedding_dim=embedding_dim,
            projection_dim=projection_dim,
            contact_hidden_dim=contact_hidden_dim,
            contact_width=contact_width,
            projection_dropout=projection_dropout,
            use_weight_matrix=use_weight_matrix,
            theta_init=theta_init,
            lambda_init=lambda_init,
            final_midpoint=final_midpoint,
            final_slope=final_slope,
            trainable_final_slope=trainable_final_slope,
            interaction_module_type=interaction_module_type,
            gamma_init=gamma_init,
            max_pooling=max_pooling,
        )
        self.contact_module = QueryPatchContactModule(
            projection_dim=projection_dim,
            num_queries=num_queries,
            num_heads=query_heads,
            num_layers=query_layers,
            dropout=query_dropout,
            contact_bias_init=contact_bias_init,
            compatibility_rank=compatibility_rank,
            compatibility_scale=compatibility_scale,
            query_gate_bias_init=query_gate_bias_init,
            compatibility_head=compatibility_head,
            contact_hidden_dim=contact_hidden_dim,
            contact_width=contact_width,
            query_evidence_mode=query_evidence_mode,
        )

    def contact_logits(
        self,
        embeddings_1: torch.Tensor,
        embeddings_2: torch.Tensor,
        mask_1: torch.Tensor | None = None,
        mask_2: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """Return query-patch contact logits through the inherited model API.

        D-SCRIPT's original contact module does not consume padding masks. This
        model-specific override passes them only to the Transformer contact
        module, keeping the baseline implementation and behavior unchanged.
        """
        contact_logits, _, _, _, _, _, _, _, _ = self.query_patch_outputs(
            embeddings_1,
            embeddings_2,
            mask_1,
            mask_2,
        )
        return contact_logits

    def query_patch_outputs(
        self,
        embeddings_1: torch.Tensor,
        embeddings_2: torch.Tensor,
        mask_1: torch.Tensor | None = None,
        mask_2: torch.Tensor | None = None,
    ) -> tuple[
        torch.Tensor,
        torch.Tensor,
        torch.Tensor,
        torch.Tensor,
        torch.Tensor,
        torch.Tensor | None,
        torch.Tensor,
        torch.Tensor,
        torch.Tensor,
    ]:
        """Return contact, query-patch, gate, compatibility, and interface tensors.

        This diagnostic method exposes which residues each learned query selects
        without changing the standard D-SCRIPT-compatible ``forward`` contract.
        Contact logits at padded residue pairs are zeroed exactly as in the base
        model; returned query masks are already zero at padded residues.
        Additive mode returns ``None`` in the gate-logit slot.
        """
        self._validate_inputs(embeddings_1, embeddings_2)
        projected_1 = self.projection(embeddings_1)
        projected_2 = self.projection(embeddings_2)
        outputs = self.contact_module(
            projected_1,
            projected_2,
            mask_1,
            mask_2,
            return_patch_details=True,
        )
        (
            contact_logits,
            patch_masks_1,
            patch_masks_2,
            patch_strengths,
            query_evidence,
            query_gate_logits,
            compatibility_contact_logits,
            interface_logits_1,
            interface_logits_2,
        ) = outputs
        pair_mask = self._pair_mask(contact_logits, mask_1, mask_2)
        contact_logits = contact_logits.masked_fill(~pair_mask, 0.0)
        if query_gate_logits is not None:
            query_gate_logits = query_gate_logits.masked_fill(~pair_mask, -20.0)
        compatibility_contact_logits = compatibility_contact_logits.masked_fill(
            ~pair_mask, -20.0
        )
        return (
            contact_logits,
            patch_masks_1,
            patch_masks_2,
            patch_strengths,
            query_evidence,
            query_gate_logits,
            compatibility_contact_logits,
            interface_logits_1,
            interface_logits_2,
        )

    def _contact_logits_with_aux(
        self,
        embeddings_1: torch.Tensor,
        embeddings_2: torch.Tensor,
        mask_1: torch.Tensor | None = None,
        mask_2: torch.Tensor | None = None,
    ) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
        """Expose branch-specific tensors without a second decoder pass."""
        (
            contact_logits,
            _,
            _,
            _,
            _,
            query_gate_logits,
            compatibility_contact_logits,
            interface_logits_1,
            interface_logits_2,
        ) = self.query_patch_outputs(embeddings_1, embeddings_2, mask_1, mask_2)
        # QP-Additive predates branch-specific auxiliary supervision. Its query
        # and compatibility branches are trained through the final logits only.
        if self.contact_module.query_evidence_mode == "additive":
            return contact_logits, {}
        return contact_logits, {
            "query_gate_logits": query_gate_logits,
            "query_gate_scale": F.softplus(self.contact_module.query_gate_raw_scale),
            "query_gate_bias": self.contact_module.query_gate_bias,
            "compatibility_contact_logits": compatibility_contact_logits,
            "interface_logits_1": interface_logits_1,
            "interface_logits_2": interface_logits_2,
        }
