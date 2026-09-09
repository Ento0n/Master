"""Regression checks for historical additive and current gated QUERY-PATCH.

Run from the repository root with ``python -m unittest discover -s tests -v``.
Small CPU tensors exercise the real model, training objective, and checkpoint
loader; no protein dataset or training run is required.
"""

import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import torch
import torch.nn.functional as F

from scripts.models.query_patch import QueryPatchContactModule, QueryPatchInteractionModel


# =============================================================================
# Contact evidence, padding, symmetry, and differentiable interaction heads
# =============================================================================


class QueryPatchModelTests(unittest.TestCase):
    def setUp(self):
        torch.manual_seed(12)
        # [batch=1, residues, projected/embedding dimensions=8]. Unequal chain
        # lengths ensure transposing a swapped pair is meaningfully exercised.
        self.protein_1 = torch.randn(1, 5, 8)
        self.protein_2 = torch.randn(1, 7, 8)
        self.contact_options = dict(
            projection_dim=8,
            num_queries=3,
            num_heads=2,
            dropout=0.0,
            contact_bias_init=-1.0,
            compatibility_rank=4,
            compatibility_scale=2.5,
            contact_hidden_dim=4,
            contact_width=3,
        )
        self.model_options = dict(
            embedding_dim=8,
            projection_dim=8,
            num_queries=3,
            query_heads=2,
            query_dropout=0.0,
            projection_dropout=0.0,
            contact_bias_init=-1.0,
            compatibility_rank=4,
            contact_hidden_dim=4,
            contact_width=3,
        )

    def test_additive_combines_query_and_compatibility_with_one_bias(self):
        for head in ("bilinear", "dscript"):
            with self.subTest(head=head):
                module = QueryPatchContactModule(
                    **self.contact_options,
                    compatibility_head=head,
                    query_evidence_mode="additive",
                ).eval()
                outputs = module(self.protein_1, self.protein_2, return_patch_details=True)
                self.assertEqual(len(outputs), 9)
                logits, masks_1, masks_2, strengths, evidence, gate, compatibility, _, _ = outputs

                # Reconstruct the historical equation from its constituent
                # masks and normalized compatibility features, including the
                # fixed bilinear scale and exactly one learned contact bias.
                expected_query = torch.einsum("bki,bk,bkj->bij", masks_1, strengths, masks_2)
                expected_query = expected_query / math.sqrt(module.num_queries)
                if head == "bilinear":
                    features_1 = F.normalize(module.compatibility_projection(module.compatibility_norm(self.protein_1)), dim=-1)
                    features_2 = F.normalize(module.compatibility_projection(module.compatibility_norm(self.protein_2)), dim=-1)
                    expected_compatibility = torch.einsum(
                        "bir,r,bjr->bij", features_1, module.compatibility_metric.tanh(), features_2
                    ) * self.contact_options["compatibility_scale"]
                else:
                    expected_compatibility = module.compatibility_module(self.protein_1, self.protein_2)
                expected_compatibility = expected_compatibility + module.contact_bias
                torch.testing.assert_close(evidence, expected_query)
                torch.testing.assert_close(compatibility, expected_compatibility)
                torch.testing.assert_close(logits, expected_query + expected_compatibility)
                self.assertIsNone(gate)
                self.assertIsNone(module.query_gate_raw_scale)
                self.assertIsNone(module.query_gate_bias)
                self.assertFalse(any("query_gate" in name for name in module.state_dict()))

    def test_default_mode_preserves_gated_parameters_and_probability_product(self):
        for head in ("bilinear", "dscript"):
            with self.subTest(head=head):
                torch.manual_seed(25)
                default = QueryPatchContactModule(**self.contact_options, compatibility_head=head).eval()
                torch.manual_seed(25)
                explicit = QueryPatchContactModule(
                    **self.contact_options, compatibility_head=head, query_evidence_mode="gated"
                ).eval()
                self.assertEqual(default.state_dict().keys(), explicit.state_dict().keys())
                for name, value in default.state_dict().items():
                    torch.testing.assert_close(value, explicit.state_dict()[name], rtol=0, atol=0)
                outputs = default(self.protein_1, self.protein_2, return_patch_details=True)
                torch.testing.assert_close(outputs[0], explicit(self.protein_1, self.protein_2))
                probability = outputs[0].sigmoid()
                gate_probability = outputs[5].sigmoid()
                torch.testing.assert_close(probability, gate_probability * outputs[6].sigmoid())
                self.assertTrue(torch.all(probability <= gate_probability).item())

    def test_zero_bilinear_scale_removes_compatibility_influence(self):
        options = {**self.contact_options, "compatibility_scale": 0.0}
        module = QueryPatchContactModule(**options, query_evidence_mode="additive").eval()
        before = module(self.protein_1, self.protein_2, return_patch_details=True)
        with torch.no_grad():
            module.compatibility_projection.weight.normal_(mean=5.0, std=2.0)
            module.compatibility_metric.neg_()
        after = module(self.protein_1, self.protein_2)
        torch.testing.assert_close(before[0], before[4] + module.contact_bias)
        torch.testing.assert_close(before[0], after, rtol=0, atol=0)

    def test_contact_maps_transpose_when_proteins_are_swapped(self):
        for mode in ("gated", "additive"):
            for head in ("bilinear", "dscript"):
                with self.subTest(mode=mode, head=head):
                    module = QueryPatchContactModule(
                        **self.contact_options, query_evidence_mode=mode, compatibility_head=head
                    ).eval()
                    original = module(self.protein_1, self.protein_2)
                    swapped = module(self.protein_2, self.protein_1)
                    torch.testing.assert_close(original, swapped.transpose(1, 2), atol=2e-6, rtol=1e-5)

    def test_padding_does_not_change_query_evidence_or_bilinear_contacts(self):
        # Padded values deliberately have large magnitude: masking must exclude
        # their content, not merely work because a collator happened to use zero.
        padded_1 = F.pad(self.protein_1, (0, 0, 0, 2), value=100.0)
        padded_2 = F.pad(self.protein_2, (0, 0, 0, 3), value=-100.0)
        mask_1 = torch.arange(7)[None, :] < 5
        mask_2 = torch.arange(10)[None, :] < 7
        pair_mask = mask_1[:, :, None] & mask_2[:, None, :]
        for mode in ("gated", "additive"):
            for head in ("bilinear", "dscript"):
                with self.subTest(mode=mode, head=head):
                    model = QueryPatchInteractionModel(
                        **self.model_options, query_evidence_mode=mode, compatibility_head=head
                    ).eval()
                    original = model.query_patch_outputs(self.protein_1, self.protein_2)
                    padded = model.query_patch_outputs(padded_1, padded_2, mask_1, mask_2)
                    torch.testing.assert_close(original[4], padded[4][:, :5, :7], atol=2e-6, rtol=1e-5)
                    self.assertTrue(torch.all(padded[1][:, :, 5:] == 0).item())
                    self.assertTrue(torch.all(padded[2][:, :, 7:] == 0).item())
                    self.assertTrue(torch.all(padded[0][~pair_mask] == 0).item())
                    # The pre-existing local D-SCRIPT CNN can mix neighboring
                    # padding into boundary cells; only its query branch has
                    # this invariance. The bilinear head is pointwise in pairs.
                    if head == "bilinear":
                        torch.testing.assert_close(original[0], padded[0][:, :5, :7], atol=2e-6, rtol=1e-5)

    def test_both_interaction_heads_train_both_evidence_branches(self):
        for mode in ("gated", "additive"):
            for head in ("bilinear", "dscript"):
                for interaction in ("normal", "max"):
                    with self.subTest(mode=mode, head=head, interaction=interaction):
                        model = QueryPatchInteractionModel(
                            **self.model_options,
                            query_evidence_mode=mode,
                            compatibility_head=head,
                            interaction_module_type=interaction,
                        ).eval()
                        logits, contacts, contact_logits, auxiliary = model(
                            self.protein_1, self.protein_2,
                            return_contact_map=True, return_contact_logits=True, return_contact_aux=True,
                        )
                        self.assertEqual(logits.shape, (1,))
                        self.assertTrue(torch.isfinite(logits).all().item())
                        torch.testing.assert_close(contacts, contact_logits.sigmoid())
                        self.assertEqual(bool(auxiliary), mode == "gated")

                        # A mixed pair/contact objective must retain gradients
                        # to learned queries and the selected compatibility head.
                        targets = torch.zeros_like(contact_logits)
                        targets[:, 1, 2] = 1.0
                        loss = F.binary_cross_entropy_with_logits(logits, torch.ones_like(logits))
                        loss = loss + F.binary_cross_entropy_with_logits(contact_logits, targets)
                        loss.backward()
                        module = model.contact_module
                        compatibility_parameter = (
                            module.compatibility_projection.weight if head == "bilinear"
                            else module.compatibility_module.contact_conv.weight
                        )
                        for parameter in (module.queries, module.patch_strength.weight, compatibility_parameter):
                            self.assertIsNotNone(parameter.grad)
                            self.assertTrue(torch.isfinite(parameter.grad).all().item())
                            self.assertGreater(parameter.grad.abs().sum().item(), 0.0)

    def test_invalid_evidence_mode_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "query_evidence_mode"):
            QueryPatchContactModule(**self.contact_options, query_evidence_mode="unknown")


# =============================================================================
# CLI, loss selection, and Lightning checkpoint reconstruction
# =============================================================================


class QueryPatchPipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from scripts.training_pipeline import execute_pipeline
        cls.pipeline = execute_pipeline

    def setUp(self):
        torch.manual_seed(31)
        self.options = dict(
            embedding_dim=8, projection_dim=8, contact_hidden_dim=4, contact_width=3,
            projection_dropout=0.0, learning_rate=0.001, interaction_loss_lambda=0.4,
            contact_threshold=0.5, loss_type="contact", interaction_module_type="max",
            trainable_final_slope=False, gamma_init=0.0, max_pooling=True,
            num_interaction_queries=3, query_heads=2, query_dropout=0.0,
            query_contact_bias_init=-1.0, compatibility_rank=4,
        )
        self.batch = dict(
            embedding_1=torch.randn(1, 3, 8), embedding_2=torch.randn(1, 4, 8),
            mask_1=torch.ones(1, 3, dtype=torch.bool), mask_2=torch.ones(1, 4, dtype=torch.bool),
            interaction_pair_mask=torch.ones(1, 3, 4, dtype=torch.bool),
            label=torch.ones(1), loss_weight=torch.ones(1),
            # -1 denotes unknown structure, not a negative contact target.
            contact_map=torch.tensor([[[1., 0., -1., 0.], [0., 0., 0., 0.], [0., -1., 0., 1.]]]),
        )

    def test_cli_defaults_to_gated_and_accepts_additive(self):
        for arguments, expected in (([], "gated"), (["--query-evidence-mode", "additive"], "additive")):
            with self.subTest(mode=expected), patch.object(sys, "argv", ["execute_pipeline", "--model", "query_patch", "--loss-type", "contact", *arguments]):
                args = self.pipeline.parse_args()
                self.assertEqual(args.query_evidence_mode, expected)

    def test_additive_shared_step_uses_only_original_task_losses(self):
        for head in ("bilinear", "dscript"):
            with self.subTest(head=head):
                module = self.pipeline.QueryPatchLightningModule(
                    **self.options, compatibility_head=head, query_evidence_mode="additive",
                    # A large auxiliary weight must still have no effect in
                    # historical mode; both extra objectives are disabled.
                    query_gate_loss_weight=100.0,
                ).eval()
                interaction, _, contact, auxiliary = module.model(
                    self.batch["embedding_1"], self.batch["embedding_2"],
                    self.batch["mask_1"], self.batch["mask_2"],
                    return_contact_map=True, return_contact_logits=True, return_contact_aux=True,
                )
                self.assertEqual(auxiliary, {})
                auxiliary_loss, metrics = module.compute_auxiliary_contact_loss(auxiliary, self.batch)
                self.assertEqual(auxiliary_loss.item(), 0.0)
                self.assertEqual(metrics, {})
                known = self.batch["contact_map"] >= 0
                expected = 0.4 * F.binary_cross_entropy_with_logits(interaction, self.batch["label"])
                expected = expected + 0.6 * F.binary_cross_entropy_with_logits(contact[known], self.batch["contact_map"][known])
                # Only Lightning logging needs a trainer; suppress it while
                # running the unchanged production loss orchestration itself.
                with patch.object(module, "log"):
                    actual = module.shared_step(self.batch, "train")
                torch.testing.assert_close(actual, expected)
                actual.backward()
                self.assertTrue(torch.isfinite(module.model.contact_module.queries.grad).all().item())

    def test_checkpoint_roundtrip_and_missing_mode_default(self):
        for mode in ("gated", "additive"):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as directory:
                module = self.pipeline.QueryPatchLightningModule(**self.options, query_evidence_mode=mode).eval()
                checkpoint = dict(
                    state_dict=module.state_dict(),
                    hyper_parameters=dict(module.hparams),
                    **{"pytorch-lightning_version": self.pipeline.pl.__version__},
                )
                self.assertEqual(checkpoint["hyper_parameters"]["query_evidence_mode"], mode)
                path = Path(directory) / "query_patch.ckpt"
                torch.save(checkpoint, path)
                restored = self.pipeline.QueryPatchLightningModule.load_from_checkpoint(path, map_location="cpu").eval()
                self.assertEqual(restored.model.contact_module.query_evidence_mode, mode)
                inputs = [self.batch[key] for key in ("embedding_1", "embedding_2", "mask_1", "mask_2")]
                torch.testing.assert_close(module(*inputs), restored(*inputs), rtol=0, atol=0)

                # Existing gated checkpoints predate the new hyperparameter.
                # Their parameter names and implicit default must still load.
                del checkpoint["hyper_parameters"]["query_evidence_mode"]
                torch.save(checkpoint, path)
                # Historical additive metadata also predates this flag and
                # requires an explicit mode override when reconstructing it.
                overrides = {} if mode == "gated" else {"query_evidence_mode": "additive"}
                restored = self.pipeline.QueryPatchLightningModule.load_from_checkpoint(
                    path, map_location="cpu", **overrides
                ).eval()
                self.assertEqual(restored.model.contact_module.query_evidence_mode, mode)
                torch.testing.assert_close(module(*inputs), restored(*inputs), rtol=0, atol=0)


if __name__ == "__main__":
    unittest.main()
