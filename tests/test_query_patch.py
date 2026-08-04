"""Focused invariants for the learned-query contact-map architecture."""

from __future__ import annotations

import unittest

import torch
import torch.nn.functional as F

from scripts.models import DScriptInteractionModel, QueryPatchInteractionModel


class QueryPatchInteractionModelTest(unittest.TestCase):
    """Check tensor contracts, masks, symmetry, and gradient connectivity."""

    @staticmethod
    def make_model(
        interaction_module_type: str = "normal",
        max_pooling: bool = False,
    ) -> QueryPatchInteractionModel:
        """Return a deterministic small model suitable for CPU tests."""
        return QueryPatchInteractionModel(
            embedding_dim=12,
            projection_dim=8,
            projection_dropout=0.0,
            use_weight_matrix=False,
            interaction_module_type=interaction_module_type,
            max_pooling=max_pooling,
            num_queries=3,
            query_heads=2,
            query_layers=1,
            query_dropout=0.0,
        )

    def setUp(self) -> None:
        """Create unequal protein lengths with explicit batch padding."""
        torch.manual_seed(7)
        self.embeddings_1 = torch.randn(2, 5, 12)
        self.embeddings_2 = torch.randn(2, 6, 12)
        self.mask_1 = torch.tensor(
            [[True, True, True, False, False], [True, True, True, True, True]]
        )
        self.mask_2 = torch.tensor(
            [[True, True, True, True, False, False], [True, True, True, True, True, True]]
        )
        self.pair_mask = self.mask_1[:, :, None] & self.mask_2[:, None, :]

    def test_public_output_contract_for_both_interaction_heads(self) -> None:
        """Both inherited PPI heads must consume the new map without API changes."""
        configurations = (("normal", False), ("normal", True), ("max", False))
        for interaction_module_type, max_pooling in configurations:
            with self.subTest(head=interaction_module_type, local_max_pool=max_pooling):
                model = self.make_model(interaction_module_type, max_pooling).eval()
                interaction_logits, contacts, contact_logits = model(
                    self.embeddings_1,
                    self.embeddings_2,
                    self.mask_1,
                    self.mask_2,
                    return_contact_map=True,
                    return_contact_logits=True,
                )

                self.assertEqual(interaction_logits.shape, (2,))
                self.assertEqual(contact_logits.shape, (2, 5, 6))
                self.assertEqual(contacts.shape, (2, 5, 6))
                self.assertTrue(torch.isfinite(interaction_logits).all())
                self.assertTrue(torch.isfinite(contact_logits).all())

                # Valid map values are sigmoid probabilities. Both raw logits
                # and probabilities are exactly zeroed outside real residues.
                torch.testing.assert_close(
                    contacts[self.pair_mask],
                    torch.sigmoid(contact_logits[self.pair_mask]),
                )
                self.assertTrue((contact_logits[~self.pair_mask] == 0).all())
                self.assertTrue((contacts[~self.pair_mask] == 0).all())

                if interaction_module_type == "max":
                    expected = contact_logits.masked_fill(~self.pair_mask, -torch.inf).amax(dim=(1, 2))
                    torch.testing.assert_close(interaction_logits, expected)

    def test_eval_swap_equivariance(self) -> None:
        """Swapping proteins must transpose contacts and preserve the PPI logit."""
        configurations = (("normal", False), ("normal", True), ("max", False))
        for interaction_module_type, max_pooling in configurations:
            with self.subTest(head=interaction_module_type, local_max_pool=max_pooling):
                model = self.make_model(interaction_module_type, max_pooling).eval()
                interaction_ab, _, logits_ab = model(
                    self.embeddings_1,
                    self.embeddings_2,
                    self.mask_1,
                    self.mask_2,
                    return_contact_map=True,
                    return_contact_logits=True,
                )
                interaction_ba, _, logits_ba = model(
                    self.embeddings_2,
                    self.embeddings_1,
                    self.mask_2,
                    self.mask_1,
                    return_contact_map=True,
                    return_contact_logits=True,
                )

                torch.testing.assert_close(logits_ab, logits_ba.transpose(1, 2), rtol=1e-5, atol=1e-6)
                torch.testing.assert_close(interaction_ab, interaction_ba, rtol=1e-5, atol=1e-6)

    def test_masked_embedding_values_cannot_affect_predictions(self) -> None:
        """Attention and pooling must ignore arbitrary values in padding rows."""
        model = self.make_model().eval()
        changed_1 = self.embeddings_1.clone()
        changed_2 = self.embeddings_2.clone()
        changed_1[~self.mask_1] = 100.0 * torch.randn_like(changed_1[~self.mask_1])
        changed_2[~self.mask_2] = 100.0 * torch.randn_like(changed_2[~self.mask_2])

        interaction_original, logits_original = model(
            self.embeddings_1,
            self.embeddings_2,
            self.mask_1,
            self.mask_2,
            return_contact_logits=True,
        )
        interaction_changed, logits_changed = model(
            changed_1,
            changed_2,
            self.mask_1,
            self.mask_2,
            return_contact_logits=True,
        )

        torch.testing.assert_close(
            logits_original[self.pair_mask],
            logits_changed[self.pair_mask],
            rtol=1e-5,
            atol=1e-6,
        )
        torch.testing.assert_close(interaction_original, interaction_changed, rtol=1e-5, atol=1e-6)

    def test_contact_and_interaction_losses_reach_all_new_components(self) -> None:
        """Joint supervision must update queries, projection, masks, and strengths."""
        model = self.make_model().train()
        interaction_logits, contact_logits = model(
            self.embeddings_1,
            self.embeddings_2,
            self.mask_1,
            self.mask_2,
            return_contact_logits=True,
        )

        contact_targets = torch.zeros_like(contact_logits)
        contact_targets[0, 0, 0] = 1.0
        contact_targets[1, 2, 3] = 1.0
        interaction_targets = torch.tensor([1.0, 0.0])
        loss = F.binary_cross_entropy_with_logits(
            contact_logits[self.pair_mask],
            contact_targets[self.pair_mask],
        ) + F.binary_cross_entropy_with_logits(interaction_logits, interaction_targets)
        loss.backward()

        parameters = {
            "queries": model.contact_module.queries,
            "shared projection": model.projection.network[0].weight,
            "query mask projection": model.contact_module.query_mask_projection.weight,
            "residue mask projection": model.contact_module.residue_mask_projection.weight,
            "patch strength": model.contact_module.patch_strength.weight,
        }
        for name, parameter in parameters.items():
            with self.subTest(parameter=name):
                self.assertIsNotNone(parameter.grad)
                self.assertTrue(torch.isfinite(parameter.grad).all())
                self.assertGreater(parameter.grad.abs().sum().item(), 0.0)

    def test_patch_diagnostics_and_default_masks(self) -> None:
        """Patch details keep their documented shapes and support unpadded calls."""
        model = self.make_model().eval()
        embeddings_1 = self.embeddings_1[:1, :3]
        embeddings_2 = self.embeddings_2[:1, :4]
        all_true_1 = torch.ones((1, 3), dtype=torch.bool)
        all_true_2 = torch.ones((1, 4), dtype=torch.bool)

        no_mask_logits = model.contact_logits(embeddings_1, embeddings_2)
        explicit_logits, masks_1, masks_2, strengths = model.query_patch_outputs(
            embeddings_1,
            embeddings_2,
            all_true_1,
            all_true_2,
        )

        torch.testing.assert_close(no_mask_logits, explicit_logits)
        self.assertEqual(masks_1.shape, (1, 3, 3))
        self.assertEqual(masks_2.shape, (1, 3, 4))
        self.assertEqual(strengths.shape, (1, 3))
        self.assertTrue(((masks_1 >= 0) & (masks_1 <= 1)).all())
        self.assertTrue(((masks_2 >= 0) & (masks_2 <= 1)).all())
        self.assertTrue((strengths > 0).all())

    def test_invalid_head_count_is_rejected(self) -> None:
        """Attention width must split evenly between transformer heads."""
        with self.assertRaises(ValueError):
            QueryPatchInteractionModel(
                embedding_dim=12,
                projection_dim=10,
                num_queries=3,
                query_heads=4,
            )

    def test_single_residue_proteins_are_finite(self) -> None:
        """The decoder and both interaction heads support the smallest maps."""
        embedding_1 = torch.randn(2, 1, 12)
        embedding_2 = torch.randn(2, 1, 12)
        for interaction_module_type in ("normal", "max"):
            with self.subTest(head=interaction_module_type):
                model = self.make_model(interaction_module_type, max_pooling=True).eval()
                interaction_logits, contact_logits = model(
                    embedding_1,
                    embedding_2,
                    return_contact_logits=True,
                )
                self.assertEqual(contact_logits.shape, (2, 1, 1))
                self.assertTrue(torch.isfinite(contact_logits).all())
                self.assertTrue(torch.isfinite(interaction_logits).all())

    def test_inherited_structural_mask_controls_query_patch_pooling(self) -> None:
        """The new map generator must retain the established interaction-mask API."""
        model = self.make_model(interaction_module_type="max").eval()
        interaction_pair_mask = torch.zeros_like(self.pair_mask)
        interaction_pair_mask[0, 1, 2] = True
        interaction_pair_mask[1, 3, 4] = True

        interaction_logits, contact_logits = model(
            self.embeddings_1,
            self.embeddings_2,
            self.mask_1,
            self.mask_2,
            interaction_pair_mask=interaction_pair_mask,
            return_contact_logits=True,
        )

        pooling_mask = self.pair_mask & interaction_pair_mask
        expected = contact_logits.masked_fill(~pooling_mask, -torch.inf).amax(dim=(1, 2))
        torch.testing.assert_close(interaction_logits, expected)

    def test_original_dscript_interaction_mask_remains_compatible(self) -> None:
        """Adding query-patch must not change D-SCRIPT's structural pooling API."""
        model = DScriptInteractionModel(
            embedding_dim=12,
            projection_dim=8,
            contact_hidden_dim=4,
            contact_width=3,
            projection_dropout=0.0,
            interaction_module_type="max",
        ).eval()
        interaction_pair_mask = torch.zeros_like(self.pair_mask)
        interaction_pair_mask[0, 0, 1] = True
        interaction_pair_mask[1, 2, 3] = True

        interaction_logits, contact_logits = model(
            self.embeddings_1,
            self.embeddings_2,
            self.mask_1,
            self.mask_2,
            interaction_pair_mask=interaction_pair_mask,
            return_contact_logits=True,
        )
        self.assertEqual(interaction_logits.shape, (2,))
        self.assertEqual(contact_logits.shape, (2, 5, 6))

        pooling_mask = self.pair_mask & interaction_pair_mask
        expected = contact_logits.masked_fill(~pooling_mask, -torch.inf).amax(dim=(1, 2))
        torch.testing.assert_close(interaction_logits, expected)


if __name__ == "__main__":
    unittest.main()
