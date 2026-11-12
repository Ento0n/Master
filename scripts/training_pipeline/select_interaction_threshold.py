#!/usr/bin/env python
"""Select a decision threshold for an already trained interaction checkpoint.

Load the original model/data implementation, score validation pairs once in
evaluation mode, and maximize validation accuracy over all distinct decisions.
Freeze that scalar threshold before reading saved test predictions. Model
weights and scores stay fixed; only probability-to-label conversion changes.

The output directory contains validation predictions, the complete threshold
curve, a validation-only threshold record, and test results at both 0.5 and the
selected threshold, overall and by dimer type. Use --source-root for runs with
an archived implementation, or --report-existing to summarize saved results.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import json
from pathlib import Path
import sys

import numpy as np


# =============================================================================
# Threshold selection: no model inference or test data are needed here
# =============================================================================

def binary_metrics(labels, probabilities, threshold: float) -> dict:
    """Score binary labels (0/1) using positive := probability >= threshold."""
    labels = np.asarray(labels)
    predictions = np.asarray(probabilities) >= threshold
    tp = int(np.sum(predictions & (labels == 1)))
    fp = int(np.sum(predictions & (labels == 0)))
    tn = int(np.sum(~predictions & (labels == 0)))
    fn = int(np.sum(~predictions & (labels == 1)))
    return {
        "n": int(labels.size),
        "accuracy": (tp + tn) / labels.size,
        "precision": tp / (tp + fp) if tp + fp else 0.0,
        "recall": tp / (tp + fn) if tp + fn else 0.0,
        "f1": 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else 0.0,
        "tp": tp, "tn": tn, "fp": fp, "fn": fn,
        "predicted_positive": tp + fp,
    }


def validate_predictions(labels, probabilities) -> tuple[np.ndarray, np.ndarray]:
    """Reject empty, malformed, one-class, or non-probability prediction data."""
    labels = np.asarray(labels)
    probabilities = np.asarray(probabilities, dtype=np.float64)
    if labels.ndim != 1 or probabilities.ndim != 1 or labels.shape != probabilities.shape:
        raise ValueError("Labels and probabilities must be matching one-dimensional arrays")
    if not labels.size or not np.all(np.isin(labels, [0, 1])):
        raise ValueError("Labels must be a nonempty array containing only 0 and 1")
    if np.unique(labels).size != 2:
        raise ValueError("Both label classes must be present")
    if not np.all(np.isfinite(probabilities)) or np.any((probabilities < 0) | (probabilities > 1)):
        raise ValueError("Probabilities must be finite values in [0, 1]")
    return labels.astype(np.int64), probabilities


def select_threshold(labels, probabilities) -> tuple[float, list[dict]]:
    """Maximize validation accuracy, resolving ties without consulting test data.

    Equal probability values must always receive the same predicted label. A
    midpoint between each pair of distinct scores covers every possible split;
    0 covers all-positive and nextafter(max, +inf) covers all-negative. The
    latter can be slightly above 1 when a probability is exactly 1. Include 0.5
    explicitly so an equally good default decision is retained. Accuracy ties
    prefer the candidate closest to 0.5, then the smaller candidate.
    """
    labels, probabilities = validate_predictions(labels, probabilities)
    scores = np.unique(probabilities)
    midpoints = scores[:-1] + (scores[1:] - scores[:-1]) / 2
    # Adjacent floating-point scores may have no representable midpoint. Using
    # the upper score still implements the intended split under the >= rule.
    midpoints = np.where(midpoints > scores[:-1], midpoints, scores[1:])
    candidates = np.unique(np.r_[0.0, 0.5, midpoints, np.nextafter(scores[-1], np.inf)])

    # Sort once, then use cumulative class counts to evaluate all candidates in
    # O(n log n). searchsorted(left) keeps values equal to the cutoff positive.
    order = np.argsort(probabilities, kind="stable")
    sorted_scores = probabilities[order]
    cumulative_positives = np.r_[0, np.cumsum(labels[order])]
    positions = np.searchsorted(sorted_scores, candidates, side="left")
    total_positive = int(labels.sum())
    total_negative = labels.size - total_positive
    curve = []
    for threshold, position in zip(candidates, positions):
        fn = int(cumulative_positives[position])
        tn = int(position) - fn
        tp, fp = total_positive - fn, total_negative - tn
        curve.append({"threshold": float(threshold), "accuracy": (tp + tn) / labels.size,
                      "tp": tp, "tn": tn, "fp": fp, "fn": fn})
    best = min(curve, key=lambda row: (-(row["tp"] + row["tn"]),
                                      abs(row["threshold"] - 0.5), row["threshold"]))
    return best["threshold"], curve


def average_precision(labels: np.ndarray, probabilities: np.ndarray) -> float:
    """Exact ranking AP with tied scores grouped, independent of the cutoff."""
    order = np.argsort(-probabilities, kind="stable")
    ends = np.r_[np.flatnonzero(np.diff(probabilities[order])), labels.size - 1]
    true_positives = np.cumsum(labels[order])[ends]
    recall = true_positives / labels.sum()
    return float(np.sum(np.diff(np.r_[0.0, recall]) * true_positives / (ends + 1)))


def sha256(path: Path) -> str:
    """Hash provenance files in bounded chunks, including model checkpoints."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


# =============================================================================
# Accuracy reporting: reuse the frozen threshold for every dimer-type subgroup
# =============================================================================

def write_accuracy_report(output_dir: Path) -> None:
    """Report overall and subgroup accuracy from completed prediction exports.

    Each subgroup includes every row of that dimer type, including labels 0
    (non-interacting) and 1 (interacting). No subgroup threshold is fitted.
    Existing global metrics must match before any report is updated. Prediction
    tables and threshold.json remain the original inference/selection records.
    """
    summary = json.loads((output_dir / "summary.json").read_text())
    threshold = summary["threshold"]
    compact = {"model": summary["model"],
               "interaction_head": summary["model_hyperparameters"]["interaction_module_type"],
               "threshold": threshold}
    rows, prediction_hashes = [], {}
    for split in ("validation", "test"):
        path = output_dir / f"{split}_predictions.tsv"
        with path.open() as handle:
            predictions = list(csv.DictReader(handle, delimiter="\t"))
        labels, probabilities = validate_predictions(
            [int(row["label"]) for row in predictions],
            [float(row["predicted_probability"]) for row in predictions])
        expected_split = "val" if split == "validation" else "test"
        if any(row["split"] != expected_split for row in predictions):
            raise ValueError(f"Wrong split in {path}")
        # Unknown/other dimer types remain visible and count toward the overall
        # accuracy. Do not silently discard rows outside the usual homo/hetero.
        dimer_types = np.asarray([row["dimer_type"].strip().lower() or "unknown" for row in predictions])
        groups = ["overall"] + sorted(set(dimer_types), key=lambda group: (group != "homo", group != "hetero", group))
        if "overall" in dimer_types:
            raise ValueError("'overall' is reserved for the complete prediction cohort")
        by_type = {}
        for group in groups:
            mask = np.ones(labels.size, dtype=bool) if group == "overall" else dimer_types == group
            metrics = {"at_0.5": binary_metrics(labels[mask], probabilities[mask], 0.5),
                       "at_selected_threshold": binary_metrics(labels[mask], probabilities[mask], threshold)}
            # This also makes report-only refreshes fail if their saved scores
            # no longer reproduce the completed evaluation's global metrics.
            if group == "overall":
                for decision, values in metrics.items():
                    if values != summary[split][decision]:
                        raise ValueError(f"Saved {split} predictions do not reproduce {decision} metrics")
            else:
                by_type[group] = metrics
            row = {"split": split, "dimer_type": group, "n": int(mask.sum()),
                   "positive": int(labels[mask].sum()), "negative": int((labels[mask] == 0).sum()),
                   "selected_threshold": threshold}
            prefix = split if group == "overall" else f"{split}_{group}"
            compact[f"{prefix}_n"] = row["n"]
            for decision, suffix in (("at_0.5", "0.5"), ("at_selected_threshold", "selected")):
                values = metrics[decision]
                row[f"correct_{suffix}"] = values["tp"] + values["tn"]
                row[f"accuracy_{suffix}"] = values["accuracy"]
                compact[f"{prefix}_accuracy_{suffix}"] = values["accuracy"]
            rows.append(row)
        summary[split]["by_dimer_type"] = by_type
        prediction_hashes[path.name] = sha256(path)
    compact["test_interaction_auprc"] = summary["test"]["interaction_auprc"]
    summary["accuracy_reporting"] = {"script_sha256": sha256(Path(__file__)),
                                     "prediction_sha256": prediction_hashes,
                                     "threshold_scope": "one frozen validation threshold per checkpoint"}
    for filename, records in (("accuracy_by_dimer_type.tsv", rows), ("summary.tsv", [compact])):
        with (output_dir / filename).open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(records)
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(compact, indent=2), flush=True)


# =============================================================================
# Restore a frozen checkpoint and its original validation-data preparation
# =============================================================================

def validate_prediction_cohort(predictions, examples, split_name: str) -> None:
    """Require the original ordered pair identities, treating missing IDs equally.

    Pandas versions differ in whether astype(str) preserves missing values or
    spells them as text. Compare numeric assembly IDs with equal_nan=True so
    missing synthetic-negative IDs compare equal without hiding real changes.
    """
    import pandas as pd

    if "split" not in predictions or not predictions["split"].eq(split_name).all():
        raise ValueError(f"Saved predictions must contain only split={split_name} rows")
    examples = examples.reset_index(drop=True)
    if len(predictions) != len(examples):
        raise ValueError("Saved prediction count differs from the original dataset")
    for column in ("new_cluster_pair", "dimer_type", "label", "pdb_id", "assembly_number"):
        actual, expected = predictions[column], examples[column]
        if column in ("label", "assembly_number"):
            actual = pd.to_numeric(actual.replace("", np.nan)).to_numpy(dtype=np.float64, na_value=np.nan)
            expected = pd.to_numeric(expected).to_numpy(dtype=np.float64, na_value=np.nan)
            matches = np.array_equal(actual, expected, equal_nan=True)
        else:
            actual = actual.fillna("").astype(str)
            expected = expected.fillna("").astype(str)
            if column == "pdb_id":
                actual, expected = actual.str.lower(), expected.str.lower()
            matches = np.array_equal(actual.to_numpy(), expected.to_numpy())
        if not matches:
            raise ValueError(f"Saved {split_name} cohort/order does not match: {column}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-existing", type=Path,
                        help="Refresh overall/subgroup reports in a completed output directory; no inference")
    parser.add_argument("--model", choices=("dscript", "query_patch"))
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--embedding-manifest", type=Path)
    parser.add_argument("--interactions", type=Path,
                        default=Path("/nfs/scratch/pdb_dimers/final_datasets/balanced_interactions_strict_cd_hit.tsv"))
    parser.add_argument("--contact-map-dir", type=Path, default=Path("/nfs/scratch/pdb_dimers/contact_maps/data"))
    parser.add_argument("--source-root", type=Path, default=Path(__file__).resolve().parents[2],
                        help="Repository or archived source containing scripts/training_pipeline/execute_pipeline.py")
    parser.add_argument("--submission-record", type=Path,
                        help="Optional submission.json; verify its source and input hashes before inference")
    parser.add_argument("--test-predictions", type=Path,
                        help="Defaults to test_predictions.tsv beside the checkpoint directory")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--reuse-validation-dir", type=Path,
                        help="Reuse a previous validation export after verifying checkpoint/data provenance")
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cpu")
    parser.add_argument("--batch-size", type=int, default=8,
                        help="Match the original evaluation batch size (ProstT5 final runs: 8)")
    parser.add_argument("--num-workers", type=int, default=4)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--verify-test-batches", type=int, default=1,
                        help="Check this many original test batches against saved scores; 0 disables the check")
    parser.add_argument("--verification-atol", type=float, default=1e-4,
                        help="Absolute tolerance for CPU/GPU score reproduction (default: 0.0001)")
    args = parser.parse_args()
    if args.report_existing:
        if any(value is not None for value in (args.model, args.checkpoint, args.embedding_manifest, args.output_dir)):
            parser.error("Use --report-existing alone; it reads the completed run's saved configuration")
        write_accuracy_report(args.report_existing.resolve(strict=True))
        return
    for name in ("model", "checkpoint", "embedding_manifest", "output_dir"):
        if getattr(args, name) is None:
            parser.error(f"--{name.replace('_', '-')} is required unless --report-existing is used")
    if args.batch_size < 1 or args.num_workers < 0 or args.verify_test_batches < 0:
        parser.error("batch-size must be positive; num-workers/verify-test-batches must be nonnegative")
    if not np.isfinite(args.verification_atol) or args.verification_atol < 0:
        parser.error("verification-atol must be finite and nonnegative")
    args.checkpoint = args.checkpoint.resolve(strict=True)
    args.source_root = args.source_root.resolve(strict=True)
    args.output_dir = args.output_dir.resolve()
    test_path = (args.test_predictions or args.checkpoint.parent.parent / "test_predictions.tsv").resolve(strict=True)
    if args.output_dir.exists():
        parser.error("output-dir already exists; choose a new directory to preserve previous results")

    # Verify archived code and the tabular input identities before loading a
    # checkpoint. No training modules are imported before selecting source-root.
    source_files = ["scripts/models/__init__.py", "scripts/models/dscript.py",
                    "scripts/models/query_patch.py", "scripts/models/fully_connected.py",
                    "scripts/training_pipeline/execute_pipeline.py",
                    "scripts/training_pipeline/losses.py", "scripts/training_pipeline/metrics.py"]
    source_hashes = {name: sha256(args.source_root / name) for name in source_files}
    input_hashes = {str(path.resolve()): sha256(path) for path in (args.interactions, args.embedding_manifest)}
    if args.submission_record:
        record = json.loads(args.submission_record.read_text())
        for name, digest in source_hashes.items():
            if record["source_sha256"].get(name) != digest:
                raise ValueError(f"Archived source hash mismatch: {name}")
        for name, digest in input_hashes.items():
            if record["input_sha256"].get(name) != digest:
                raise ValueError(f"Original input hash mismatch: {name}")

    sys.path.insert(0, str(args.source_root))
    pipeline = importlib.import_module("scripts.training_pipeline.execute_pipeline")
    if Path(pipeline.__file__).resolve() != args.source_root / source_files[4]:
        raise RuntimeError("The imported training pipeline does not match --source-root")
    import pandas as pd
    import torch

    pipeline.pl.seed_everything(args.seed, workers=True)
    torch.set_float32_matmul_precision("high")
    if args.device == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but is unavailable; run in a GPU allocation")
    device = torch.device(args.device)
    model_class = pipeline.DScriptLightningModule if args.model == "dscript" else pipeline.QueryPatchLightningModule
    model = model_class.load_from_checkpoint(str(args.checkpoint), map_location="cpu")
    model.eval().to(device)
    model.requires_grad_(False)
    data = pipeline.DScriptDataModule(
        interaction_path=args.interactions, embedding_manifest=args.embedding_manifest,
        contact_map_dir=args.contact_map_dir, batch_size=args.batch_size,
        num_workers=args.num_workers, loss_mode="basic", limit_rows=None,
        negative_contact_maps=False,
    )
    # The original collator retains unknown-residue and padding masks. Contact
    # maps are needed for these interaction masks even though no loss is run.
    data.setup("fit")
    if data.val_dataset is None or not len(data.val_dataset):
        raise ValueError("No validation examples are available")
    provenance = {
        "model": args.model, "checkpoint": str(args.checkpoint), "checkpoint_sha256": sha256(args.checkpoint),
        "source_root": str(args.source_root), "source_sha256": source_hashes, "input_sha256": input_hashes,
        "model_hyperparameters": dict(model.hparams), "batch_size": args.batch_size, "seed": args.seed,
        "contact_map_dir": str(args.contact_map_dir.resolve()), "shuffle": False,
        "float32_matmul_precision": torch.get_float32_matmul_precision(),
        "threshold_script_sha256": sha256(Path(__file__)),
        "submission_record": str(args.submission_record.resolve()) if args.submission_record else None,
        "submission_record_sha256": sha256(args.submission_record) if args.submission_record else None,
        "device": str(device), "torch_version": torch.__version__, "numpy_version": np.__version__,
    }
    if args.reuse_validation_dir:
        previous = json.loads((args.reuse_validation_dir / "threshold.json").read_text())
        # Recovery after an export/check failure must reuse the same checkpoint,
        # masks, inputs, source implementation, and evaluation configuration.
        for key in ("model", "checkpoint_sha256", "source_sha256", "input_sha256",
                    "model_hyperparameters", "batch_size", "seed", "contact_map_dir",
                    "shuffle", "float32_matmul_precision", "device", "torch_version"):
            if previous[key] != provenance[key]:
                raise ValueError(f"Cannot reuse validation predictions with different {key}")
        validation_path = args.reuse_validation_dir / "validation_predictions.tsv"
        validation = pd.read_csv(validation_path, sep="\t", float_precision="round_trip")
        validate_prediction_cohort(validation, data.val_dataset.examples, "val")
        provenance["reused_validation_predictions"] = str(validation_path.resolve())
        provenance["reused_validation_predictions_sha256"] = sha256(validation_path)
        print(f"Reusing verified validation predictions from {validation_path}", flush=True)
    else:
        print(f"Scoring {len(data.val_dataset)} validation pairs with {args.checkpoint.name}", flush=True)
        validation = pipeline.collect_split_predictions("val", data.val_dataset, data, model, device)
    labels, probabilities = validate_predictions(validation["label"], validation["predicted_probability"])
    threshold, curve = select_threshold(labels, probabilities)
    validation_metrics = {
        "at_0.5": binary_metrics(labels, probabilities, 0.5),
        "at_selected_threshold": binary_metrics(labels, probabilities, threshold),
        "interaction_auprc": average_precision(labels, probabilities),
    }
    if args.reuse_validation_dir and (threshold != previous["threshold"] or validation_metrics != previous["validation"]):
        raise ValueError("Reused validation predictions do not reproduce the recorded threshold/metrics")
    threshold_record = {
        **provenance, "threshold": threshold, "selection_split": "val", "selection_metric": "accuracy",
        "tie_break": "closest to 0.5, then smaller threshold", "validation": validation_metrics,
    }
    args.output_dir.mkdir(parents=True)
    # Commit the decision record before reading any saved test probabilities or
    # labels. Subsequent test metrics cannot influence threshold selection.
    (args.output_dir / "threshold.json").write_text(json.dumps(threshold_record, indent=2) + "\n")
    validation["predicted_label_selected_threshold"] = (probabilities >= threshold).astype(np.int64)
    validation.to_csv(args.output_dir / "validation_predictions.tsv", sep="\t", index=False)
    with (args.output_dir / "validation_threshold_curve.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(curve[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(curve)
    print(f"Frozen validation threshold: {threshold:.17g}; accuracy "
          f"{validation_metrics['at_0.5']['accuracy']:.4%} -> "
          f"{validation_metrics['at_selected_threshold']['accuracy']:.4%}", flush=True)

    # Saved test outputs must describe the original test cohort in the same
    # ordered batches. Read identifiers as strings to preserve missing values
    # and identifiers that resemble numbers (e.g. a PDB ID).
    test = pd.read_csv(test_path, sep="\t", dtype=str, keep_default_na=False)
    examples = data.test_dataset.examples.reset_index(drop=True)
    validate_prediction_cohort(test, examples, "test")
    test_labels, test_probabilities = validate_predictions(
        pd.to_numeric(test["label"]).to_numpy(), pd.to_numeric(test["predicted_probability"]).to_numpy())
    if not np.array_equal(pd.to_numeric(test["predicted_label"]).to_numpy(), test_probabilities >= 0.5):
        raise ValueError("Saved test labels do not match the original >=0.5 rule")

    verification = {"n": 0}
    if args.verify_test_batches:
        count = min(args.verify_test_batches * args.batch_size, len(test))
        subset = pipeline.DScriptDataset(examples.iloc[:count].copy())
        check = pipeline.collect_split_predictions("test", subset, data, model, device)
        reproduced = check["predicted_probability"].to_numpy()
        differences = np.abs(reproduced - test_probabilities[:count])
        # CPU and GPU kernels need not be bitwise identical, particularly when
        # the original GPU run used matmul precision='high'. Preserve every
        # probed score and the tolerance, including on a failed verification.
        check["original_predicted_probability"] = test_probabilities[:count]
        check["absolute_probability_error"] = differences
        check.to_csv(args.output_dir / "test_score_verification.tsv", sep="\t", index=False)
        verification = {
            "n": count, "max_absolute_probability_error": float(differences.max()),
            "mean_absolute_probability_error": float(differences.mean()),
            "absolute_tolerance": args.verification_atol, "relative_tolerance": 2e-5,
            "classification_disagreements_at_0.5": int(np.sum((reproduced >= 0.5) != (test_probabilities[:count] >= 0.5))),
            "classification_disagreements_at_selected_threshold": int(np.sum((reproduced >= threshold) != (test_probabilities[:count] >= threshold))),
        }
        if not np.allclose(reproduced, test_probabilities[:count], atol=args.verification_atol, rtol=2e-5):
            raise ValueError(f"Checkpoint does not reproduce saved test probabilities: {verification}")

    test_metrics = {
        "at_0.5": binary_metrics(test_labels, test_probabilities, 0.5),
        "at_selected_threshold": binary_metrics(test_labels, test_probabilities, threshold),
        "interaction_auprc": average_precision(test_labels, test_probabilities),
    }
    summary = {**threshold_record, "test": test_metrics, "test_predictions": str(test_path),
               "test_predictions_sha256": sha256(test_path), "test_score_verification": verification}
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    test["predicted_label_selected_threshold"] = (test_probabilities >= threshold).astype(np.int64)
    test.to_csv(args.output_dir / "test_predictions.tsv", sep="\t", index=False)
    write_accuracy_report(args.output_dir)


if __name__ == "__main__":
    main()
