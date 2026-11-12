#!/usr/bin/env python3
"""Select one contact cutoff per frozen checkpoint using validation macro F1.

Restore the original model and padded batches, evaluate every known cell of
real structural maps, and search all distinct validation decisions. Save the
cutoff before test inference. Report per-map and pooled metrics at 0.5 and at
the frozen cutoff, alongside threshold-independent contact AP. No training or
interaction-threshold changes occur. Raw known-cell scores are cached so the
contact analysis is reproducible without another model forward pass.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import importlib
import json
from pathlib import Path
import sys

import numpy as np

# Import adjacent analysis helpers without importing the repository's scripts
# package: the archived training package must be selected by --source-root.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from contact_threshold_metrics import (ExactContactThresholdSearch,
                                       aggregate_contact_metrics, contact_average_precision,
                                       contact_metrics)
from select_interaction_threshold import sha256, validate_prediction_cohort


# =============================================================================
# Frozen-checkpoint inference and reproducible known-cell exports
# =============================================================================

def collect_contacts(pipeline, data, model, dataset, device, output_dir, split,
                     threshold=None, sample_count=5):
    """Return real-map known labels/scores plus identities in original order.

    The forward pass uses ALL pairs in cohort order at the original batch size,
    including negatives, because padding can change convolution outputs. Only real structural maps
    with known cells enter metrics. Labels -1 (unknown/padding) are excluded;
    known labels are 0/1 and contact scores are probabilities in [0,1].
    """
    import torch

    maps, rows, interaction_scores = [], [], []
    candidates = np.flatnonzero(dataset.examples["contact_map_path"].notna().to_numpy())
    sampled = set(np.random.default_rng(42).choice(
        candidates, size=min(sample_count, len(candidates)), replace=False).tolist())
    skipped_unknown = skipped_nonstructural = offset = 0
    loader = pipeline.make_prediction_dataloader(data, dataset)
    with torch.inference_mode():
        for batch_number, batch in enumerate(loader):
            batch = pipeline.move_tensor_batch_to_device(batch, device)
            logits, probabilities = model.model(
                batch["embedding_1"], batch["embedding_2"], batch["mask_1"], batch["mask_2"],
                interaction_pair_mask=batch["interaction_pair_mask"], return_contact_map=True)
            interaction_scores.extend(torch.sigmoid(logits).cpu().tolist())
            # These arrays have shape [B,L1,L2]; keep float32 scores exactly as
            # emitted, but threshold comparisons use float64 in the metric code.
            scores = probabilities.cpu().numpy()
            targets = batch["contact_map"].cpu().numpy()
            structural = batch["has_true_contact_map"].cpu().numpy()
            for index in range(len(scores)):
                if not structural[index]:
                    skipped_nonstructural += 1
                    continue
                target = targets[index]
                known = target >= 0
                if not known.any():
                    skipped_unknown += 1
                    continue
                labels = target[known].copy()
                score = scores[index][known].copy()
                # Validate before converting labels to uint8, which could hide
                # malformed fractional values. Also rejects NaNs/infinite scores.
                contact_metrics(labels, score, 0.5)
                labels = labels.astype(np.uint8)
                maps.append((labels, score))
                split_index = offset + index
                example = dataset.examples.iloc[split_index]
                map_path = Path(example["contact_map_path"])
                rows.append({"split_index": split_index, "split": split,
                             "pdb_id": str(example["pdb_id"]),
                             "assembly_number": int(example["assembly_number"]),
                             "new_cluster_pair": str(example["new_cluster_pair"]),
                             "dimer_type": str(example["dimer_type"]),
                             "contact_map_path": str(map_path), "contact_map_sha256": sha256(map_path),
                             "known_cells": int(labels.size), "positive_cells": int(labels.sum())})
                if threshold is not None and split_index in sampled:
                    # Export the same deterministic five complexes for all runs,
                    # using the evaluated batch outputs, not singleton re-inference.
                    height = int(batch["mask_1"][index].sum().item())
                    width = int(batch["mask_2"][index].sum().item())
                    truth = target[:height, :width].astype(np.int8)
                    full_scores = scores[index, :height, :width].astype(np.float64)
                    stem = f"test_{split_index:04d}_{example['pdb_id']}_assembly{int(example['assembly_number'])}"
                    sample_dir = output_dir / "contact_maps"
                    for name, cutoff in (("0.5", 0.5), ("selected", threshold)):
                        binary = np.where(truth >= 0, full_scores >= cutoff, -1).astype(np.int8)
                        pipeline.save_contact_map_png(binary, sample_dir / f"{stem}_{name}.png",
                                                      f"Predicted contacts (threshold >= {cutoff:.6g})")
                    pipeline.save_contact_map_png(truth, sample_dir / f"{stem}_true.png", "True contacts")
                    np.savez_compressed(sample_dir / f"{stem}.npz", probabilities=full_scores,
                                        labels=truth, selected_threshold=threshold)
            offset += len(scores)
            if (batch_number + 1) % 25 == 0:
                print(f"{split}: {offset}/{len(dataset)} pairs; {len(maps)} eligible maps", flush=True)
    if offset != len(dataset) or not maps:
        raise ValueError(f"Incomplete or empty {split} contact cohort")
    # Offsets delimit maps in flattened arrays; TSV split_index maps them back
    # to the original cohort. Each map preserves row-major known-cell order.
    np.savez_compressed(output_dir / f"{split}_contacts.npz",
                        labels=np.concatenate([item[0] for item in maps]),
                        probabilities=np.concatenate([item[1] for item in maps]),
                        offsets=np.r_[0, np.cumsum([len(item[0]) for item in maps])],
                        split_indices=np.array([row["split_index"] for row in rows]))
    exclusions = {"total_pairs": offset, "without_real_map": skipped_nonstructural,
                  "real_maps_without_known_cells": skipped_unknown}
    return maps, rows, np.asarray(interaction_scores), exclusions


def summarize_contacts(maps, rows, threshold, output_path):
    """Report exact macro AP and 4096-bin micro AP, plus both binary cutoffs.

    Macro metrics average map-level ratios equally. Micro metrics use a single
    split-wide confusion matrix, unlike historical batch-averaged diagnostics.
    Undefined precision/recall/F1 and AP for all-negative maps are zero.
    """
    decisions = {"at_0.5": [], "at_selected_threshold": []}
    positive_bins = np.zeros(4096, dtype=np.int64)
    negative_bins = np.zeros(4096, dtype=np.int64)
    report_rows = []
    for (labels, scores), metadata in zip(maps, rows):
        row = dict(metadata)
        row["auprc"] = contact_average_precision(labels, scores)
        for name, cutoff in (("at_0.5", 0.5), ("at_selected_threshold", threshold)):
            metrics = contact_metrics(labels, scores, cutoff)
            decisions[name].append(metrics)
            row.update({f"{name}_{key}": value for key, value in metrics.items()})
        report_rows.append(row)
        # Match BinnedBinaryAveragePrecision's float32 bin assignment, including
        # clamping probability 1 into the final bin. Unknown cells were removed.
        bins = np.minimum((scores * 4096).astype(np.int64), 4095)
        positive_bins += np.bincount(bins[labels == 1], minlength=4096)
        negative_bins += np.bincount(bins[labels == 0], minlength=4096)
    tp = np.cumsum(positive_bins[::-1], dtype=np.float64)
    predicted = tp + np.cumsum(negative_bins[::-1])
    precision = np.divide(tp, predicted, out=np.zeros_like(tp), where=predicted > 0)
    micro_ap = float(np.sum(precision * positive_bins[::-1]) / max(1, positive_bins.sum()))
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(report_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(report_rows)
    return {**{name: aggregate_contact_metrics(values) for name, values in decisions.items()},
            "contact_auprc_macro": float(np.mean([row["auprc"] for row in report_rows])),
            "contact_auprc_micro": micro_ap}


# =============================================================================
# CLI orchestration: validation decision is persisted before test evaluation
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", choices=("dscript", "query_patch"), required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--embedding-manifest", type=Path, required=True)
    parser.add_argument("--interactions", type=Path,
                        default=Path("/nfs/scratch/pdb_dimers/final_datasets/balanced_interactions_strict_cd_hit.tsv"))
    parser.add_argument("--contact-map-dir", type=Path, default=Path("/nfs/scratch/pdb_dimers/contact_maps/data"))
    parser.add_argument("--source-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--submission-record", type=Path)
    parser.add_argument("--test-predictions", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cpu")
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--num-workers", type=int, default=4)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--selection-metric", choices=("macro_f1", "micro_f1", "micro_accuracy"), default="macro_f1")
    parser.add_argument("--sample-count", type=int, default=5)
    parser.add_argument("--verification-atol", type=float, default=1e-4)
    args = parser.parse_args()
    if args.batch_size < 1 or args.num_workers < 0 or args.sample_count < 0:
        parser.error("batch-size must be positive; num-workers and sample-count must be nonnegative")
    if not np.isfinite(args.verification_atol) or args.verification_atol < 0:
        parser.error("verification-atol must be finite and nonnegative")
    args.checkpoint = args.checkpoint.resolve(strict=True)
    args.source_root = args.source_root.resolve(strict=True)
    args.output_dir = args.output_dir.resolve()
    if args.output_dir.exists():
        parser.error("output-dir already exists; choose a new directory to preserve previous results")
    source_files = ["scripts/models/__init__.py", "scripts/models/dscript.py",
                    "scripts/models/query_patch.py", "scripts/models/fully_connected.py",
                    "scripts/training_pipeline/execute_pipeline.py",
                    "scripts/training_pipeline/losses.py", "scripts/training_pipeline/metrics.py"]
    source_hashes = {name: sha256(args.source_root / name) for name in source_files}
    input_hashes = {str(path.resolve()): sha256(path) for path in (args.interactions, args.embedding_manifest)}
    if args.submission_record:
        record = json.loads(args.submission_record.read_text())
        for key, values in (("source_sha256", source_hashes), ("input_sha256", input_hashes)):
            for name, digest in values.items():
                if record[key].get(name) != digest:
                    raise ValueError(f"Original {key} mismatch: {name}")

    sys.path.insert(0, str(args.source_root))
    pipeline = importlib.import_module("scripts.training_pipeline.execute_pipeline")
    if Path(pipeline.__file__).resolve() != args.source_root / source_files[4]:
        raise RuntimeError("Imported training pipeline does not match --source-root")
    import pandas as pd
    import torch

    pipeline.pl.seed_everything(args.seed, workers=True)
    torch.set_float32_matmul_precision("high")
    if args.device == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable; run in a GPU allocation")
    device = torch.device(args.device)
    model_class = pipeline.DScriptLightningModule if args.model == "dscript" else pipeline.QueryPatchLightningModule
    model = model_class.load_from_checkpoint(str(args.checkpoint), map_location="cpu")
    # Preserve the checkpoint's parameter gradient flags. inference_mode already
    # disables autograd; additionally calling requires_grad_(False) changes
    # Query Patch's numeric execution on this runtime. A two-path GPU check
    # reproduced original scores exactly only with the original flags intact.
    model.eval().to(device)
    data = pipeline.DScriptDataModule(
        interaction_path=args.interactions, embedding_manifest=args.embedding_manifest,
        contact_map_dir=args.contact_map_dir, batch_size=args.batch_size,
        num_workers=args.num_workers, loss_mode="basic", limit_rows=None, negative_contact_maps=False)
    data.setup("fit")
    provenance = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "model": args.model, "seed": args.seed, "checkpoint": str(args.checkpoint),
        "parameter_gradient_flags": "preserved from checkpoint restoration; autograd disabled by inference_mode",
        "checkpoint_sha256": sha256(args.checkpoint), "model_hyperparameters": dict(model.hparams),
        "source_root": str(args.source_root), "source_sha256": source_hashes, "input_sha256": input_hashes,
        "embedding_manifest": str(args.embedding_manifest.resolve()),
        "contact_map_dir": str(args.contact_map_dir.resolve()),
        "submission_record": str(args.submission_record.resolve()) if args.submission_record else None,
        "submission_record_sha256": sha256(args.submission_record) if args.submission_record else None,
        "analysis_sha256": {name: sha256(Path(__file__).with_name(name)) for name in
                            ("select_contact_threshold.py", "contact_threshold_metrics.py", "select_interaction_threshold.py")},
        "batch_size": args.batch_size, "shuffle": False, "device": str(device),
        "validation_ordering_note": "Ordered post-hoc inference; historical training validation loader shuffled batches",
        "float32_matmul_precision": torch.get_float32_matmul_precision(),
        "torch_version": torch.__version__, "numpy_version": np.__version__,
        "eligibility": "Real structural maps with known cells; -1 unknown/padding and synthetic maps excluded",
        "zero_division": 0, "micro_ap_bins": 4096,
        "historical_contact_map_hashes_available": False,
    }
    args.output_dir.mkdir(parents=True)
    validation_maps, validation_rows, _, exclusions = collect_contacts(
        pipeline, data, model, data.val_dataset, device, args.output_dir, "validation")
    search = ExactContactThresholdSearch()
    for labels, scores in validation_maps:
        search.update(labels, scores)
    threshold, curve, search_metadata = search.select(args.selection_metric)
    validation = summarize_contacts(validation_maps, validation_rows, threshold,
                                    args.output_dir / "validation_maps.tsv")
    aggregation, metric = args.selection_metric.split("_")
    actual_objective = validation["at_selected_threshold"][aggregation][metric]
    if not np.isclose(actual_objective, search_metadata["selected_objective_value"], rtol=0, atol=1e-12):
        raise RuntimeError("Selected search objective disagrees with independent fixed-cutoff evaluation")
    validation["exclusions"] = exclusions
    threshold_record = {**provenance, "threshold": threshold, "selection_split": "val",
                        "selection_metric": args.selection_metric, "search": search_metadata,
                        "tie_break": "closest to 0.5, then smaller threshold",
                        "validation": validation}
    # The immutable validation-only record precedes ALL test inference and
    # metric calculations. Test labels cannot influence the cutoff or checkpoint.
    (args.output_dir / "threshold.json").write_text(json.dumps(threshold_record, indent=2) + "\n")
    with (args.output_dir / "validation_threshold_curve.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(curve[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(curve)
    print(f"Frozen contact threshold {threshold:.17g}; validation macro F1 "
          f"{validation['at_0.5']['macro']['f1']:.6f} -> "
          f"{validation['at_selected_threshold']['macro']['f1']:.6f}", flush=True)
    del search, validation_maps

    test_path = (args.test_predictions or args.checkpoint.parent.parent / "test_predictions.tsv").resolve(strict=True)
    original = pd.read_csv(test_path, sep="\t", dtype=str, keep_default_na=False)
    validate_prediction_cohort(original, data.test_dataset.examples, "test")
    test_maps, test_rows, reproduced, exclusions = collect_contacts(
        pipeline, data, model, data.test_dataset, device, args.output_dir, "test", threshold, args.sample_count)
    expected = pd.to_numeric(original["predicted_probability"]).to_numpy()
    differences = np.abs(reproduced - expected)
    verification = {"n": len(expected), "absolute_tolerance": args.verification_atol,
                    "relative_tolerance": 2e-5, "max_absolute_probability_error": float(differences.max()),
                    "mean_absolute_probability_error": float(differences.mean())}
    original["reproduced_probability"] = reproduced
    original["absolute_probability_error"] = differences
    original.to_csv(args.output_dir / "test_score_verification.tsv", sep="\t", index=False)
    if not np.allclose(reproduced, expected, atol=args.verification_atol, rtol=2e-5):
        raise ValueError(f"Restored checkpoint fails saved interaction-score verification: {verification}")
    test = summarize_contacts(test_maps, test_rows, threshold, args.output_dir / "test_maps.tsv")
    test["exclusions"] = exclusions
    summary = {**threshold_record, "test": test, "test_score_verification": verification,
               "test_predictions": str(test_path), "test_predictions_sha256": sha256(test_path),
               "output_sha256": {path.name: sha256(path) for path in args.output_dir.iterdir() if path.is_file()}}
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Test macro F1 {test['at_0.5']['macro']['f1']:.6f} -> "
          f"{test['at_selected_threshold']['macro']['f1']:.6f}; "
          f"macro contact AP {test['contact_auprc_macro']:.6f}", flush=True)


if __name__ == "__main__":
    main()
