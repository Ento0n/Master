#!/usr/bin/env python3
"""Verify contact-threshold analyses and export thesis test metrics as CSVs.

Each checkpoint chooses one contact cutoff using validation mean per-map F1.
This report verifies frozen cutoffs, exported map-table hashes and confusion
counts, then exports test precision, recall and F1 at that fitted cutoff, plus
cutoff-independent test AUPRC. Every metric is calculated separately for each
structural map and averaged with equal map weight. Unknown and padded cells are
excluded. The per-run CSV retains the fitted cutoff and its validation F1;
the aggregate CSV reports only test-metric means and sample standard deviations
across training seeds. All metrics and cutoffs are on a 0-1 scale.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path


DECISIONS = ("at_0.5", "at_selected_threshold")
METRICS = ("accuracy", "precision", "recall", "f1")
COUNTS = ("tp", "tn", "fp", "fn")
MODELS = ("dscript", "query_patch")
EMBEDDINGS = ("esm2", "prostt5")


# =============================================================================
# Verification: the exported per-map confusion counts are the reporting source
# =============================================================================

def metrics_from_counts(counts: dict[str, int]) -> dict[str, float]:
    """Use positive := score >= cutoff; zero denominators yield zero."""
    tp, tn, fp, fn = (counts[key] for key in COUNTS)
    return {
        "accuracy": (tp + tn) / max(tp + tn + fp + fn, 1),
        "precision": tp / max(tp + fp, 1),
        "recall": tp / max(tp + fn, 1),
        "f1": 2 * tp / max(2 * tp + fp + fn, 1),
    }


def require_close(actual: object, expected: float, context: str) -> None:
    """Allow text/float serialization error without accepting metric changes."""
    if not math.isclose(float(actual), expected, rel_tol=1e-7, abs_tol=1e-9):
        raise ValueError(f"{context}: recorded {actual}, recomputed {expected}")


def verify_saved_interaction_scores(summary: dict, summary_path: Path) -> tuple[dict, list]:
    """Audit every restored interaction score using its recorded tolerances.

    The tolerance is checked separately per score: atol + rtol * abs(original).
    Comparing only the largest absolute difference against atol would reject
    legitimate relative-tolerance matches and cannot verify the other rows.
    """
    verification_path = summary_path.parent / "test_score_verification.tsv"
    original_path = Path(summary["test_predictions"])
    for path, digest in ((verification_path, summary["output_sha256"].get(verification_path.name)),
                         (original_path, summary["test_predictions_sha256"])):
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Interaction-score table hash mismatch: {path}")
    with original_path.open(newline="") as handle:
        originals = list(csv.DictReader(handle, delimiter="\t"))
    with verification_path.open(newline="") as handle:
        reproduced_rows = list(csv.DictReader(handle, delimiter="\t"))
    verification = summary["test_score_verification"]
    if not originals or len(originals) != len(reproduced_rows) or len(originals) != int(verification["n"]):
        raise ValueError(f"Interaction-score verification count mismatch: {summary_path}")
    if len(originals) != int(summary["test"]["exclusions"]["total_pairs"]):
        raise ValueError(f"Interaction verification does not cover every test pair: {summary_path}")
    atol, rtol = (float(verification[key]) for key in ("absolute_tolerance", "relative_tolerance"))
    if any(not math.isfinite(value) or value < 0 for value in (atol, rtol)):
        raise ValueError(f"Invalid interaction-score verification tolerance: {summary_path}")
    errors = []
    for index, (original, reproduced) in enumerate(zip(originals, reproduced_rows)):
        # All original metadata and probability text must survive the export.
        # This also verifies the ordered cohort and source scores independently.
        if any(reproduced.get(key) != value for key, value in original.items()):
            raise ValueError(f"Original test row changed: {verification_path}: row {index}")
        expected = float(original["predicted_probability"])
        observed = float(reproduced["reproduced_probability"])
        if any(not math.isfinite(value) or not 0 <= value <= 1 for value in (expected, observed)):
            raise ValueError(f"Invalid interaction probability: {verification_path}: row {index}")
        error = abs(observed - expected)
        if error > atol + rtol * abs(expected):
            raise ValueError(f"Interaction-score reproduction exceeds recorded tolerance: {verification_path}: row {index}")
        # The inference CLI parsed expected probabilities with pandas, whose
        # decimal conversion can differ from Python's by one float64 ULP.
        # Allow that serialization error when checking stored error summaries.
        require_close(reproduced["absolute_probability_error"], error, f"{verification_path}: row {index} error")
        errors.append(error)
    require_close(verification["max_absolute_probability_error"], max(errors), f"{verification_path}: maximum error")
    require_close(verification["mean_absolute_probability_error"], statistics.fmean(errors), f"{verification_path}: mean error")
    result = {
        "score_verification_n": len(originals),
        "score_verification_max_absolute_error": float(verification["max_absolute_probability_error"]),
        "score_verification_mean_absolute_error": float(verification["mean_absolute_probability_error"]),
        "score_verification_absolute_tolerance": atol,
        "score_verification_relative_tolerance": rtol,
    }
    cohort = [tuple(row[key] for key in ("split", "new_cluster_pair", "cluster_1", "cluster_2",
                                        "pdb_id", "assembly_number", "dimer_type", "label")) for row in originals]
    return result, cohort


def load_verified_run(summary_path: Path) -> tuple[dict, dict]:
    """Validate both splits and return verified metrics plus cohort identities."""
    summary = json.loads(summary_path.read_text())
    frozen_path = summary_path.parent / "threshold.json"
    frozen = json.loads(frozen_path.read_text())
    if summary["output_sha256"].get("threshold.json") != hashlib.sha256(frozen_path.read_bytes()).hexdigest():
        raise ValueError(f"Frozen threshold hash mismatch: {frozen_path}")
    # Archived analyses retain their original keys; CSV output names use the
    # project's single convention of averaging metrics equally across maps.
    if summary.get("selection_metric") != "macro_f1" or summary.get("selection_split") != "val":
        raise ValueError(f"Expected a validation mean per-map F1 cutoff: {summary_path}")
    if "test" in frozen:
        raise ValueError(f"Frozen threshold record unexpectedly contains test metrics: {summary_path}")
    # The frozen file may predate validation output hashes, so compare every
    # field it does contain against the final summary, including provenance.
    for key, value in frozen.items():
        if key not in summary or summary[key] != value:
            raise ValueError(f"Frozen threshold field changed: {summary_path}: {key}")
    for key in ("threshold", "selection_metric", "selection_split", "checkpoint_sha256", "validation"):
        if key not in frozen:
            raise ValueError(f"Missing frozen threshold field: {summary_path}: {key}")
    threshold = float(summary["threshold"])
    if not math.isfinite(threshold) or not 0 <= threshold <= math.nextafter(1.0, math.inf):
        raise ValueError(f"Invalid selected threshold: {summary_path}")

    # Final-run directories identify the embedding independently of optional
    # provenance fields; the validated checkpoint belongs to the same run.
    checkpoint = Path(summary["checkpoint"])
    embeddings = [name for name in EMBEDDINGS if name in checkpoint.parts]
    if len(embeddings) != 1 or summary["model"] not in MODELS:
        raise ValueError(f"Unrecognized final checkpoint/model: {summary_path}")
    embedding = embeddings[0]
    row = {
        "embedding": embedding, "model": summary["model"], "seed": int(summary["seed"]),
        "threshold": threshold, "checkpoint": str(checkpoint),
        "checkpoint_sha256": summary["checkpoint_sha256"],
        "summary": str(summary_path.resolve()),
    }
    cohorts = {}
    for split in ("validation", "test"):
        table_name = f"{split}_maps.tsv"
        table_path = summary_path.parent / table_name
        expected_hash = summary["output_sha256"].get(table_name)
        actual_hash = hashlib.sha256(table_path.read_bytes()).hexdigest()
        if expected_hash != actual_hash:
            raise ValueError(f"Map-table hash mismatch: {table_path}")
        with table_path.open(newline="") as handle:
            maps = list(csv.DictReader(handle, delimiter="\t"))
        if not maps:
            raise ValueError(f"No evaluable structural maps: {table_path}")
        # Split indices preserve original cohort order, including gaps where
        # synthetic negatives or maps without known cells were excluded.
        indices = [int(item["split_index"]) for item in maps]
        if len(set(indices)) != len(indices) or indices != sorted(indices):
            raise ValueError(f"Duplicate or unordered map indices: {table_path}")
        cohorts[split] = [tuple(item.get(key, "") for key in
                               ("split_index", "dimer_type", "pdb_id", "assembly_number",
                                "new_cluster_pair", "contact_map_sha256",
                                "known_cells", "positive_cells")) for item in maps]
        known_cells = sum(int(item["known_cells"]) for item in maps)
        positive_cells = sum(int(item["positive_cells"]) for item in maps)
        split_summary = summary[split]
        for decision in DECISIONS:
            metric_values = {key: [] for key in METRICS}
            for item in maps:
                counts = {key: int(item[f"{decision}_{key}"]) for key in COUNTS}
                known, positive = int(item["known_cells"]), int(item["positive_cells"])
                if (known <= 0 or min(counts.values()) < 0 or sum(counts.values()) != known
                        or counts["tp"] + counts["fn"] != positive):
                    raise ValueError(f"Invalid map confusion counts: {table_path}: {item['split_index']}")
                metrics = metrics_from_counts(counts)
                for key in METRICS:
                    require_close(item[f"{decision}_{key}"], metrics[key], f"{table_path}: {decision}/{key}")
                    metric_values[key].append(metrics[key])
            recorded = split_summary[decision]
            for key, value in (("n_maps", len(maps)), ("known_cells", known_cells),
                               ("positive_cells", positive_cells)):
                if int(recorded[key]) != value:
                    raise ValueError(f"Aggregate count mismatch: {table_path}: {decision}/{key}")
                row[f"{split}_{key}"] = value
            # Average each map's metric, so large maps do not receive extra
            # weight. Keep the legacy input key only for result verification.
            for key, values in metric_values.items():
                value = statistics.fmean(values)
                require_close(recorded["macro"][key], value, f"{table_path}: {decision}/{key}")
                row[f"{split}_{decision}_{key}"] = value
        mean_ap = statistics.fmean(float(item["auprc"]) for item in maps)
        contact_ap = float(split_summary["contact_auprc_macro"])
        require_close(contact_ap, mean_ap, f"{table_path}: mean per-map AP")
        if not math.isfinite(contact_ap) or not 0 <= contact_ap <= 1:
            raise ValueError(f"Invalid contact AP: {table_path}")
        row[f"{split}_contact_auprc"] = contact_ap
    if row["validation_at_selected_threshold_f1"] + 1e-9 < row["validation_at_0.5_f1"]:
        raise ValueError(f"Selected validation F1 is below the included 0.5 baseline: {summary_path}")
    verification, cohorts["test_interactions"] = verify_saved_interaction_scores(summary, summary_path)
    row.update(verification)
    return row, cohorts


# =============================================================================
# Reporting: seeds are repeated training runs on the same held-out map cohort
# =============================================================================

def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", type=Path, action="append", required=True,
                        help="Directory containing completed contact analyses; repeat to combine model/embedding arrays")
    parser.add_argument("--output-prefix", type=Path, required=True,
                        help="Output prefix for the per-run .csv and _aggregate.csv files")
    args = parser.parse_args()
    analysis_roots = [root.resolve(strict=True) for root in args.analysis_root]
    summaries = []
    for root in analysis_roots:
        found = sorted(root.rglob("summary.json")) if root.is_dir() else []
        if not found:
            raise ValueError(f"Each requested analysis root must contain at least one completed run summary: {root}")
        summaries.extend(found)
    rows, reference_cohorts = [], None
    for path in summaries:
        row, cohorts = load_verified_run(path)
        if reference_cohorts is not None and cohorts != reference_cohorts:
            raise ValueError(f"Structural map cohort differs across final runs: {path}")
        reference_cohorts = cohorts
        rows.append(row)
    # Require all six model/seed combinations for each included embedding.
    # Reject duplicate roots from separate analyses of the same checkpoint.
    combinations = [(row["embedding"], row["model"], row["seed"]) for row in rows]
    embeddings = sorted({row["embedding"] for row in rows})
    expected = {(embedding, model, seed) for embedding in embeddings
                for model in MODELS for seed in (42, 123, 456)}
    if len(combinations) != len(set(combinations)) or set(combinations) != expected:
        raise ValueError("Each embedding requires exactly D-SCRIPT/QP-Conv-IC × seeds 42/123/456")
    rows.sort(key=lambda row: (row["embedding"], row["model"], row["seed"]))

    # Verification above still covers both splits and cutoffs. Keep the one
    # validation score used to select each cutoff alongside the thesis test
    # results; AP uses continuous scores and is independent of the cutoff.
    report_rows = [{
        "embedding": row["embedding"], "model": row["model"], "seed": row["seed"],
        "threshold": row["threshold"],
        "validation_f1": row["validation_at_selected_threshold_f1"],
        "test_auprc": row["test_contact_auprc"],
        "test_precision": row["test_at_selected_threshold_precision"],
        "test_recall": row["test_at_selected_threshold_recall"],
        "test_f1": row["test_at_selected_threshold_f1"],
    } for row in rows]

    # Each row summarizes repeated training seeds on the same test cohort.
    # Cutoffs remain per-run settings, so only test metrics are aggregated.
    aggregates = []
    for embedding in embeddings:
        for model in MODELS:
            group = [row for row in report_rows if row["embedding"] == embedding and row["model"] == model]
            aggregate = {"embedding": embedding, "model": model, "n_seeds": len(group)}
            for key in ("test_auprc", "test_precision", "test_recall", "test_f1"):
                values = [row[key] for row in group]
                aggregate[f"{key}_mean"] = statistics.fmean(values)
                aggregate[f"{key}_sample_sd"] = statistics.stdev(values)
            aggregates.append(aggregate)

    prefix = args.output_prefix.resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = prefix.with_suffix(".csv")
    aggregate_path = prefix.with_name(prefix.name + "_aggregate.csv")
    write_csv(csv_path, report_rows)
    write_csv(aggregate_path, aggregates)
    print(f"Verified {len(rows)} runs; wrote {csv_path} and {aggregate_path}")


if __name__ == "__main__":
    main()
