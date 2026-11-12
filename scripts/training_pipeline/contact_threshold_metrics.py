"""Fit one validation contact cutoff, then evaluate it on known residue cells.

Each ``update`` receives a single structural map flattened after removing padding
and unknown labels. Each metric weights nonempty maps equally. An all-negative
structural map is real evaluation data and has F1/AP zero. Empty maps contribute
neither cells nor an averaging denominator.

The exact search stores one event per distinct score within each map, not a
Python object per residue cell. Grouped score events update global confusion
counts and the sum of per-map F1 values. Consequently every attainable decision
is considered before the reported curve is thinned, even for large datasets.
Only NumPy is required; model inference and validation/test separation belong to
the caller.
"""

from __future__ import annotations

import numpy as np


# =============================================================================
# Known-cell validation and metrics at a threshold that is already frozen
# =============================================================================

def _validate_known_cells(labels, scores) -> tuple[np.ndarray, np.ndarray]:
    """Return matching 1-D binary targets and float64 probabilities in [0, 1]."""
    labels = np.asarray(labels)
    scores = np.asarray(scores, dtype=np.float64)
    if labels.ndim != 1 or scores.ndim != 1 or labels.shape != scores.shape:
        raise ValueError("Labels and scores must be matching one-dimensional arrays")
    if not np.all((labels == 0) | (labels == 1)):
        raise ValueError("Known contact labels must contain only 0 and 1")
    if not np.all(np.isfinite(scores)) or np.any((scores < 0) | (scores > 1)):
        raise ValueError("Contact scores must be finite probabilities in [0, 1]")
    return labels.astype(np.int8, copy=False), scores


def _metrics_from_counts(tp: int, tn: int, fp: int, fn: int) -> dict:
    """Use zero for undefined precision/recall/F1 and for empty-map accuracy."""
    n = tp + tn + fp + fn
    return {
        "accuracy": (tp + tn) / n if n else 0.0,
        "precision": tp / (tp + fp) if tp + fp else 0.0,
        "recall": tp / (tp + fn) if tp + fn else 0.0,
        "f1": 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else 0.0,
        "tp": tp, "tn": tn, "fp": fp, "fn": fn,
        "known_cells": n, "positive_cells": tp + fn,
        "predicted_positive": tp + fp,
    }


def contact_metrics(labels, scores, threshold: float) -> dict:
    """Evaluate one known-cell map using positive := float64(score) >= cutoff.

    Float64 comparison is deliberate: a cutoff between adjacent float32 model
    outputs must not be rounded back onto the lower output. A finite cutoff
    slightly above 1 is valid and represents the all-negative decision when an
    original probability equals 1. No threshold is fitted by this function.
    """
    labels, scores = _validate_known_cells(labels, scores)
    if not np.isfinite(threshold):
        raise ValueError("The contact threshold must be finite")
    predicted = scores >= float(threshold)
    positive = labels == 1
    tp = int(np.count_nonzero(predicted & positive))
    fp = int(np.count_nonzero(predicted & ~positive))
    fn = int(np.count_nonzero(~predicted & positive))
    tn = int(labels.size) - tp - fp - fn
    return _metrics_from_counts(tp, tn, fp, fn)


def aggregate_contact_metrics(per_map_metrics) -> dict:
    """Average map-level ratios equally and retain summed cell counts.

    Empty maps are excluded even if the caller passes their zero-valued metric
    records. All-negative nonempty maps remain in the averaging denominator.
    The input may be a generator, allowing inference to discard each full map
    immediately. Cell counts describe the cohort; they do not set metric weights.
    """
    totals = {key: 0 for key in ("tp", "tn", "fp", "fn")}
    metric_sums = {key: 0.0 for key in ("accuracy", "precision", "recall", "f1")}
    n_maps = 0
    for metrics in per_map_metrics:
        if not metrics["known_cells"]:
            continue
        n_maps += 1
        for key in totals:
            totals[key] += int(metrics[key])
        for key in metric_sums:
            metric_sums[key] += float(metrics[key])
    if not n_maps:
        raise ValueError("No maps with known contact cells are available")
    return {
        **{key: value / n_maps for key, value in metric_sums.items()},
        **totals,
        "n_maps": n_maps,
        "known_cells": sum(totals.values()),
        "positive_cells": totals["tp"] + totals["fn"],
        "predicted_positive": totals["tp"] + totals["fp"],
    }


def contact_average_precision(labels, scores) -> float:
    """Compute exact non-interpolated ranking AP, grouping equal scores.

    This threshold-independent value is zero for a map with no positive known
    labels (including empty maps). Exclude empty maps when averaging over maps.
    """
    labels, scores = _validate_known_cells(labels, scores)
    positives = int(labels.sum())
    if not positives:
        return 0.0
    order = np.argsort(-scores, kind="stable")
    sorted_scores = scores[order]
    ends = np.r_[np.flatnonzero(sorted_scores[1:] != sorted_scores[:-1]), labels.size - 1]
    true_positives = np.cumsum(labels[order], dtype=np.int64)[ends]
    added_positives = np.diff(np.r_[0, true_positives])
    return float(np.sum((added_positives / positives) * (true_positives / (ends + 1))))


# =============================================================================
# Exact validation search using one removal event per distinct map score
# =============================================================================

class ExactContactThresholdSearch:
    """Accumulate structural maps and select one global contact threshold.

    Memory is O(number of distinct within-map scores), approximately 40 bytes
    per stored event on platforms with extended-precision NumPy long doubles.
    Selection temporarily allocates sorting/grouping arrays. There are no Python
    loops over individual residue cells or candidate cutoffs. Extended precision
    limits accumulated rounding in F1 deltas across many large maps.
    """

    def __init__(self) -> None:
        self._events = []
        self._initial_f1_sum = np.longdouble(0)
        self.n_maps = 0
        self.excluded_empty_maps = 0
        self.known_cells = 0
        self.positive_cells = 0

    def update(self, labels, scores) -> None:
        """Add one map; the caller must first mask unknown and padded cells."""
        labels, scores = _validate_known_cells(labels, scores)
        if not labels.size:
            self.excluded_empty_maps += 1
            return
        n = int(labels.size)
        positives = int(labels.sum())
        order = np.argsort(scores, kind="stable")
        sorted_scores = scores[order]
        starts = np.r_[0, np.flatnonzero(sorted_scores[1:] != sorted_scores[:-1]) + 1]
        removed_positive = np.add.reduceat(labels[order].astype(np.int64), starts)
        block_sizes = np.diff(np.r_[starts, n])
        removed_negative = block_sizes - removed_positive

        # Increasing the cutoff past a tied score removes its entire block from
        # predicted positives. Per-map F1 = 2*TP / (true positives + predicted
        # positives), so its change requires only these cumulative counts.
        remaining_tp = positives - np.cumsum(removed_positive)
        remaining_predicted = n - np.cumsum(block_sizes)
        denominator = positives + remaining_predicted
        after_f1 = np.zeros(starts.size, dtype=np.longdouble)
        np.divide(2 * remaining_tp.astype(np.longdouble), denominator,
                  out=after_f1, where=denominator != 0)
        initial_f1 = np.longdouble(2 * positives) / (positives + n)
        f1_delta = np.diff(np.r_[initial_f1, after_f1])
        self._events.append((sorted_scores[starts], removed_positive, removed_negative, f1_delta))
        self._initial_f1_sum += initial_f1
        self.n_maps += 1
        self.known_cells += n
        self.positive_cells += positives

    def select(self, objective: str = "f1", max_curve_points: int = 2001) -> tuple[float, list[dict], dict]:
        """Optimize every score split; thin only the exported diagnostic curve.

        Candidates match the interaction-threshold convention: 0 (all positive),
        midpoints of adjacent distinct scores, nextafter(max_score, +inf) (all
        negative), and 0.5. With adjacent floating-point scores, the upper score
        replaces an unrepresentable midpoint under the >= comparison rule.
        Ties prefer the candidate closest to 0.5, then the smaller candidate.
        F1 ties allow only extended-precision roundoff (64 machine eps).

        At most ``max_curve_points`` rows are returned, always including both
        boundaries, the selected cutoff, and 0.5. Selection requires at least
        one positive known contact in the complete validation cohort.
        """
        if objective != "f1":
            raise ValueError("Objective must be f1")
        if not isinstance(max_curve_points, (int, np.integer)) or max_curve_points < 4:
            raise ValueError("max_curve_points must be an integer of at least 4")
        if not self.n_maps or not self.positive_cells:
            raise ValueError("Threshold selection requires maps with positive known contact labels")

        # Each column is an array, not millions of Python records. Global score
        # grouping prevents equal probabilities from receiving different labels,
        # including equal scores belonging to different structural maps.
        scores, removed_positive, removed_negative, f1_delta = (
            np.concatenate([event[column] for event in self._events]) for column in range(4)
        )
        order = np.argsort(scores, kind="stable")
        scores = scores[order]
        starts = np.r_[0, np.flatnonzero(scores[1:] != scores[:-1]) + 1]
        scores = scores[starts]
        removed_positive = np.add.reduceat(removed_positive[order], starts)
        removed_negative = np.add.reduceat(removed_negative[order], starts)
        f1_delta = np.add.reduceat(f1_delta[order], starts)
        del order, starts

        # State index k means that precisely the k lowest distinct score blocks
        # have been removed. The zero state and final state are constant labels.
        fn = np.r_[0, np.cumsum(removed_positive)]
        tn = np.r_[0, np.cumsum(removed_negative)]
        f1 = np.r_[self._initial_f1_sum, self._initial_f1_sum + np.cumsum(f1_delta)] / self.n_maps
        f1[-1] = 0
        np.clip(f1, 0, 1, out=f1)
        del removed_positive, removed_negative, f1_delta
        thresholds = np.empty(scores.size + 1, dtype=np.float64)
        thresholds[0] = 0.0
        thresholds[-1] = np.nextafter(scores[-1], np.inf)
        thresholds[1:-1] = scores[:-1] + (scores[1:] - scores[:-1]) / 2
        thresholds[1:-1] = np.where(thresholds[1:-1] > scores[:-1], thresholds[1:-1], scores[1:])
        default_state = int(np.searchsorted(scores, 0.5, side="left"))
        default_is_candidate = bool(np.any(thresholds == 0.5))

        # F1 uses the accumulated changes in each map's ratio. The explicit
        # default may represent an already covered state but is an additional
        # candidate because the distance-to-0.5 tie break matters.
        best_value = f1.max()
        tie_atol = 64 * np.finfo(np.longdouble).eps
        tied = np.flatnonzero(np.abs(f1 - best_value) <= tie_atol)
        # Candidates are sorted, so argmin also implements the smaller-cutoff
        # preference for two candidates equally distant from 0.5.
        best_state = int(tied[np.argmin(np.abs(thresholds[tied] - 0.5))])
        threshold = float(thresholds[best_state])
        if abs(f1[default_state] - best_value) <= tie_atol:
            best_state, threshold = default_state, 0.5

        # Keep the exact winner and boundaries while limiting curve serialization
        # to a few thousand Python dictionaries. Reserve a row for 0.5 even when
        # it lies inside an existing decision interval.
        required = {0, thresholds.size - 1}
        if threshold != 0.5:
            required.add(best_state)
        default_index = int(np.searchsorted(thresholds, 0.5))
        if default_is_candidate:
            required.add(default_index)
        budget = max_curve_points - (not default_is_candidate)
        if thresholds.size <= budget:
            sampled = np.arange(thresholds.size)
        else:
            extra = np.linspace(0, thresholds.size - 1, budget - len(required), dtype=np.int64)
            sampled = np.unique(np.r_[list(required), extra])

        def curve_row(state: int, cutoff: float) -> dict:
            false_negative, true_negative = int(fn[state]), int(tn[state])
            true_positive = self.positive_cells - false_negative
            false_positive = self.known_cells - self.positive_cells - true_negative
            return {"threshold": float(cutoff), "f1": float(f1[state]),
                    "tp": true_positive, "tn": true_negative, "fp": false_positive, "fn": false_negative}

        curve = [curve_row(int(state), thresholds[state]) for state in sampled]
        if not default_is_candidate:
            curve.append(curve_row(default_state, 0.5))
            curve.sort(key=lambda row: row["threshold"])
        metadata = {
            "objective": objective, "selected_objective_value": float(f1[best_state]),
            "n_maps": self.n_maps, "excluded_empty_maps": self.excluded_empty_maps,
            "known_cells": self.known_cells, "positive_cells": self.positive_cells,
            "unique_scores": int(scores.size),
            "evaluated_thresholds": int(thresholds.size + (not default_is_candidate)),
            "reported_thresholds": len(curve),
            "selection_is_exact": True,
            "tie_break": "closest to 0.5, then smaller threshold",
            "objective_tie_atol": float(tie_atol),
        }
        return threshold, curve, metadata
