#!/usr/bin/env python3
"""Compare additive and final QUERY-PATCH on three previously selected maps.

Read the frozen ProstT5 seed-42 reference (1QC9, 1SEZ, 7QAA) and the completed additive
checkpoint's contact-threshold analysis. Match maps by identity, reconstruct
binary predictions at the validation-selected cutoff, and verify metrics before
exporting full contact maps. This performs no inference, training, or example
selection. Weak/Moderate/Strong describe the FINAL model's existing selection;
they do not rank the additive model or establish a general model advantage.
Only PNGs go to figures/; reference inputs, PDF, and metadata live with the run.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

# Direct execution puts only scripts/figures on the import path. Add the
# repository root so the shared training helpers also resolve outside the repo.
if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.training_pipeline import contact_threshold_metrics

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ANALYSIS = (REPO_ROOT / "jobs/training_pipeline/results/prostt5_additive_contact_figure"
                    / "7171489/qp_additive_normal_seed_42/contact_threshold_selection")
DEFAULT_ARTIFACTS = DEFAULT_ANALYSIS.parent / "figure_1sez"
DEFAULT_REFERENCE = DEFAULT_ARTIFACTS / "reference"


# =============================================================================
# Input integrity, identity matching, and independently checked map metrics
# =============================================================================

def sha256(path: Path) -> str:
    """Hash large prediction caches without reading them all into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_metrics(actual: dict, record: dict, prefix: str = "") -> None:
    """Check exact confusion counts and full-precision scalar metric exports."""
    for key, value in actual.items():
        expected = float(record[prefix + key])
        matches = value == expected if isinstance(value, int) else np.isclose(
            value, expected, rtol=0, atol=1e-12)
        if not matches:
            raise ValueError(f"Metric mismatch for {prefix + key}: {value} != {expected}")


def load_examples(analysis: Path, reference: Path) -> tuple[list, dict]:
    """Validate sources and return complete [-1,0,1] maps plus their metrics.

    The additive cache contains flattened known cells only, in row-major order.
    Its offsets are indexed by split_indices, independently of the TSV order.
    Frozen full maps restore shape and unknown cells (-1); unknown cells are
    masked identically in all three panels and never enter AP or F1.
    """
    input_hashes = {}
    manifest_path = reference / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    expected_files = {"selected_examples.csv", "provenance.json",
                      "1QC9_assembly1.npz", "1SEZ_assembly1.npz", "7QAA_assembly1.npz"}
    recorded_hashes = manifest["file_sha256"]
    if not expected_files.issubset(recorded_hashes):
        raise ValueError("Reference manifest must hash the selection, provenance, and all three maps")
    for name, expected_hash in recorded_hashes.items():
        path = (reference / name).resolve(strict=True)
        if Path(name).is_absolute() or not path.is_relative_to(reference):
            raise ValueError(f"Reference manifest entry escapes its directory: {name}")
        actual_hash = sha256(path)
        if actual_hash != expected_hash:
            raise ValueError(f"Frozen reference hash mismatch: {name}")
        input_hashes[str(path)] = actual_hash
    input_hashes[str(manifest_path)] = sha256(manifest_path)
    reference_provenance = json.loads((reference / "provenance.json").read_text())
    if reference_provenance["embedding"] != "ProstT5" or reference_provenance["seed"] != 42:
        raise ValueError("The reference must use ProstT5 embeddings and seed 42")
    with (reference / "selected_examples.csv").open(newline="") as handle:
        selected = [row for row in csv.DictReader(handle)
                    if row["model"] == "query_patch" and int(row["seed"]) == 42]
    expected_order = [("Weak", "1QC9", 1), ("Moderate", "1SEZ", 1), ("Strong", "7QAA", 1)]
    if [(row["selection_role"], row["pdb_id"], int(row["assembly"]))
            for row in selected] != expected_order:
        raise ValueError("Reference selection must use Weak 1QC9, Moderate 1SEZ, Strong 7QAA")
    if len({int(row["split_index"]) for row in selected}) != 3:
        raise ValueError("Reference split indices must be unique")

    summary_path = analysis / "summary.json"
    summary = json.loads(summary_path.read_text())
    input_hashes[str(summary_path)] = sha256(summary_path)
    hyperparameters = summary["model_hyperparameters"]
    required_hyperparameters = {"query_evidence_mode": "additive",
                                "compatibility_head": "bilinear", "embedding_dim": 1024}
    if (summary["model"] != "query_patch" or summary["seed"] != 42
            or any(hyperparameters.get(key) != value
                   for key, value in required_hyperparameters.items())
            or "prostt5" not in summary["embedding_manifest"].lower()):
        raise ValueError("Analysis must describe the ProstT5 seed-42 additive/bilinear model")
    if summary["selection_split"] != "val" or summary["selection_metric"] != "macro_f1":
        raise ValueError("Additive cutoff must be selected on validation macro F1")
    threshold = float(summary["threshold"])
    if not np.isfinite(threshold):
        raise ValueError("The validation-selected threshold must be finite")
    for name in ("test_maps.tsv", "test_contacts.npz"):
        path = analysis / name
        actual_hash = sha256(path)
        if actual_hash != summary["output_sha256"][name]:
            raise ValueError(f"Additive analysis hash mismatch: {name}")
        input_hashes[str(path)] = actual_hash
    with (analysis / "test_maps.tsv").open(newline="") as handle:
        table = list(csv.DictReader(handle, delimiter="\t"))
    rows_by_index = {int(row["split_index"]): row for row in table}
    if len(rows_by_index) != len(table) or any(row["split"] != "test" for row in table):
        raise ValueError("Additive table must contain unique test split indices")
    with np.load(analysis / "test_contacts.npz", allow_pickle=False) as cache:
        labels, scores = cache["labels"], cache["probabilities"]
        offsets, indices = cache["offsets"], cache["split_indices"]
    if (labels.ndim != 1 or scores.shape != labels.shape or indices.ndim != 1
            or offsets.ndim != 1 or offsets.size != indices.size + 1
            or not np.issubdtype(offsets.dtype, np.integer)
            or not np.issubdtype(indices.dtype, np.integer)
            or offsets[0] != 0 or offsets[-1] != labels.size
            or np.any(np.diff(offsets) <= 0)):
        raise ValueError("Malformed additive known-cell cache shapes or offsets")
    cache_positions = {int(index): position for position, index in enumerate(indices)}
    if len(cache_positions) != len(indices) or set(cache_positions) != set(rows_by_index):
        raise ValueError("Additive cache and table must contain identical unique split indices")
    for split_index, position in cache_positions.items():
        if offsets[position + 1] - offsets[position] != int(rows_by_index[split_index]["known_cells"]):
            raise ValueError(f"Cache/table known-cell count mismatch for split index {split_index}")

    examples = []
    for selected_row in selected:
        split_index = int(selected_row["split_index"])
        row = rows_by_index[split_index]
        identity = (selected_row["pdb_id"], int(selected_row["assembly"]))
        if (row["pdb_id"], int(row["assembly_number"])) != identity:
            raise ValueError(f"Additive/reference identity mismatch at split index {split_index}")
        map_name = f"{identity[0]}_assembly{identity[1]}.npz"
        with np.load(reference / map_name, allow_pickle=False) as cache:
            truth, final_binary = cache["labels"], cache["query_patch"]
        if (truth.ndim != 2 or not all(truth.shape) or final_binary.shape != truth.shape
                or not np.all(np.isin(truth, [-1, 0, 1]))
                or not np.all(np.isin(final_binary, [-1, 0, 1]))
                or not np.array_equal(truth < 0, final_binary < 0)):
            raise ValueError(f"Malformed full reference map or unknown-cell mask: {map_name}")
        if truth.shape != (int(selected_row["length_1"]), int(selected_row["length_2"])):
            raise ValueError(f"Reference CSV/map shape mismatch: {map_name}")
        # Validate the original ground truth's serialized hash as well as its
        # actual labels. NPZ container hashes cannot replace the source NPY hash.
        truth_path = Path(row["contact_map_path"]).resolve(strict=True)
        input_hashes[str(truth_path)] = sha256(truth_path)
        if input_hashes[str(truth_path)] != row["contact_map_sha256"]:
            raise ValueError(f"Source contact-map hash mismatch: {truth_path}")
        if not np.array_equal(np.load(truth_path, allow_pickle=False), truth):
            raise ValueError(f"Additive analysis and reference have different ground truth: {map_name}")
        known = truth >= 0
        position = cache_positions[split_index]
        left, right = offsets[position:position + 2]
        known_labels = labels[left:right]
        known_scores = scores[left:right].astype(np.float64)
        if not np.array_equal(known_labels, truth[known]):
            raise ValueError(f"Additive cache labels do not match full reference map: {map_name}")
        additive_metrics = contact_threshold_metrics.contact_metrics(known_labels, known_scores, threshold)
        verify_metrics(additive_metrics, row, "at_selected_threshold_")
        additive_metrics["auprc"] = contact_threshold_metrics.contact_average_precision(known_labels, known_scores)
        verify_metrics({"auprc": additive_metrics["auprc"],
                        "known_cells": additive_metrics["known_cells"],
                        "positive_cells": additive_metrics["positive_cells"]}, row)
        additive_binary = np.full(truth.shape, -1, dtype=np.int8)
        additive_binary[known] = (known_scores >= threshold).astype(np.int8)

        # Final F1/counts are independently recoverable from frozen binary maps.
        # Final AP requires continuous scores, so preserve the original numeric
        # CSV value and cross-check it against the hashed reference provenance.
        final_metrics = contact_threshold_metrics.contact_metrics(truth[known], final_binary[known], 0.5)
        verify_metrics(final_metrics, selected_row)
        final_metrics["auprc"] = float(selected_row["auprc"])
        final_threshold = float(selected_row["threshold"])
        if (not 0 <= final_metrics["auprc"] <= 1 or not np.isfinite(final_threshold)
                or final_threshold != float(reference_provenance["thresholds"]["query_patch"])):
            raise ValueError("Invalid final reference AP or inconsistent final threshold")
        reference_records = [record for record in reference_provenance["examples"]
                             if record["model"] == "query_patch" and int(record["seed"]) == 42
                             and int(record["split_index"]) == split_index]
        if len(reference_records) != 1:
            raise ValueError(f"Missing or duplicate reference provenance: {map_name}")
        verify_metrics(final_metrics, reference_records[0])
        if (int(selected_row["unknown_cells"]) != int((~known).sum())
                or not np.isclose(float(selected_row["known_fraction"]), known.mean(), rtol=0, atol=1e-12)):
            raise ValueError(f"Reference annotation coverage mismatch: {map_name}")
        examples.append({"selection": selected_row, "truth": truth,
                         "predictions": {"qp_additive": additive_binary, "qp_conv_ic": final_binary},
                         "metrics": {"qp_additive": additive_metrics, "qp_conv_ic": final_metrics},
                         "thresholds": {"qp_additive": threshold, "qp_conv_ic": final_threshold},
                         "contact_map_sha256": row["contact_map_sha256"]})
    provenance = {
        "embedding": "ProstT5", "seed": 42,
        "selection": reference_provenance["selection"],
        "selection_note": "The examples match the current contact-map comparison: 1QC9, 1SEZ, 7QAA. "
                          "The middle example was updated to 1SEZ after additive training. "
                          "Weak/Moderate/Strong refer to FINAL QUERY-PATCH Conv-IC F1, not additive F1. "
                          "These illustrative examples do not estimate an average model advantage.",
        "reference_selection_note": reference_provenance["selection_note"],
        "display": "Complete maps with square cells and no crops or zooms; red=contact (1), "
                   "white=no contact (0), gray=unannotated (-1). Unknown cells are excluded from metrics.",
        "threshold_selection": "Each model uses its own frozen validation-macro-F1 threshold.",
        "thresholds": examples[0]["thresholds"],
        "final_ap_note": "Final AP is preserved from the frozen CSV/provenance; binary-map counts and F1 "
                         "are independently recomputed. Additive AP is recomputed from continuous scores.",
        "additive_checkpoint": summary["checkpoint"],
        "additive_checkpoint_sha256": summary["checkpoint_sha256"],
        "additive_model_hyperparameters": hyperparameters,
        "reference_source_summaries": reference_provenance["source_summaries"],
        "input_sha256": input_hashes,
    }
    return examples, provenance


# =============================================================================
# Full-extent red/white/gray panels and export metadata
# =============================================================================

def render_figure(examples: list, output: Path, artifacts: Path) -> int:
    """Preserve full maps; export PNGs separately from the PDF and metadata."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.patches import Patch

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    gray, white, red = "#bdbdbd", "#ffffff", "#d00000"
    cmap = ListedColormap([gray, white, red])
    norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], 3)
    panel_width, column_gap = 2.0, 0.18
    left_gutter, right_margin = 0.70, 0.08
    header_height, metrics_height, row_gap, footer_height = 0.35, 0.24, 0.18, 0.55
    row_heights = [panel_width * item["truth"].shape[0] / item["truth"].shape[1] for item in examples]
    figure_width = left_gutter + 3 * panel_width + 2 * column_gap + right_margin
    figure_height = header_height + sum(row_heights) + 3 * metrics_height + 2 * row_gap + footer_height
    figure = plt.figure(figsize=(figure_width, figure_height), facecolor="white")
    titles = ("Ground truth", "QUERY-PATCH Additive", "QUERY-PATCH Conv-IC")
    for column, title in enumerate(titles):
        center = left_gutter + column * (panel_width + column_gap) + panel_width / 2
        figure.text(center / figure_width, 1 - 0.10 / figure_height, title,
                    ha="center", va="top", fontsize=9, fontweight="bold")
    top = figure_height - header_height
    for example, height in zip(examples, row_heights):
        truth, row = example["truth"], example["selection"]
        bottom = top - height
        figure.text(0.035 / figure_width, (top - 0.01) / figure_height,
                    f"{row['selection_role']}\n{row['pdb_id']}", ha="left", va="top", fontsize=9)
        matrices = (truth, example["predictions"]["qp_additive"], example["predictions"]["qp_conv_ic"])
        for column, matrix in enumerate(matrices):
            left = left_gutter + column * (panel_width + column_gap)
            axis = figure.add_axes([left / figure_width, bottom / figure_height,
                                    panel_width / figure_width, height / figure_height])
            axis.imshow(matrix, cmap=cmap, norm=norm, interpolation="nearest", aspect="equal", rasterized=True)
            axis.set_xticks([])
            axis.set_yticks([])
            for spine in axis.spines.values():
                spine.set_color("#999999")
                spine.set_linewidth(0.4)
            if column == 0:
                annotation = f"{truth.shape[0]} × {truth.shape[1]} residues"
            else:
                metrics = example["metrics"][("qp_additive", "qp_conv_ic")[column - 1]]
                annotation = f"F1 {metrics['f1']:.3f}   ·   AUPRC {metrics['auprc']:.3f}"
            figure.text((left + panel_width / 2) / figure_width, (bottom - 0.055) / figure_height,
                        annotation, ha="center", va="top", fontsize=8)
        top = bottom - metrics_height - row_gap
    figure.text(0.5, 0.39 / figure_height, "Row labels refer to final QUERY-PATCH Conv-IC performance.",
                ha="center", va="center", fontsize=8)
    figure.legend(handles=[Patch(facecolor=red, label="Contact"),
                           Patch(facecolor=white, edgecolor="#999999", label="No contact"),
                           Patch(facecolor=gray, label="Unannotated")],
                  loc="lower center", bbox_to_anchor=(0.5, 0.006), ncol=3, frameon=False, fontsize=8)
    # One or more export pixels per residue cell prevents large maps from losing
    # individual contacts in the scientific export; preview resolution is lower.
    dpi = max(500, int(np.ceil(max(item["truth"].shape[1] for item in examples) / panel_width)))
    figure.savefig(artifacts / "contact_additive_comparison.pdf", dpi=dpi)
    figure.savefig(output / "contact_additive_comparison.png", dpi=dpi)
    figure.savefig(output / "contact_additive_comparison_preview.png", dpi=160)
    plt.close(figure)
    return dpi


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--additive-analysis", default=DEFAULT_ANALYSIS, type=Path,
                        help="Completed select_contact_threshold.py output directory")
    parser.add_argument("--reference-dir", default=DEFAULT_REFERENCE, type=Path,
                        help="Frozen three-example figure inputs with manifest.json")
    parser.add_argument("--output-dir", default=REPO_ROOT / "figures", type=Path,
                        help="Directory for PNG figures only; reruns replace these exports")
    parser.add_argument("--artifacts-dir", default=DEFAULT_ARTIFACTS, type=Path,
                        help="Separate run-results directory for PDF, metrics, and provenance")
    args = parser.parse_args()
    analysis = args.additive_analysis.resolve(strict=True)
    reference = args.reference_dir.resolve(strict=True)
    output = args.output_dir.resolve()
    artifacts = args.artifacts_dir.resolve()
    examples, provenance = load_examples(analysis, reference)
    output.mkdir(parents=True, exist_ok=True)
    artifacts.mkdir(parents=True, exist_ok=True)
    provenance["dpi"] = render_figure(examples, output, artifacts)
    records = []
    for example in examples:
        for model in ("qp_additive", "qp_conv_ic"):
            records.append({**example["selection"], "selection_metric": "final_query_patch_f1",
                            "model": model, "threshold": example["thresholds"][model],
                            "contact_map_sha256": example["contact_map_sha256"],
                            **example["metrics"][model]})
    with (artifacts / "selected_examples.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    provenance["examples"] = records
    provenance["input_sha256"][str(Path(__file__).resolve())] = sha256(Path(__file__).resolve())
    # Record the helper actually imported, now located in training_pipeline.
    helper = Path(contact_threshold_metrics.__file__).resolve()
    provenance["input_sha256"][str(helper)] = sha256(helper)
    # Explicit paths keep the two output locations verifiable without including
    # an old provenance file or other figures in this export's own hashes.
    output_paths = (output / "contact_additive_comparison.png",
                    output / "contact_additive_comparison_preview.png",
                    artifacts / "contact_additive_comparison.pdf",
                    artifacts / "selected_examples.csv")
    provenance["output_paths"] = {path.name: str(path) for path in output_paths}
    provenance["output_sha256"] = {path.name: sha256(path) for path in output_paths}
    (artifacts / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps({"output": str(output), "artifacts": str(artifacts), "examples": records}, indent=2))


if __name__ == "__main__":
    main()
