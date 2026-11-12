#!/usr/bin/env python3
"""Create contact_maps.png: ground truth, D-SCRIPT, and QUERY-PATCH Conv-IC.

Run with the training environment's Python; only NumPy and Matplotlib are needed:
    python scripts/figures/plot_contact_map_comparison.py

Use the same frozen 1QC9, 1SEZ, and 7QAA assembly-1 test examples as
plot_additive_contact_comparison.py. Predictions come from the final ProstT5
seed-42 checkpoints (D-SCRIPT epoch 7; QUERY-PATCH Conv-IC epoch 5), each at its
own validation F1 cutoff. No inference, training, or example selection is
repeated. Verify the reference, checkpoint identities, and exported metrics
before rendering. Only PNGs go to figures/; PDF and provenance stay with the run.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

# Resolve repository imports both from the repo root and by absolute script path.
if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.figures import _contact_map_figure
from scripts.figures._contact_map_figure import (
    DEFAULT_REFERENCE, REPO_ROOT, render_figure, sha256, verify_metrics,
)
from scripts.training_pipeline import contact_threshold_metrics
from scripts.misc.contact_reporting import normalize_contact_summary

DEFAULT_ARTIFACTS = (REPO_ROOT / "jobs/training_pipeline/results/final_comparison"
                     / "prostt5/6956426/contact_map_figure")
MODELS = (("dscript", "D-SCRIPT"), ("query_patch", "QUERY-PATCH Conv-IC"))
EXAMPLES = (("Weak", "1QC9", 1), ("Moderate", "1SEZ", 1), ("Strong", "7QAA", 1))


# =============================================================================
# Frozen examples, checkpoint provenance, and independently checked map metrics
# =============================================================================

def load_examples(reference: Path) -> tuple[list, dict]:
    """Load full [-1,0,1] maps; exclude unknown cells from every metric.

    Frozen binary maps reproduce the original threshold decisions without GPU
    inference. F1 and confusion counts are recomputed from those maps. AUPRC
    needs continuous scores, so it is taken from the hashed analysis table and
    checked against both frozen CSV and provenance records, never from binaries.
    """
    reference = reference.resolve(strict=True)
    manifest_path = reference / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    expected_files = {"selected_examples.csv", "provenance.json"} | {
        f"{pdb_id}_assembly{assembly}.npz" for _, pdb_id, assembly in EXAMPLES}
    if not expected_files.issubset(manifest["file_sha256"]):
        raise ValueError("Reference manifest must hash the selection, provenance, and all three maps")
    input_hashes = {str(manifest_path): sha256(manifest_path)}
    for name, expected_hash in manifest["file_sha256"].items():
        path = (reference / name).resolve(strict=True)
        if Path(name).is_absolute() or not path.is_relative_to(reference):
            raise ValueError(f"Reference manifest entry escapes its directory: {name}")
        input_hashes[str(path)] = sha256(path)
        if input_hashes[str(path)] != expected_hash:
            raise ValueError(f"Frozen reference hash mismatch: {name}")
    provenance = json.loads((reference / "provenance.json").read_text())
    if provenance["embedding"] != "ProstT5" or provenance["seed"] != 42:
        raise ValueError("The reference must use ProstT5 embeddings and seed 42")
    with (reference / "selected_examples.csv").open(newline="") as handle:
        selected = list(csv.DictReader(handle))
    expected_keys = {(model, pdb_id, assembly) for model, _ in MODELS
                     for _, pdb_id, assembly in EXAMPLES}
    records = {(row["model"], row["pdb_id"], int(row["assembly"])): row for row in selected}
    if len(records) != len(selected) or set(records) != expected_keys:
        raise ValueError("Reference must contain exactly the three examples for both final models")

    # Verify source summaries and actual checkpoint files. The saved per-map
    # table links each binary map to a test identity and its validation cutoff.
    tables, checkpoints = {}, {}
    for model, _ in MODELS:
        source = provenance["source_summaries"][model]
        summary_path = Path(source["path"]).resolve(strict=True)
        input_hashes[str(summary_path)] = sha256(summary_path)
        if input_hashes[str(summary_path)] != source["sha256"]:
            raise ValueError(f"Source summary hash mismatch: {summary_path}")
        summary = normalize_contact_summary(json.loads(summary_path.read_text()))
        parameters = summary["model_hyperparameters"]
        if (summary["model"] != model or summary["seed"] != 42
                or "prostt5" not in summary["embedding_manifest"].lower()
                or parameters["embedding_dim"] != 1024
                or parameters["interaction_module_type"] != "normal"):
            raise ValueError(f"Expected final ProstT5 seed-42 normal-head model: {model}")
        if model == "query_patch" and (parameters["query_evidence_mode"] != "gated"
                                       or parameters["compatibility_head"] != "dscript"):
            raise ValueError("QUERY-PATCH must be the final gated Conv-IC model")
        threshold = float(provenance["thresholds"][model])
        if (summary["selection_split"] != "val" or summary["selection_metric"] != "f1"
                or not np.isfinite(threshold) or threshold != float(summary["threshold"])):
            raise ValueError(f"Inconsistent validation F1 cutoff for {model}")
        checkpoint = Path(summary["checkpoint"]).resolve(strict=True)
        input_hashes[str(checkpoint)] = sha256(checkpoint)
        if input_hashes[str(checkpoint)] != summary["checkpoint_sha256"]:
            raise ValueError(f"Checkpoint hash mismatch: {checkpoint}")
        checkpoints[model] = {"path": str(checkpoint), "sha256": input_hashes[str(checkpoint)],
                              "model_hyperparameters": parameters}
        table_path = summary_path.parent / "test_maps.tsv"
        input_hashes[str(table_path)] = sha256(table_path)
        if input_hashes[str(table_path)] != summary["output_sha256"][table_path.name]:
            raise ValueError(f"Analysis table hash mismatch: {table_path}")
        with table_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        tables[model] = {int(row["split_index"]): row for row in rows}
        if len(tables[model]) != len(rows) or any(row["split"] != "test" for row in rows):
            raise ValueError(f"Analysis table must contain unique test split indices: {model}")

    examples = []
    for role, pdb_id, assembly in EXAMPLES:
        selection = records["query_patch", pdb_id, assembly]
        split_index = int(selection["split_index"])
        with np.load(reference / f"{pdb_id}_assembly{assembly}.npz", allow_pickle=False) as cache:
            truth = cache["labels"]
            predictions = {model: cache[model] for model, _ in MODELS}
        if (truth.ndim != 2 or not all(truth.shape) or not np.isin(truth, (-1, 0, 1)).all()):
            raise ValueError(f"Malformed ground truth: {pdb_id}")
        # Retain the full sequence axes, including unannotated residues. The
        # identical unknown mask is gray in all columns and excluded from F1/AP.
        known = truth >= 0
        metrics = {}
        truth_path = Path(tables["query_patch"][split_index]["contact_map_path"]).resolve(strict=True)
        truth_hash = sha256(truth_path)
        input_hashes[str(truth_path)] = truth_hash
        if not np.array_equal(np.load(truth_path, allow_pickle=False), truth):
            raise ValueError(f"Source ground truth differs from the frozen map: {pdb_id}")
        for model, _ in MODELS:
            row = records[model, pdb_id, assembly]
            source_row = tables[model][split_index]
            predicted = predictions[model]
            if (int(row["seed"]) != 42 or row["selection_role"] != role
                    or int(row["split_index"]) != split_index
                    or (source_row["pdb_id"], int(source_row["assembly_number"])) != (pdb_id, assembly)
                    or source_row["contact_map_sha256"] != truth_hash
                    or float(row["threshold"]) != float(provenance["thresholds"][model])
                    or truth.shape != (int(row["length_1"]), int(row["length_2"]))):
                raise ValueError(f"Inconsistent identity, shape, threshold, or ground truth: {pdb_id}/{model}")
            if (predicted.shape != truth.shape or not np.isin(predicted, (-1, 0, 1)).all()
                    or not np.array_equal(predicted < 0, ~known)):
                raise ValueError(f"Malformed binary map or unknown mask: {pdb_id}/{model}")
            if (int(row["unknown_cells"]) != int((~known).sum())
                    or not np.isclose(float(row["known_fraction"]), known.mean(), rtol=0, atol=1e-12)):
                raise ValueError(f"Reference annotation coverage mismatch: {pdb_id}/{model}")
            metrics[model] = contact_threshold_metrics.contact_metrics(truth[known], predicted[known], 0.5)
            verify_metrics(metrics[model], source_row, "at_selected_threshold_")
            metrics[model]["auprc"] = float(source_row["auprc"])
            if not 0 <= metrics[model]["auprc"] <= 1:
                raise ValueError(f"Invalid source AUPRC: {pdb_id}/{model}")
            verify_metrics(metrics[model], row)
            matching = [r for r in provenance["examples"] if r["model"] == model
                        and r["pdb_id"] == pdb_id and int(r["assembly"]) == assembly
                        and int(r["seed"]) == 42 and int(r["split_index"]) == split_index]
            if len(matching) != 1:
                raise ValueError(f"Missing or duplicate reference provenance: {pdb_id}/{model}")
            verify_metrics(metrics[model], matching[0])
        examples.append({"selection": selection, "truth": truth, "predictions": predictions,
                         "metrics": metrics, "thresholds": provenance["thresholds"],
                         "contact_map_sha256": truth_hash})
    provenance.update(input_sha256=input_hashes, checkpoints=checkpoints,
                      threshold_selection="Each model uses its frozen validation F1 threshold.",
                      ap_note="AUPRC is preserved from the hashed continuous-score analysis tables; "
                              "binary-map confusion counts and F1 are independently recomputed.")
    return examples, provenance


# =============================================================================
# Standalone export; the additive script shares the same renderer
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", default=DEFAULT_REFERENCE, type=Path,
                        help="Frozen three-example inputs shared with the additive figure")
    parser.add_argument("--output-dir", default=REPO_ROOT / "figures", type=Path,
                        help="Directory for PNG figures only; reruns replace these exports")
    parser.add_argument("--artifacts-dir", default=DEFAULT_ARTIFACTS, type=Path,
                        help="Separate run-results directory for PDF, metrics, and provenance")
    args = parser.parse_args()
    examples, provenance = load_examples(args.reference_dir)
    output, artifacts = args.output_dir.resolve(), args.artifacts_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    artifacts.mkdir(parents=True, exist_ok=True)
    provenance["dpi"] = render_figure(examples, output, artifacts, stem="contact_maps", columns=MODELS)
    records = [{**example["selection"], "model": model,
                "selection_metric": "final_query_patch_f1", "threshold": example["thresholds"][model],
                "contact_map_sha256": example["contact_map_sha256"], **example["metrics"][model]}
               for example in examples for model, _ in MODELS]
    with (artifacts / "selected_examples.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    provenance["examples"] = records
    for module_path in (Path(__file__), Path(_contact_map_figure.__file__),
                        Path(contact_threshold_metrics.__file__)):
        path = module_path.resolve()
        provenance["input_sha256"][str(path)] = sha256(path)
    output_paths = (output / "contact_maps.png",
                    artifacts / "contact_maps.pdf", artifacts / "selected_examples.csv")
    provenance["output_paths"] = {path.name: str(path) for path in output_paths}
    provenance["output_sha256"] = {path.name: sha256(path) for path in output_paths}
    (artifacts / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps({"output": str(output), "artifacts": str(artifacts), "examples": records}, indent=2))


if __name__ == "__main__":
    main()
