"""Translate saved contact-analysis summaries to the current metric names.

Legacy names are handled only here so reporting and plotting can read existing
results. This adapter copies the relevant dictionaries without changing saved
files, recalculating metrics, or checking result integrity.
"""


def normalize_contact_summary(summary: dict) -> dict:
    """Flatten archived per-map metrics and rename their F1/AP fields.

    Current-format records retain their values. Unsupported selection objectives
    remain unchanged so callers can reject them rather than relabel their scores.
    """
    objective_names = {"macro_f1": "f1"}
    result = dict(summary)
    result["selection_metric"] = objective_names.get(summary["selection_metric"], summary["selection_metric"])
    if "search" in summary:
        search = dict(summary["search"])
        search["objective"] = objective_names.get(search["objective"], search["objective"])
        result["search"] = search

    for split in ("validation", "test"):
        if split not in summary:
            continue
        metrics = dict(summary[split])
        if "contact_auprc_macro" in metrics:
            metrics["contact_auprc"] = metrics.pop("contact_auprc_macro")
        metrics.pop("contact_auprc_micro", None)
        for decision in ("at_0.5", "at_selected_threshold"):
            record = dict(metrics[decision])
            # Archived records store the equal-map scores in a nested section.
            # Keep their cohort counts and discard the obsolete pooled scores.
            record.update(record.pop("macro", {}))
            record.pop("micro", None)
            metrics[decision] = record
        result[split] = metrics
    return result
