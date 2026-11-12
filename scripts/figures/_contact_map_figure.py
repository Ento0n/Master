"""Shared integrity checks and full-map rendering for the two thesis figures.

Inputs are unpadded [chain_1_length, chain_2_length] arrays with -1/0/1
encoding unknown/noncontact/contact. Both figures keep the original unknown
mask and use square cells; only their two prediction columns differ.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_REFERENCE = (REPO_ROOT / "jobs/training_pipeline/results/prostt5_additive_contact_figure"
                     / "7171489/qp_additive_normal_seed_42/figure_1sez/reference")


# =============================================================================
# Source integrity and per-map metric verification
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


# =============================================================================
# Shared three-column layout: ground truth and two checkpoint predictions
# =============================================================================

def render_figure(examples: list, output: Path, artifacts: Path, *,
                  stem: str, columns: tuple[tuple[str, str], ...]) -> int:
    """Preserve full maps; export PNGs separately from the PDF and metadata."""
    if len(columns) != 2 or len(examples) != 3:
        raise ValueError("The thesis layout requires three examples and two prediction columns")
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
    titles = ("Ground truth", *(title for _, title in columns))
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
        matrices = (truth, *(example["predictions"][model] for model, _ in columns))
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
                metrics = example["metrics"][columns[column - 1][0]]
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
    # individual contacts in the scientific export.
    dpi = max(500, int(np.ceil(max(item["truth"].shape[1] for item in examples) / panel_width)))
    figure.savefig(artifacts / f"{stem}.pdf", dpi=dpi)
    figure.savefig(output / f"{stem}.png", dpi=dpi)
    plt.close(figure)
    return dpi
