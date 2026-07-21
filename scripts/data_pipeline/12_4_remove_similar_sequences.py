#!/usr/bin/env python3
"""Filter validation and test interactions using CD-HIT keep-lists.

The CD-HIT outputs are FASTA-like files with headers such as:
    >0_test|fix_14530_100

For rows in the interaction TSV with split == "test", both cluster ids in
new_cluster_pair must be present in train_test.out.  For split == "val", both
cluster ids must be present in train_val.out.  Other splits are copied through.
"""

from __future__ import annotations

import argparse
from collections.abc import Iterable
from pathlib import Path

import pandas as pd


DEFAULT_INPUT = Path("/nfs/scratch/pdb_dimers/dataset_iterations/06_removed_suspicious_contacts.tsv")
DEFAULT_OUTPUT = Path("/nfs/scratch/pdb_dimers/dataset_iterations/07_cdhit.tsv")
DEFAULT_VAL_CLUSTERS = Path("/nfs/scratch/pdb_dimers/CD-HIT/train_val.out")
DEFAULT_TEST_CLUSTERS = Path("/nfs/scratch/pdb_dimers/CD-HIT/train_test.out")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Remove val/test interaction rows whose new_cluster_pair contains "
            "cluster ids not kept by the matching CD-HIT output."
        )
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Interaction TSV to filter.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Filtered TSV to write.")
    parser.add_argument(
        "--val-clusters",
        type=Path,
        default=DEFAULT_VAL_CLUSTERS,
        help="CD-HIT output containing allowed validation cluster ids.",
    )
    parser.add_argument(
        "--test-clusters",
        type=Path,
        default=DEFAULT_TEST_CLUSTERS,
        help="CD-HIT output containing allowed test cluster ids.",
    )
    return parser.parse_args()


def load_allowed_clusters(path: Path) -> set[str]:
    """Read cluster ids from FASTA headers in a CD-HIT output file."""
    clusters: set[str] = set()
    bad_headers: list[str] = []

    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or not line.startswith(">"):
                continue

            header = line[1:].split(None, 1)[0]
            if "|" not in header:
                bad_headers.append(line)
                continue

            cluster_id = header.rsplit("|", 1)[1].strip()
            if cluster_id:
                clusters.add(cluster_id)
            else:
                bad_headers.append(line)

    if bad_headers:
        examples = ", ".join(bad_headers[:5])
        raise ValueError(f"{path} contains FASTA headers without cluster ids after '|': {examples}")
    if not clusters:
        raise ValueError(f"No FASTA headers with cluster ids were found in {path}")

    return clusters


def split_cluster_pair(value: str, row_index: object) -> tuple[str, str]:
    """Split one new_cluster_pair value and validate that it has two endpoints."""
    parts = [part.strip() for part in str(value).split(",")]
    if len(parts) != 2 or not parts[0] or not parts[1]:
        raise ValueError(
            f"new_cluster_pair must contain two comma-separated cluster ids; "
            f"bad value at row {row_index}: {value!r}"
        )
    return parts[0], parts[1]


def row_allowed(row: pd.Series, allowed_by_split: dict[str, set[str]]) -> bool:
    """Return whether an interaction row should be kept."""
    split = str(row["split"]).strip()
    allowed_clusters = allowed_by_split.get(split)
    if allowed_clusters is None:
        return True

    cluster_a, cluster_b = split_cluster_pair(row["new_cluster_pair"], row.name)
    return cluster_a in allowed_clusters and cluster_b in allowed_clusters


def print_split_counts(label: str, splits: Iterable[str], df: pd.DataFrame) -> None:
    counts = df["split"].value_counts().to_dict()
    print(label)
    for split in splits:
        print(f"  {split}: {counts.get(split, 0)}")


def read_interactions(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

    required_columns = {"new_cluster_pair", "split"}
    missing_columns = sorted(required_columns - set(df.columns))
    if missing_columns:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing_columns)}")

    # Preserve a leading blank index column from older pandas outputs when round-tripping.
    for column in df.columns:
        if column.startswith("Unnamed:"):
            df = df.rename(columns={column: ""})
            break

    return df


def main() -> None:
    args = parse_args()

    allowed_by_split = {
        "val": load_allowed_clusters(args.val_clusters),
        "test": load_allowed_clusters(args.test_clusters),
    }
    print(f"Allowed val clusters: {len(allowed_by_split['val'])}")
    print(f"Allowed test clusters: {len(allowed_by_split['test'])}")

    interactions = read_interactions(args.input)
    split_order = list(interactions["split"].drop_duplicates())
    print_split_counts("Before filtering:", split_order, interactions)

    keep_mask = interactions.apply(row_allowed, axis=1, allowed_by_split=allowed_by_split)
    filtered = interactions.loc[keep_mask].copy()

    print_split_counts("After filtering:", split_order, filtered)
    print(f"Removed rows: {len(interactions) - len(filtered)}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    filtered.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote filtered interactions to {args.output}")


if __name__ == "__main__":
    main()
