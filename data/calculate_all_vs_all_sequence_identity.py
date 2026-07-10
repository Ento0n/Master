from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import edlib
import numpy as np
import pandas as pd


DATA_DIR = "/nfs/scratch/pdb_dimers"
IDENTITY_DTYPE = np.float32
LONG_FORMAT_FILE = "all_vs_all_sequence_identity_long_format.tsv"
CLUSTER_ORDER_FILE = "all_vs_all_sequence_identity_cluster_order.tsv"

_CIGAR_RE = re.compile(r"(\d+)([=XIDM])")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate all-vs-all sequence identity with compact triangular storage."
    )
    parser.add_argument(
        "--local",
        action="store_true",
        help="Calculate local/HW CD-HIT-style sequence identity.",
    )
    parser.add_argument(
        "--global",
        dest="calculate_global",
        action="store_true",
        help="Calculate global/NW sequence identity.",
    )
    parser.add_argument(
        "--data-dir",
        default=DATA_DIR,
        help=f"Pipeline data directory. Default: {DATA_DIR}",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory. Default: <data-dir>/own_all_vs_all",
    )
    parser.add_argument(
        "--no-long-format",
        action="store_true",
        help="Skip the streamed long-format TSV output.",
    )
    parser.add_argument(
        "--write-square-tsv",
        action="store_true",
        help="Also write old-style dense square TSV matrices from the triangular arrays.",
    )
    args = parser.parse_args()

    if not args.local and not args.calculate_global:
        parser.error("At least one of --local or --global must be set.")

    return args


def triangular_size(n_items: int) -> int:
    return n_items * (n_items - 1) // 2


def triangular_index(i: int, j: int, n_items: int) -> int:
    if i == j:
        raise ValueError("diagonal entries are not stored in the triangular matrix")
    if i > j:
        i, j = j, i
    return i * n_items - i * (i + 1) // 2 + (j - i - 1)


def count_matches_from_cigar(cigar: str) -> int | None:
    if not cigar:
        return None

    matches = 0
    parsed_chars = 0
    for match in _CIGAR_RE.finditer(cigar):
        parsed_chars += len(match.group(0))
        length = int(match.group(1))
        operation = match.group(2)
        if operation == "=":
            matches += length
        elif operation == "M":
            # "M" does not distinguish matches from mismatches, so fall back.
            return None

    if parsed_chars != len(cigar):
        return None

    return matches


def calc_local_identity(seq1: str, seq2: str) -> float:
    """
    CD-HIT-style identity via edlib:
      - mode="HW": alignment must cover the whole query (the shorter sequence),
        but can start/end anywhere on the target (the longer sequence).
      - Count identical aligned columns (excluding gaps) and divide by len(shorter).
    """
    if len(seq1) <= len(seq2):
        query, target = seq1, seq2
    else:
        query, target = seq2, seq1

    # The original local metric is matches / len(shorter), excluding gap
    # columns. editDistance alone cannot recover that value when the local
    # alignment contains target insertions, because insertions increase edit
    # distance but do not consume query residues.
    res = edlib.align(query, target, mode="HW", task="path")
    matches = count_matches_from_cigar(res.get("cigar"))

    if matches is None:
        nice = edlib.getNiceAlignment(res, query, target)
        matches = sum(
            1
            for qc, tc in zip(nice["query_aligned"], nice["target_aligned"])
            if qc != "-" and tc != "-" and qc == tc
        )

    return matches / len(query)


def calc_global_identity(seq1: str, seq2: str) -> float:
    result = edlib.align(seq1, seq2, mode="NW", task="distance")
    edits = result["editDistance"]
    return 1 - (edits / max(len(seq1), len(seq2)))


def selected_metric_names(calculate_local: bool, calculate_global: bool) -> list[str]:
    metric_names = []
    if calculate_local:
        metric_names.append("local")
    if calculate_global:
        metric_names.append("global")
    return metric_names


def identity_array_path(output_dir: Path, metric_name: str) -> Path:
    return output_dir / f"{metric_name}_all_vs_all_sequence_identity.triu.npy"


def square_tsv_path(output_dir: Path, metric_name: str) -> Path:
    return output_dir / f"{metric_name}_all_vs_all_sequence_identity.tsv"


def open_identity_arrays(output_dir: Path, metric_names: list[str], n_items: int):
    shape = (triangular_size(n_items),)
    arrays = {}
    for metric_name in metric_names:
        arrays[metric_name] = np.lib.format.open_memmap(
            identity_array_path(output_dir, metric_name),
            mode="w+",
            dtype=IDENTITY_DTYPE,
            shape=shape,
        )
    return arrays


def format_identity(value: float) -> str:
    return f"{float(value):.8g}"


def load_sequences(data_dir: Path) -> pd.DataFrame:
    seqs_df = pd.read_csv(
        data_dir / "unique_sequences.tsv",
        sep="\t",
        usecols=["new_cluster_id", "sequence"],
    )
    int_df = pd.read_csv(
        data_dir / "final_filtered_interactions.tsv",
        sep="\t",
        usecols=["new_cluster_pair"],
    )
    used_clusters = set(int_df["new_cluster_pair"].str.split(",").explode().unique())
    seqs_df = seqs_df[seqs_df["new_cluster_id"].isin(used_clusters)].reset_index(drop=True)

    missing_clusters = used_clusters.difference(seqs_df["new_cluster_id"])
    if missing_clusters:
        preview = ", ".join(sorted(missing_clusters)[:10])
        raise ValueError(
            f"{len(missing_clusters)} clusters from final_filtered_interactions.tsv "
            f"are missing in unique_sequences.tsv. Examples: {preview}"
        )

    return seqs_df


def write_cluster_order(output_dir: Path, clusters: list[str]) -> None:
    with (output_dir / CLUSTER_ORDER_FILE).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["index", "new_cluster_id"])
        for index, cluster_id in enumerate(clusters):
            writer.writerow([index, cluster_id])


def open_long_format_writer(output_dir: Path, metric_names: list[str]):
    handle = (output_dir / LONG_FORMAT_FILE).open("w", newline="")
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["cluster1", "cluster2", *(f"{name}_identity" for name in metric_names)])
    return handle, writer


def write_square_tsv(
    output_path: Path,
    values: np.memmap,
    clusters: list[str],
    diagonal_value: float = 1.0,
) -> None:
    n_items = len(clusters)
    with output_path.open("w", newline="") as handle:
        handle.write("\t")
        handle.write("\t".join(clusters))
        handle.write("\n")

        for i, cluster_id in enumerate(clusters):
            row_values = []
            for j in range(n_items):
                if i == j:
                    value = diagonal_value
                else:
                    value = values[triangular_index(i, j, n_items)]
                row_values.append(format_identity(value))
            handle.write(cluster_id)
            handle.write("\t")
            handle.write("\t".join(row_values))
            handle.write("\n")


def calculate_identities(
    seqs: list[str],
    clusters: list[str],
    arrays: dict[str, np.memmap],
    long_writer: csv.writer | None,
) -> None:
    n_items = len(seqs)
    total_pairs = triangular_size(n_items)
    completed_pairs = 0

    for i in range(n_items - 1):
        seq1 = seqs[i]
        cluster1 = clusters[i]

        for j in range(i + 1, n_items):
            index = triangular_index(i, j, n_items)
            row = [cluster1, clusters[j]] if long_writer is not None else None

            if "local" in arrays:
                local_identity = IDENTITY_DTYPE(calc_local_identity(seq1, seqs[j]))
                arrays["local"][index] = local_identity
                if row is not None:
                    row.append(format_identity(local_identity))

            if "global" in arrays:
                global_identity = IDENTITY_DTYPE(calc_global_identity(seq1, seqs[j]))
                arrays["global"][index] = global_identity
                if row is not None:
                    row.append(format_identity(global_identity))

            if long_writer is not None and row is not None:
                long_writer.writerow(row)

            completed_pairs += 1

        if i % 100 == 0 or i == n_items - 2:
            percent = completed_pairs / total_pairs * 100 if total_pairs else 100
            print(
                f"Completed row {i + 1}/{n_items - 1}: "
                f"{completed_pairs}/{total_pairs} pairs ({percent:.2f}%)",
                flush=True,
            )


def main():
    args = parse_args()
    data_dir = Path(args.data_dir)
    output_dir = Path(args.output_dir) if args.output_dir else data_dir / "own_all_vs_all"
    output_dir.mkdir(parents=True, exist_ok=True)

    metric_names = selected_metric_names(args.local, args.calculate_global)
    seqs_df = load_sequences(data_dir)
    seqs = seqs_df["sequence"].tolist()
    clusters = seqs_df["new_cluster_id"].tolist()

    print(
        f"Calculating {', '.join(metric_names)} all-vs-all sequence identity "
        f"for {len(seqs)} sequences ({triangular_size(len(seqs))} unique pairs)",
        flush=True,
    )

    write_cluster_order(output_dir, clusters)
    arrays = open_identity_arrays(output_dir, metric_names, len(seqs))

    long_handle = None
    long_writer = None
    if not args.no_long_format:
        long_handle, long_writer = open_long_format_writer(output_dir, metric_names)

    try:
        calculate_identities(seqs, clusters, arrays, long_writer)
    finally:
        if long_handle is not None:
            long_handle.close()
        for values in arrays.values():
            values.flush()

    if args.write_square_tsv:
        for metric_name, values in arrays.items():
            output_path = square_tsv_path(output_dir, metric_name)
            print(f"Writing dense square TSV: {output_path}", flush=True)
            write_square_tsv(output_path, values, clusters)

    for metric_name in metric_names:
        print(f"Wrote compact triangular array: {identity_array_path(output_dir, metric_name)}", flush=True)
    print(f"Wrote cluster order: {output_dir / CLUSTER_ORDER_FILE}", flush=True)
    if not args.no_long_format:
        print(f"Wrote streamed long-format TSV: {output_dir / LONG_FORMAT_FILE}", flush=True)


if __name__ == "__main__":
    main()
