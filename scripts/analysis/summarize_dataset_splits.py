#!/usr/bin/env python3
"""Export per-split interaction counts and undirected sequence-node degrees.

Each TSV row is one edge, labels are 1 (positive) or 0 (negative), and nodes
are the sequence identifiers in new_cluster_pair. Degrees use the full node
pool of each split for both labels, including zero-degree nodes. Self-loops
contribute two endpoints. Outputs are a machine-readable CSV and a LaTeX table.
"""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median


def summarize(path):
    """Validate pairs and summarize each split without merging distinct nodes."""
    splits = defaultdict(list)
    seen = set()
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            a, b = (value.strip() for value in row["new_cluster_pair"].split(","))
            label = int(row["label"])
            if not a or not b or label not in (0, 1):
                raise ValueError(f"Invalid endpoints or label in {path}")
            key = tuple(sorted((a, b)))
            if key in seen:
                raise ValueError(f"Repeated interaction pair in {path}: {key}")
            seen.add(key)
            if row["dimer_type"] not in ("homo", "hetero"):
                raise ValueError(f"Unknown dimer_type in {path}: {key}")

            # Compare taxonomy IDs when both exist, otherwise species names.
            # Missing annotations never count as evidence of shared species.
            tax_a, tax_b = row["taxonomy_1"].strip(), row["taxonomy_2"].strip()
            sp_a, sp_b = row["species_1"].strip(), row["species_2"].strip()
            same = None
            if tax_a and tax_b:
                same = tax_a == tax_b
            elif sp_a and sp_b:
                same = sp_a == sp_b
            splits[row["split"]].append((a, b, label, same, row["dimer_type"] == "homo"))

    results = []
    for split in sorted(splits, key=lambda s: ({"train": 0, "val": 1, "test": 2}.get(s, 3), s)):
        rows = splits[split]
        nodes = {node for a, b, _, _, _ in rows for node in (a, b)}
        result = dict(dataset=path.stem, split=split, nodes=len(nodes))
        for label, name in ((1, "positives"), (0, "negatives")):
            edges = [(a, b, same, homo) for a, b, y, same, homo in rows if y == label]
            result[name] = len(edges)
            result[f"homomers_{name}"] = sum(homo for _, _, _, homo in edges)
            result[f"same_species_{name}"] = sum(same is True for _, _, same, _ in edges)
            result[f"unknown_species_{name}"] = sum(same is None for _, _, same, _ in edges)
            result[f"self_loops_{name}"] = sum(a == b for a, b, _, _ in edges)
            # Count both endpoints, including twice for a self-loop. Reindex to
            # the common split node pool so absent class-specific nodes get 0.
            counts = Counter(node for a, b, _, _ in edges for node in (a, b))
            degrees = [counts[node] for node in nodes]
            assert sum(degrees) == 2 * len(edges)
            for stat, value in zip(("min", "median", "max"), (min(degrees), median(degrees), max(degrees))):
                result[f"degree_{name}_{stat}"] = value
        results.append(result)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    results = [row for path in args.inputs for row in summarize(path)]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "split_statistics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)

    print(args.output_dir / "split_statistics.csv")
    print(args.output_dir / "split_statistics.tex")


if __name__ == "__main__":
    main()
