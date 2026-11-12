import os
from collections import defaultdict


DATA_DIR = "/nfs/scratch/pdb_dimers"
PARTITIONS_FILE = os.path.join(DATA_DIR, "KaHIP", "seq_ident_partitions_strong_output.txt")
SEQUENCE_IDENTITY_GRAPH = os.path.join(DATA_DIR, "MMseqs2", "sequence_identity.graph")
OUTPUT_FILE = os.path.join(DATA_DIR, "KaHIP", "partition_sequence_identity_summary_strong.tsv")


def load_partitions():
    with open(PARTITIONS_FILE) as f:
        return [int(line.strip()) for line in f if line.strip()]


def summarize_graph(partitions):
    summary = defaultdict(lambda: {"sum": 0.0, "max": 0.0, "count": 0})

    with open(SEQUENCE_IDENTITY_GRAPH) as f:
        n_nodes = int(f.readline().split()[0])
        if n_nodes != len(partitions):
            raise ValueError(
                f"Graph has {n_nodes} nodes, but partition file has {len(partitions)} rows"
            )

        for node_id, line in enumerate(f, start=1):
            source_partition = partitions[node_id - 1]
            fields = line.split()

            for neighbor, weight in zip(fields[0::2], fields[1::2]):
                neighbor_id = int(neighbor)
                if node_id >= neighbor_id:
                    continue

                target_partition = partitions[neighbor_id - 1]
                if source_partition == target_partition:
                    continue

                partition_pair = tuple(sorted((source_partition, target_partition)))
                identity = int(weight) / 100

                summary[partition_pair]["sum"] += identity
                summary[partition_pair]["max"] = max(summary[partition_pair]["max"], identity)
                summary[partition_pair]["count"] += 1

    return summary


def write_summary(summary):
    with open(OUTPUT_FILE, "w") as f:
        f.write("partition_a\tpartition_b\tmax_identity\tavg_identity\tcount\n")

        for (partition_a, partition_b), values in sorted(summary.items()):
            avg_identity = values["sum"] / values["count"]
            f.write(
                f"{partition_a}\t{partition_b}\t{values['max']}\t"
                f"{avg_identity}\t{values['count']}\n"
            )


def main():
    partitions = load_partitions()
    summary = summarize_graph(partitions)
    write_summary(summary)
    print(f"Wrote summary to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
