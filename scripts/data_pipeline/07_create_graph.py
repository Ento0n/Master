import os
from datetime import datetime
from typing import Dict
import argparse

import numpy as np
import pandas as pd

DATA_DIR = "/nfs/scratch/pdb_dimers"
CLUSTER_ORDER_FILE = "all_vs_all_sequence_identity_cluster_order.tsv"
INTERACTION_FILE = os.path.join("dataset_iterations", "04_removed_duplicates.tsv")


def triangular_size(n_items: int) -> int:
    return n_items * (n_items - 1) // 2


def triangular_index(i: int, j: int, n_items: int) -> int:
    if i == j:
        raise ValueError("diagonal entries are not stored in the triangular matrix")
    if i > j:
        i, j = j, i
    return i * n_items - i * (i + 1) // 2 + (j - i - 1)


class TriangularSequenceIdentity:
    def __init__(self, values_path: str, cluster_order_path: str, diagonal_value: float = 1.0):
        cluster_order = pd.read_csv(cluster_order_path, sep="\t")
        required_columns = {"index", "new_cluster_id"}
        missing_columns = required_columns.difference(cluster_order.columns)
        if missing_columns:
            raise ValueError(
                f"{cluster_order_path} is missing required columns: "
                f"{', '.join(sorted(missing_columns))}"
            )

        cluster_order = cluster_order.sort_values("index").reset_index(drop=True)
        expected_index = pd.RangeIndex(len(cluster_order))
        if not cluster_order["index"].equals(expected_index.to_series()):
            raise ValueError(
                f"{cluster_order_path} must contain contiguous zero-based indices."
            )

        values = np.load(values_path, mmap_mode="r")
        if values.ndim != 1:
            raise ValueError(f"{values_path} must be a 1D compact triangular array.")

        n_items = len(cluster_order)
        expected_values = triangular_size(n_items)
        if len(values) != expected_values:
            raise ValueError(
                f"{values_path} contains {len(values)} values, but {n_items} clusters "
                f"require {expected_values} upper-triangular values."
            )

        self.values = values
        self.n_items = n_items
        self.diagonal_value = diagonal_value
        self.cluster_to_index = dict(
            zip(cluster_order["new_cluster_id"], cluster_order["index"])
        )

    def index_for_cluster(self, cluster_id: str) -> int:
        try:
            return self.cluster_to_index[cluster_id]
        except KeyError as exc:
            raise KeyError(
                f"Cluster {cluster_id!r} is missing from the sequence identity cluster order."
            ) from exc

    def identity_by_index(self, i: int, j: int) -> float:
        if i == j:
            return self.diagonal_value
        return float(self.values[triangular_index(i, j, self.n_items)])


def parse_cluster_pair_indices(
    cluster_pairs: pd.Series,
    seq_identity: TriangularSequenceIdentity,
) -> list[tuple[int, int]]:
    endpoint_indices = []
    missing_clusters = set()

    for pair in cluster_pairs:
        cluster1, cluster2 = pair.split(",")
        try:
            endpoint_indices.append(
                (
                    seq_identity.index_for_cluster(cluster1),
                    seq_identity.index_for_cluster(cluster2),
                )
            )
        except KeyError:
            missing_clusters.update(
                cluster
                for cluster in (cluster1, cluster2)
                if cluster not in seq_identity.cluster_to_index
            )

    if missing_clusters:
        preview = ", ".join(sorted(missing_clusters)[:10])
        raise ValueError(
            f"{len(missing_clusters)} clusters from {INTERACTION_FILE} "
            f"are missing from the sequence identity cluster order. Examples: {preview}"
        )

    return endpoint_indices


def create_metis_interaction_graph(data_dir: str):
    # Load interacion dataframe
    interaction_df = pd.read_csv(os.path.join(data_dir, INTERACTION_FILE), sep="\t")
    graph_dir = os.path.join(data_dir, "graphs")
    os.makedirs(graph_dir, exist_ok=True)

    # Create a integer cluster ID mapping
    unique_clusters = set()
    for pair in interaction_df["new_cluster_pair"]:
        cluster1, cluster2 = pair.split(",")
        if cluster1 == cluster2:
            continue  # Skip self-connections
        unique_clusters.add(cluster1)
        unique_clusters.add(cluster2)
    cluster_id_mapping = {cluster: idx+1 for idx, cluster in enumerate(sorted(unique_clusters))}
    print(f"Total unique clusters: {len(cluster_id_mapping)}")
    print(f"Test: {list(cluster_id_mapping.items())[:5]}")

    # Save the cluster ID mapping to a file
    with open(os.path.join(data_dir, "KaHIP", "cluster_id_mapping.tsv"), "w") as f:
        f.write("cluster\tid\n")
        for cluster, idx in cluster_id_mapping.items():
            f.write(f"{cluster}\t{idx}\n")

    # Extract which cluster is connected to which other clusters
    cluster_connections: Dict[int, set] = {}
    for _, row in interaction_df.iterrows():
        cluster1, cluster2 = row["new_cluster_pair"].split(",")
        if cluster1 == cluster2:
            continue  # Skip self-connections
        cluster_connections.setdefault(cluster_id_mapping[cluster1], set()).add(cluster_id_mapping[cluster2])
        cluster_connections.setdefault(cluster_id_mapping[cluster2], set()).add(cluster_id_mapping[cluster1])

    # Sort cluster connections by key
    cluster_connections = dict(sorted(cluster_connections.items()))
    
    print(f"Total clusters: {len(cluster_connections)}")
    print(f"Test: {list(cluster_connections.items())[:5]}")

    # Save the graph in METIS format
    with open(os.path.join(graph_dir, "cluster.graph"), "w") as f:
        # Write header, no of nodes and no of edges
        num_edges = sum(len(neighbors) for neighbors in cluster_connections.values()) // 2
        f.write(f"{len(cluster_connections)} {num_edges}\n")
        for neighbors in cluster_connections.values():
            for neighbor in neighbors:
                f.write(f"{neighbor} ")
            f.write("\n")
    
    print(f"METIS graph saved to: {os.path.join(graph_dir, 'cluster.graph')}")


def interaction_weight(
    pair1_indices: tuple[int, int],
    pair2_indices: tuple[int, int],
    seq_identity: TriangularSequenceIdentity,
) -> float:
    a, b = pair1_indices
    c, d = pair2_indices

    return max(
        seq_identity.identity_by_index(a, c),
        seq_identity.identity_by_index(a, d),
        seq_identity.identity_by_index(b, c),
        seq_identity.identity_by_index(b, d),
    )


def create_metis_sequence_identity_graph(data_dir: str, identity_metric: str):
    int_df = pd.read_csv(os.path.join(data_dir, INTERACTION_FILE), sep="\t")
    identity_dir = os.path.join(data_dir, "own_all_vs_all")
    seq_identity = TriangularSequenceIdentity(
        values_path=os.path.join(
            identity_dir,
            f"{identity_metric}_all_vs_all_sequence_identity.triu.npy",
        ),
        cluster_order_path=os.path.join(identity_dir, CLUSTER_ORDER_FILE),
        diagonal_value=1.0,
    )
    cluster_pairs = int_df["new_cluster_pair"]
    pair_indices = parse_cluster_pair_indices(cluster_pairs, seq_identity)

    n_nodes = len(pair_indices)
    n_edges = n_nodes * (n_nodes - 1) // 2
    output_path = os.path.join(data_dir, "graphs", "sequence_identity.graph")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    print(
        f"Writing complete sequence-identity graph with {n_nodes} nodes "
        f"and {n_edges} edges from {identity_metric} triangular identities.",
        flush=True,
    )

    with open(output_path, "w") as f:
        f.write(f"{n_nodes} {n_edges} 001\n")  # 001 indicates weighted graph

        for i, pair1_indices in enumerate(pair_indices):
            row_parts = []
            for j in range(i):
                pair2_indices = pair_indices[j]
                weight = interaction_weight(pair1_indices, pair2_indices, seq_identity)
                int_weight = int(weight * 100)  # Convert to integer weight for METIS
                row_parts.append(f"{j + 1} {int_weight} ")

            for j in range(i + 1, n_nodes):
                pair2_indices = pair_indices[j]
                weight = interaction_weight(pair1_indices, pair2_indices, seq_identity)
                int_weight = int(weight * 100)  # Convert to integer weight for METIS
                row_parts.append(f"{j + 1} {int_weight} ")

            f.write("".join(row_parts))
            f.write("\n")

            if (i + 1) % 100 == 0 or i == n_nodes - 1:
                print(f"Wrote adjacency row {i + 1}/{n_nodes}", flush=True)

    print(f"METIS graph saved to: {output_path}")


def create_nx_graph_from_metis(data_dir: str):
    import networkx as nx

    G = nx.Graph()

    graph_dir = os.path.join(data_dir, "graphs")
    with open(os.path.join(graph_dir, "cluster.graph"), "r") as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith("%")]

    header = lines[0].split()
    n_nodes = int(header[0])

    for i in range(1, n_nodes + 1):
        G.add_node(i)

    for i, line in enumerate(lines[1:], start=1):
        neighbors = list(map(int, line.split()))

        for j in neighbors:
            # METIS uses 1-based node IDs
            if i < j:  # avoid adding undirected edges twice
                G.add_edge(i, j)

    nx.write_gexf(G, os.path.join(graph_dir, "cluster.graph.gexf"))
    print(f"NetworkX graph saved to: {os.path.join(graph_dir, 'cluster.graph.gexf')}")

def main():
    parser = argparse.ArgumentParser(description="Create a graph from the final filtered interactions and save it in METIS format.")
    parser.add_argument("--type", choices=["interaction", "sequence_identity"], help="Type of graph to create")
    parser.add_argument(
        "--data-dir",
        default=DATA_DIR,
        help=f"Pipeline data directory. Default: {DATA_DIR}",
    )
    parser.add_argument(
        "--identity-metric",
        choices=["local", "global"],
        default="local",
        help="Triangular identity metric to use for the sequence-identity graph.",
    )
    args = parser.parse_args()

    # print timestamp
    print(f"Starting graph creation at {datetime.now().isoformat(timespec='seconds')}")

    if args.type == "interaction":
        create_metis_interaction_graph(args.data_dir)
    else:
        create_metis_sequence_identity_graph(args.data_dir, args.identity_metric)

    print(f"Graph creation completed at {datetime.now().isoformat(timespec='seconds')}")

    if args.type == "interaction":
        create_nx_graph_from_metis(args.data_dir)
        print(f"NetworkX graph created from METIS file at {datetime.now().isoformat(timespec='seconds')}")

if __name__ == "__main__":
    main()
