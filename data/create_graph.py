import os
from datetime import datetime
from typing import Dict
import argparse
from itertools import combinations
from collections import defaultdict

import pandas as pd
import networkx as nx

DATA_DIR = "/nfs/scratch/pdb_dimers/"

def create_metis_interaction_graph():
    # Load interacion dataframe
    interaction_df = pd.read_csv(os.path.join(DATA_DIR, "final_filtered_interactions.tsv"), sep="\t")

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
    with open(os.path.join(DATA_DIR, "KaHIP", "cluster_id_mapping.tsv"), "w") as f:
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
    with open(os.path.join(DATA_DIR, "KaHIP", "cluster.graph"), "w") as f:
        # Write header, no of nodes and no of edges
        num_edges = sum(len(neighbors) for neighbors in cluster_connections.values()) // 2
        f.write(f"{len(cluster_connections)} {num_edges}\n")
        for neighbors in cluster_connections.values():
            for neighbor in neighbors:
                f.write(f"{neighbor} ")
            f.write("\n")
    
    print(f"METIS graph saved to: {os.path.join(DATA_DIR, 'KaHIP', 'cluster.graph')}")

def interaction_weight(pair1, pair2, seq_id_df):
    a, b = pair1.split(",")
    c, d = pair2.split(",")

    return max(
        seq_id_df.loc[a, c],
        seq_id_df.loc[a, d],
        seq_id_df.loc[b, c],
        seq_id_df.loc[b, d]
    )

def create_metis_sequence_identity_graph():
    int_df = pd.read_csv(os.path.join(DATA_DIR, "final_filtered_interactions.tsv"), sep="\t")
    seq_id_df = pd.read_csv(os.path.join(DATA_DIR, "own_all_vs_all", "local_all_vs_all_sequence_identity.tsv"), sep="\t", index_col=0)
    cluster_pairs = int_df["new_cluster_pair"]

    node_id_mapping = {pair: i + 1 for i, pair in enumerate(cluster_pairs)}

    # extract sequence identity for each cluster pair
    adjacency_list = defaultdict(list)
    for pair1, pair2 in combinations(cluster_pairs, 2):
        weight = interaction_weight(pair1, pair2, seq_id_df)
        int_weight = int(weight * 100)  # Convert to integer weight for METIS

        i = node_id_mapping[pair1]
        j = node_id_mapping[pair2]

        adjacency_list[i].append((j, int_weight))
        adjacency_list[j].append((i, int_weight))

    n_nodes = len(cluster_pairs)
    n_edges = sum(len(neighbors) for neighbors in adjacency_list.values()) // 2

    with open(os.path.join(DATA_DIR, "MMseqs2", "sequence_identity.graph"), "w") as f:
        f.write(f"{n_nodes} {n_edges} 001\n")  # 001 indicates weighted graph
        for i in range(1, n_nodes + 1):
            neighbors = adjacency_list[i]
            for j, weight in neighbors:
                f.write(f"{j} {weight} ")
            f.write("\n")

def create_nx_graph_from_metis():
    G = nx.Graph()

    with open(os.path.join(DATA_DIR, "KaHIP", "cluster.graph"), "r") as f:
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

    nx.write_gexf(G, os.path.join(DATA_DIR, "KaHIP", "cluster.graph.gexf"))
    print(f"NetworkX graph saved to: {os.path.join(DATA_DIR, 'KaHIP', 'cluster.graph.gexf')}")

def main():
    parser = argparse.ArgumentParser(description="Create a graph from the final filtered interactions and save it in METIS format.")
    parser.add_argument("--type", choices=["interaction", "sequence_identity"], help="Type of graph to create")
    args = parser.parse_args()

    # print timestamp
    print(f"Starting graph creation at {datetime.now().isoformat(timespec='seconds')}")

    if args.type == "interaction":
        create_metis_interaction_graph()
    else:
        create_metis_sequence_identity_graph()

    print(f"Graph creation completed at {datetime.now().isoformat(timespec='seconds')}")

    if args.type == "interaction":
        create_nx_graph_from_metis()
        print(f"NetworkX graph created from METIS file at {datetime.now().isoformat(timespec='seconds')}")

if __name__ == "__main__":
    main()