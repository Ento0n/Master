import os
import sys
from datetime import datetime
from typing import Dict

import pandas as pd

DATA_DIR = "/nfs/scratch/pdb_dimers/"

def main():
    # print timestamp
    print(f"Starting graph creation at {datetime.now().isoformat(timespec='seconds')}")

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




    


if __name__ == "__main__":
    main()