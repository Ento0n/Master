import os
import sys
from datetime import datetime

import pandas as pd

DATA_DIR = "/nfs/scratch/pdb_dimers/"

def main():
    # Load the dataframes
    seqs_df = pd.read_csv(os.path.join(DATA_DIR, "entity_sequences.tsv"), sep="\t")
    int_df = pd.read_csv(os.path.join(DATA_DIR, "filtered_interactions.tsv"), sep="\t")

    # Only keep rows in seqs_df that are present in int_df
    valid_entity_names = set()
    for _, row in int_df.iterrows():
        entity_pair = row["entity_pair"].split(",")
        valid_entity_names.update([entity.strip() for entity in entity_pair])
    seqs_df = seqs_df[seqs_df["entity_name"].isin(valid_entity_names)]

    # Group by unique Sequence and aggregate entity names, cluster IDs
    sequence2clusterids_df = (
        seqs_df.groupby("sequence")
            .agg({
                "entity_name": lambda x: list(x.unique()),
                "cluster_id": lambda x: list(x.unique()),
            })
            .reset_index()
    )

    # Is there any cluster with more than one unique sequence? If so, save the connected PDB IDs, print a Warning and the corresponding cluster ID, entity names, and sequences
    falsy_pdb_ids = []
    for i, row in sequence2clusterids_df.iterrows():
        if len(row["cluster_id"]) > 1:
            falsy_pdb_ids.extend([entity.split("_")[0] for entity in row["entity_name"]])
            print(f"Warning: Row {i} has more than one unique cluster ID. Entity names: {row['entity_name']}, Cluster IDs: {row['cluster_id']}, Sequence: {row['sequence']}", file=sys.stderr)
    
    # Remove rows with falsy cluster assignments from the interaction dataframe
    print(f"Removing interactions involving PDB IDs: {set(falsy_pdb_ids)}")
    print(f"Number of interactions before filtering: {len(int_df)}")
    int_df = int_df[~int_df["pdb_id"].isin(falsy_pdb_ids)]
    print(f"Number of interactions after filtering: {len(int_df)}")

    # Give each row a new cluster ID based on the sequence, e.g. "fix_1_100", "fix_2_100", etc.
    sequence2clusterids_df["new_cluster_id"] = sequence2clusterids_df.index.map(lambda x: f"fix_{x+1}_100")

    print(f"Check generated cluster names:\b{sequence2clusterids_df['new_cluster_id'].head()}")

    # Merge the new cluster IDs back to the original seqs_df
    seqs_df = seqs_df.merge(sequence2clusterids_df[["sequence", "new_cluster_id"]], on="sequence", how="left")

    # Add the new cluster IDs to the interaction dataframe
    new_cluster_pairs = []
    for i, row in int_df.iterrows():
        entity_pair = row["entity_pair"].split(",")
        cluster_ids = []
        for entity in entity_pair:
            cluster_id = seqs_df[seqs_df["entity_name"] == entity]["new_cluster_id"].values
            if len(cluster_id) > 0:
                cluster_ids.append(cluster_id[0])
            else:
                print(f"Warning: Entity {entity} not found in seqs_df", file=sys.stderr)
                cluster_ids.append("unknown_cluster")
        new_cluster_pairs.append(",".join(cluster_ids))

    int_df.insert(
        int_df.columns.get_loc("cluster_pair_100pct") + 1,
        "new_cluster_pair",
        new_cluster_pairs,
    )

    # Save the new interaction dataframe and the new sequence-to-cluster mapping with the new cluster pairs to a new TSV file
    int_df.to_csv(os.path.join(DATA_DIR, "filtered_interactions_with_clusters.tsv"), sep="\t", index=False)
    sequence2clusterids_df.to_csv(os.path.join(DATA_DIR, "unique_sequences.tsv"), sep="\t", index=False)

if __name__ == "__main__":
    main()