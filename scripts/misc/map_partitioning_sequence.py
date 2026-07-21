import sys
import os

import pandas as pd

DATA_DIR = "/nfs/scratch/pdb_dimers"

def main():
    # Load the datasets
    seqs_df = pd.read_csv(os.path.join(DATA_DIR, "unique_sequences.tsv"), sep="\t")
    mapping_df = pd.read_csv(os.path.join(DATA_DIR, "KaHIP", "cluster_id_mapping.tsv"), sep="\t")
    partition_df = pd.read_csv(os.path.join(DATA_DIR, "KaHIP", "output.txt"), header=None, names=["partition"])

    # Add the cluster number to the partition dataframe
    partition_df["cluster_number"] = partition_df.index.map(lambda x: x + 1)

    print(partition_df.head())

    # Merge the cluster number with the mapping dataframe to get the sequence IDs
    merged_df = pd.merge(mapping_df, partition_df, left_on="id", right_on="cluster_number", how="left")

    print(merged_df.head())

    # Merge the cluster column now with the sequences dataframe to append the cluster information to the sequences
    final_df = pd.merge(seqs_df, merged_df[["partition", "cluster"]], left_on="new_cluster_id", right_on="cluster", how="left")

    print(final_df.head())

    # Save the new dataframe with the cluster information
    final_df.to_csv(os.path.join(DATA_DIR, "unique_sequences_with_partitioning.tsv"), sep="\t", index=False)

if __name__ == "__main__":
    main()