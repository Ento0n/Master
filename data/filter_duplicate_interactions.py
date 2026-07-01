import sys
import os

import pandas as pd

DATA_DIR = "/nfs/scratch/pdb_dimers/"

def main():
    # Load the interaction dataframe
    int_df = pd.read_csv(os.path.join(DATA_DIR, "filtered_interactions_with_clusters.tsv"), sep="\t")

    # Identify duplicate interactions based on clusterpair in new_cluster_pair
    int_df["tmp_cluster_pair"] = int_df["new_cluster_pair"].apply(lambda x: ",".join(sorted(x.split(","))))

    # Filter out duplicate interactions, keeping only one representative for each cluster pair
    filtered_int_df = int_df.drop_duplicates(
        subset=["tmp_cluster_pair"],
        keep="first",
    ).copy()

    print(f"Number of interactions before filtering duplicates: {len(int_df)}")
    print(f"Number of interactions after filtering duplicates: {len(filtered_int_df)}")

    # Save the filtered interactions to a new TSV file
    filtered_int_df.drop(columns=["tmp_cluster_pair"], inplace=True)
    filtered_int_df.to_csv(os.path.join(DATA_DIR, "final_filtered_interactions.tsv"), sep="\t", index=False)

if __name__ == "__main__":
    main()