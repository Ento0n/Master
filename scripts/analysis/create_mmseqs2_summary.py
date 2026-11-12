import sys
import os

import pandas as pd

DATA_DIR = "/nfs/scratch/pdb_dimers"

def main():
    # Load the MMseqs2 hits file
    hits_df = pd.read_csv(os.path.join(DATA_DIR, "MMseqs2", "hits.tsv"), sep="\t", header=None, names=["query", "target", "pident", "nindent", "alnlen", "qcov", "tcov", "evalue", "bits"])
    
    # Remove self hits
    hits_df = hits_df[hits_df["query"] != hits_df["target"]]

    # Add partition information as sole columns
    hits_df["query_partition"] = hits_df["query"].str.split("|").str[0].astype(int)
    hits_df["target_partition"] = hits_df["target"].str.split("|").str[0].astype(int)

    # Group by query and target and calculate the average pident for each pair
    summary = (
        hits_df.groupby(["query_partition", "target_partition"])["pident"]
        .agg(
            max_identity="max",
            avg_identity="mean",
            count="count"
        )
        .reset_index()
    )

    # Save the summary to a new TSV file
    summary.to_csv(os.path.join(DATA_DIR, "MMseqs2", "partition_identity_summary.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
