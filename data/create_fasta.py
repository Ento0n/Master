import os
import sys
from datetime import datetime

import pandas as pd

DATA_DIR = "/nfs/scratch/pdb_dimers/"

def main():
    # print timestamp
    print(f"Starting FASTA file creation at {datetime.now().isoformat(timespec='seconds')}")

    # Load the sequence file
    seqs_df = pd.read_csv(os.path.join(DATA_DIR, "entity_sequences.tsv"), sep="\t")

    # Group by unique cluster IDs and aggregate entity names, sequences should be identical therefore no list
    fasta_df = (
        seqs_df.groupby("sequence")
            .agg({
                "entity_name": lambda x: list(x.unique()),
                "cluster_id": lambda x: list(x.unique()),
            })
            .reset_index()
    )

    # Is there any cluster with more than one unique sequence? If so, print a Warning and the corresponding cluster ID, entity names, and sequences
    for i, row in fasta_df.iterrows():
        if len(row["cluster_id"]) > 1:
            print(f"Warning: Row {i} has more than one unique cluster ID. Entity names: {row['entity_name']}, Cluster IDs: {row['cluster_id']}, Sequence: {row['sequence']}", file=sys.stderr)

    # Write the FASTA file
    with open(os.path.join(DATA_DIR, "output.fasta"), "w") as fasta_file:
        for i, row in fasta_df.iterrows():
            cluster_ids = row["cluster_id"]
            entity_names = row["entity_name"]
            sequence = row["sequence"]

            # Create a FASTA header with the cluster ID and entity names
            header = f">{i};{' '.join(cluster_ids)};{' '.join(entity_names)}"
            fasta_file.write(f"{header}\n")
            fasta_file.write(f"{sequence}\n")

    print(f"FASTA file created at {os.path.join(DATA_DIR, 'output.fasta')}")
    print(f"Finished FASTA file creation at {datetime.now().isoformat(timespec='seconds')}")


if __name__ == "__main__":
    main()