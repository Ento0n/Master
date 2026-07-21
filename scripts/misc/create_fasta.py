import os
import sys
from datetime import datetime

import pandas as pd

DATA_DIR = "/nfs/scratch/pdb_dimers/"

def main():
    # print timestamp
    print(f"Starting FASTA file creation at {datetime.now().isoformat(timespec='seconds')}")

    # Load the sequence file
    seqs_df = pd.read_csv(os.path.join(DATA_DIR, "unique_sequences_with_partitioning.tsv"), sep="\t")

    # Create the FASTA file
    with open(os.path.join(DATA_DIR, "unique_sequences.fasta"), "w") as fasta_file:
        for _, row in seqs_df.iterrows():
            if pd.isna(row["partition"]) or pd.isna(row["cluster"]):
                continue  # Skip rows with missing partition or cluster information

            fasta_file.write(f">{int(row['partition'])}|{row['cluster']}\n")
            fasta_file.write(f"{row['sequence']}\n")

    print(f"FASTA file created at {os.path.join(DATA_DIR, 'unique_sequences.fasta')}")
    print(f"Finished FASTA file creation at {datetime.now().isoformat(timespec='seconds')}")


if __name__ == "__main__":
    main()