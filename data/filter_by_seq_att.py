import sys
import os

import pandas as pd

DATA_DIR = "/nfs/scratch/pdb_dimers/"

def main():
    # load dataframes
    interaction_df = pd.read_csv(os.path.join(DATA_DIR, "all_candidate_assemblies.tsv"), sep="\t")
    seqs_df = pd.read_csv(os.path.join(DATA_DIR, "entity_sequences.tsv"), sep="\t")

    # filter interactions by sequence length
    seqs_df["length"] = seqs_df["sequence"].apply(len)

    print(f"Interactions before filtering long sequences: {len(interaction_df)}")

    # Remove sequences longer than x residues
    x = 1250
    long_entities = seqs_df[seqs_df["length"] > x]["entity_name"].tolist()
    filtered_interaction_df = interaction_df[
        ~interaction_df["entity_pair"].apply(lambda pair: any(entity in long_entities for entity in pair.split(",")))
    ]

    print(f"Interactions after filtering long sequences: {len(filtered_interaction_df)}")

    # Remove sequences shorter than y residues
    y = 50
    short_entities = seqs_df[seqs_df["length"] < y]["entity_name"].tolist()
    filtered_interaction_df = filtered_interaction_df[
        ~filtered_interaction_df["entity_pair"].apply(lambda pair: any(entity in short_entities for entity in pair.split(",")))
    ]

    print(f"Interactions after filtering short sequences: {len(filtered_interaction_df)}")

    # Remove sequences with more than z% unknown residues
    z = 15
    unknown_entities = seqs_df[
        seqs_df["sequence"].apply(lambda seq: seq.count("X") / len(seq) * 100 > z)
    ]["entity_name"].tolist()
    filtered_interaction_df = filtered_interaction_df[
        ~filtered_interaction_df["entity_pair"].apply(lambda pair: any(entity in unknown_entities for entity in pair.split(",")))
    ]

    print(f"Interactions after filtering sequences with high unknown residues: {len(filtered_interaction_df)}")

    # Save the filtered interactions to a new TSV file
    filtered_interaction_df.to_csv(os.path.join(DATA_DIR, "filtered_interactions.tsv"), sep="\t", index=False)

if __name__ == "__main__":
    main()