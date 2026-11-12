import pandas as pd

SEQUENCE_PATH = "/nfs/scratch/pdb_dimers/sequences/entity_sequences.tsv"
INTERACTION_PATH = "/nfs/scratch/pdb_dimers/final_datasets/balanced_interactions_keep_homomers_cd_hit.tsv"

def main():
    sequence_df = pd.read_csv(SEQUENCE_PATH, sep="\t")
    int_df = pd.read_csv(INTERACTION_PATH, sep="\t")

    entity_names = int_df["entity_pair"].str.split(",").explode().str.strip()
    sequence_df = sequence_df[sequence_df["entity_name"].isin(entity_names)].copy()

    sequence_df["counts"] = sequence_df["binary_mask"].str.replace("0", "").str.len()

    print(f"Minimal count of structural annotated AAs fro a sequece is: {sequence_df['counts'].min()}")

if __name__ == "__main__":
    main()