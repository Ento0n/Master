import pandas as pd

INTERACTION_FILE = "/nfs/scratch/pdb_dimers/final_filtered_interactions_with_partitions.tsv"
SUMMARY_FILE = "/nfs/scratch/pdb_dimers/contact_maps/contact_map_summary.tsv"

def main():
    int_df = pd.read_csv(INTERACTION_FILE, sep="\t")
    summary_df = pd.read_csv(SUMMARY_FILE, sep="\t")

    # 2 PDB IDs with max sparsity, look completely wrong...
    blacklist = ["8UG2", "8UG0"]

    summary_df["sparsity"] = summary_df["contacts"] / (summary_df["residues_1"] * summary_df["residues_2"] - summary_df["unknown_pairs"])
    blacklist.extend(summary_df[summary_df["sparsity"] == 0.0]["pdb_id"])

    print(f"Complete blacklist: {blacklist}")
    print(f"Size before blacklisting: {len(int_df)}")

    int_df = int_df[~int_df["pdb_id"].isin(blacklist)]

    print(f"Size after blacklisting: {len(int_df)}")

    int_df.to_csv("/nfs/scratch/pdb_dimers/final_final_filtered_interactions_with_partitions.tsv", sep="\t")


if __name__ == "__main__":
    main()
