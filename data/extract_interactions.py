import pandas as pd

def main():
    interactions = pd.read_csv("/nfs/scratch/pdb_dimers/final_final_filtered_interactions_with_partitions.tsv", sep="\t")

    output = interactions["new_cluster_pair"].str.split(",", expand=True)
    output.rename(columns={0:"protein1", 1:"protein2"}, inplace=True)

    output.to_csv("/nfs/scratch/pdb_dimers/pinder/new_interactions.tsv", sep="\t", index=None)

if __name__ == "__main__":
    main()