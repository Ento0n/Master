import pandas as pd

def main():
    interactions = pd.read_csv("/nfs/scratch/pdb_dimers/final_final_filtered_interactions_with_partitions.tsv", sep="\t")

    output = interactions[["uniprot_1", "uniprot_2"]].copy()

    output["uniprot_1"] = output["uniprot_1"].str.split(";").str[0]
    output["uniprot_2"] = output["uniprot_2"].str.split(";").str[0]

    output.rename(columns={"uniprot_1":"protein1", "uniprot_2":"protein2"}, inplace=True)

    output.dropna(subset=["protein1", "protein2"], inplace=True)

    output.to_csv("/nfs/scratch/pdb_dimers/pinder/pdb_interactions.tsv", sep="\t", index=None)

if __name__ == "__main__":
    main()