import os
from datetime import datetime

import pandas as pd

INTERACTION_FILE = "/nfs/scratch/pdb_dimers/dataset_iterations/06_removed_suspicious_contacts.tsv"
SEQUENCE_FILE = "/nfs/scratch/pdb_dimers/sequences/unique_sequences.tsv"

OUTPUT_DIR = "/nfs/scratch/pdb_dimers/CD-HIT"


SPLITS = ("train", "val", "test")


def cluster_ids_for_split(interaction_df, split):
    split_df = interaction_df.loc[interaction_df["split"] == split, ["cluster1", "cluster2"]]
    return split_df.stack().dropna().astype(str).str.strip().drop_duplicates().tolist()


def write_split_fasta(sequence_df, split, cluster_ids):
    output_path = os.path.join(OUTPUT_DIR, f"{split}.fasta")
    split_sequence_df = (
        pd.DataFrame({"new_cluster_id": cluster_ids})
        .merge(sequence_df, on="new_cluster_id", how="left")
    )

    missing_clusters = split_sequence_df.loc[split_sequence_df["sequence"].isna(), "new_cluster_id"].tolist()
    if missing_clusters:
        raise ValueError(
            f"{len(missing_clusters)} clusters from split {split!r} are missing in {SEQUENCE_FILE}."
            f"Examples: {missing_clusters[:10]}"
        )

    with open(output_path, "w") as fasta_file:
        for index, row in enumerate(split_sequence_df.itertuples(index=False)):
            fasta_file.write(f">{index}_{split}|{row.new_cluster_id}\n")
            fasta_file.write(f"{row.sequence}\n")

    print(f"Wrote {len(split_sequence_df)} sequences to {output_path}")


def main():
    print(f"Starting split FASTA creation at {datetime.now().isoformat(timespec='seconds')}")

    interaction_df = pd.read_csv(INTERACTION_FILE, sep="\t")
    sequence_df = pd.read_csv(SEQUENCE_FILE, sep="\t")

    cluster_columns = interaction_df["new_cluster_pair"].str.split(",", n=1, expand=True)
    interaction_df["cluster1"] = cluster_columns[0].astype(str).str.strip()
    interaction_df["cluster2"] = cluster_columns[1].astype(str).str.strip()

    sequence_df = sequence_df[["new_cluster_id", "sequence"]].copy()
    sequence_df["new_cluster_id"] = sequence_df["new_cluster_id"].astype(str).str.strip()
    sequence_df["sequence"] = sequence_df["sequence"].astype(str).str.strip()
    sequence_df = sequence_df.drop_duplicates(subset="new_cluster_id", keep="first")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for split in SPLITS:
        cluster_ids = cluster_ids_for_split(interaction_df, split)
        write_split_fasta(sequence_df, split, cluster_ids)

    print(f"Finished split FASTA creation at {datetime.now().isoformat(timespec='seconds')}")


if __name__ == "__main__":
    main()
