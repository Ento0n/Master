import sys
import os

import pandas as pd
import edlib

DATA_DIR = "/nfs/scratch/pdb_dimers"

def calc_local_identity(seq1: str, seq2: str):
    """
    CD-HIT–style identity via edlib:
      - mode="HW": alignment must cover the whole query (the shorter sequence),
        but can start/end anywhere on the target (the longer sequence).
      - Count identical aligned columns (excluding gaps) and divide by len(shorter).
      - Return (identity, distance).
    """

    # Ensure query is the shorter sequence
    if len(seq1) <= len(seq2):
        query, target = seq1, seq2
    else:
        query, target = seq2, seq1

    # Align: path required so we can reconstruct matches
    res = edlib.align(query, target, mode="HW", task="path")

    # Reconstruct aligned strings
    nice = edlib.getNiceAlignment(res, query, target)
    q_aln = nice["query_aligned"]
    t_aln = nice["target_aligned"]

    # Count exact matches in aligned columns
    matches = 0
    for qc, tc in zip(q_aln, t_aln):
        if qc != "-" and tc != "-" and qc == tc:
            matches += 1

    identity = matches / len(query)

    return identity

def calc_global_identity(seq1: str, seq2: str):
    # calculate sequence identity using edlib
    result = edlib.align(seq1, seq2, mode="NW", task="path")
    edits = result["editDistance"]  # total edits is edit distance
    identity = 1 - (edits / max(len(seq1), len(seq2)))  # identity is 1 - (edits divided by length of longer sequence)

    return identity


def main():
    seqs_df = pd.read_csv(os.path.join(DATA_DIR, "unique_sequences_with_partitioning.tsv"), sep="\t")

    # filter out sequences with no partition or cluster
    seqs_df = seqs_df[~seqs_df["partition"].isna()]
    seqs_df = seqs_df[~seqs_df["cluster"].isna()]

    print(f"Calculating all vs all sequence identity for {len(seqs_df)} sequences")

    # get an all vs all matrix frame
    seqs = seqs_df["sequence"].tolist()
    clusters = seqs_df["cluster"].tolist()

    local_identity_matrix = [[0.0] * len(seqs) for _ in range(len(seqs))]
    global_identity_matrix = [[0.0] * len(seqs) for _ in range(len(seqs))]

    # calculate sequence identity for all pairs
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            seq1 = seqs[i]
            seq2 = seqs[j]

            local_identity = calc_local_identity(seq1, seq2)
            local_identity_matrix[i][j] = local_identity
            local_identity_matrix[j][i] = local_identity  # symmetric matrix

            global_identity = calc_global_identity(seq1, seq2)
            global_identity_matrix[i][j] = global_identity
            global_identity_matrix[j][i] = global_identity  # symmetric matrix
    
    # save the identity matrix to a file
    local_identity_df = pd.DataFrame(local_identity_matrix, index=clusters, columns=clusters)
    local_identity_df.to_csv(os.path.join(os.path.join(DATA_DIR, "own_all_vs_all", "local_all_vs_all_sequence_identity.tsv")), sep="\t")
    global_identity_df = pd.DataFrame(global_identity_matrix, index=clusters, columns=clusters)
    global_identity_df.to_csv(os.path.join(os.path.join(DATA_DIR, "own_all_vs_all", "global_all_vs_all_sequence_identity.tsv")), sep="\t")

    # Create a long format DataFrame for easier analysis
    long_format = []
    partitions = seqs_df["partition"].tolist()
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            long_format.append({
                "cluster1": clusters[i],
                "cluster2": clusters[j],
                "partition1": partitions[i],
                "partition2": partitions[j],
                "local_identity": local_identity_matrix[i][j],
                "global_identity": global_identity_matrix[i][j]
            })
    long_format_df = pd.DataFrame(long_format)
    long_format_df.to_csv(os.path.join(os.path.join(DATA_DIR, "own_all_vs_all", "all_vs_all_sequence_identity_long_format.tsv")), sep="\t", index=False)

    # Save a summary per partition for average and max sequence identity between the partitions
    summary = (
        long_format_df.groupby(["partition1", "partition2"]).agg(
            max_local_identity=("local_identity", "max"),
            avg_local_identity=("local_identity", "mean"),
            max_global_identity=("global_identity", "max"),
            avg_global_identity=("global_identity", "mean"),
            count=("local_identity", "count")
        ).reset_index()
    )
    summary.to_csv(os.path.join(os.path.join(DATA_DIR, "own_all_vs_all", "all_vs_all_sequence_identity_summary.tsv")), sep="\t", index=False)

if __name__ == "__main__":
    main()
