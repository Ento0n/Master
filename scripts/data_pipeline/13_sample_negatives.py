#!/usr/bin/env python3
"""Sample balanced negative protein-protein interactions per split.

The input table contains positive interactions.  The script creates negative
interactions independently for each split, because train/val/test must each keep
their own balance.

Default sampling strategy inside each split:
1. Keep all different-species heteromers.  For same-species heteromers, keep a
   species group only when matching negative heteromers can reproduce the exact
   sequence degree within that group.
2. Create as many negative homomers as possible from sequences that remain in
   the kept positive endpoint pool and are not known positive self-interactors.
   Sequence degree is ignored here.
3. Remove overflowing positive homomers so positive and negative homomer counts
   match.
4. Sample heteromer negatives for the retained heteromer positives.

Every sampled negative pair is checked against all known positive pairs from the
full input table, and against negative pairs already sampled in earlier splits.

Use --keep-positive-homomers to keep all positive homomers.  In that mode the
missing negative rows are filled with extra same-species heteromer negatives,
using positive-vs-negative sequence degree deficits as the guide.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


DATASET_DIR = Path("/nfs/scratch/pdb_dimers/dataset_iterations")
DEFAULT_INPUT = DATASET_DIR / "06_removed_suspicious_contacts.tsv"
CD_HIT_INPUT = DATASET_DIR / "07_cdhit.tsv"


def default_output_path(input_file: Path, keep_positive_homomers: bool) -> Path:
    """Choose the canonical final filename for the selected input and mode."""
    mode = "keep_homomers" if keep_positive_homomers else "strict"
    cd_hit_suffix = "_cd_hit" if input_file.name == CD_HIT_INPUT.name else ""
    # Negative-inclusive datasets are final variants, so their names
    # intentionally have no numeric iteration prefix.
    return input_file.parent / f"balanced_interactions_{mode}{cd_hit_suffix}.tsv"


def canonical_pair(cluster_a: str, cluster_b: str) -> tuple[str, str]:
    """Return an order-independent pair key."""
    return (cluster_a, cluster_b) if cluster_a <= cluster_b else (cluster_b, cluster_a)


def species_key(taxonomy: str, species: str) -> str:
    """Prefer taxonomy ids for species matching and fall back to species names."""
    taxonomy = str(taxonomy).strip()
    species = str(species).strip()
    if taxonomy:
        return taxonomy
    if species:
        return species
    return "unknown"


def most_common_value(values: pd.Series) -> str:
    """Pick the most common non-empty metadata value for a cluster."""
    clean_values = [str(value).strip() for value in values if str(value).strip()]
    if not clean_values:
        return ""
    return Counter(clean_values).most_common(1)[0][0]


def split_pair_values(values: pd.Series, column_name: str) -> tuple[pd.Series, pd.Series]:
    """Split a comma-separated pair column into two stripped endpoint columns."""
    split_values = values.astype(str).str.split(",", n=1, expand=True)
    if split_values.shape[1] == 1:
        split_values[1] = ""

    endpoint_a = split_values[0].fillna("").str.strip()
    endpoint_b = split_values[1].fillna("").str.strip()
    if (endpoint_a == "").any() or (endpoint_b == "").any():
        bad_rows = values.index[(endpoint_a == "") | (endpoint_b == "")].tolist()[:5]
        raise ValueError(f"{column_name} must contain two comma-separated values; bad rows: {bad_rows}")

    return endpoint_a, endpoint_b


def add_helper_columns(df: pd.DataFrame, reset_order: bool = False) -> pd.DataFrame:
    """Add the columns used for balancing without changing the original columns."""
    df = df.copy()

    if reset_order:
        df["_row_order"] = range(len(df))
    elif "_row_order" not in df.columns:
        df["_row_order"] = range(len(df))

    if df.empty:
        df["_cluster_a"] = pd.Series(dtype=str)
        df["_cluster_b"] = pd.Series(dtype=str)
        df["_entity_a"] = pd.Series(dtype=str)
        df["_entity_b"] = pd.Series(dtype=str)
        df["_cluster_100pct_a"] = pd.Series(dtype=str)
        df["_cluster_100pct_b"] = pd.Series(dtype=str)
        df["_pair_key"] = pd.Series(dtype=object)
        df["_species_key_a"] = pd.Series(dtype=str)
        df["_species_key_b"] = pd.Series(dtype=str)
        df["_is_homo"] = pd.Series(dtype=bool)
        df["_species_relation"] = pd.Series(dtype=str)
        return df

    df["_cluster_a"], df["_cluster_b"] = split_pair_values(df["new_cluster_pair"], "new_cluster_pair")
    df["_entity_a"], df["_entity_b"] = split_pair_values(df["entity_pair"], "entity_pair")
    df["_cluster_100pct_a"], df["_cluster_100pct_b"] = split_pair_values(
        df["cluster_pair_100pct"],
        "cluster_pair_100pct",
    )
    df["_pair_key"] = [canonical_pair(a, b) for a, b in zip(df["_cluster_a"], df["_cluster_b"])]
    df["_species_key_a"] = [species_key(t, s) for t, s in zip(df["taxonomy_1"], df["species_1"])]
    df["_species_key_b"] = [species_key(t, s) for t, s in zip(df["taxonomy_2"], df["species_2"])]
    df["_is_homo"] = df["_cluster_a"] == df["_cluster_b"]

    df["_species_relation"] = "different_species"
    df.loc[df["_species_key_a"] == df["_species_key_b"], "_species_relation"] = "same_species"
    df.loc[df["_is_homo"], "_species_relation"] = "homo"
    return df


def load_and_prepare_interactions(input_file: Path) -> tuple[pd.DataFrame, list[str]]:
    """Read the TSV with pandas and add helper columns for balancing."""
    df = pd.read_csv(input_file, sep="\t", dtype=str, keep_default_na=False)

    # Pandas renames an empty first column header to "Unnamed: 0"; keep it round-trippable.
    renamed_columns = {}
    for column in df.columns:
        if column.startswith("Unnamed:"):
            renamed_columns[column] = ""
            break
    if renamed_columns:
        df = df.rename(columns=renamed_columns)

    required_columns = {
        "entity_pair",
        "cluster_pair_100pct",
        "new_cluster_pair",
        "species_1",
        "species_2",
        "taxonomy_1",
        "taxonomy_2",
        "partitions",
        "split",
    }
    missing_columns = sorted(required_columns - set(df.columns))
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")

    original_columns = list(df.columns)
    return add_helper_columns(df), original_columns


def concat_frames(frames: list[pd.DataFrame], template: pd.DataFrame) -> pd.DataFrame:
    """Concatenate non-empty frames while preserving columns for empty results."""
    non_empty_frames = [frame for frame in frames if not frame.empty]
    if non_empty_frames:
        return pd.concat(non_empty_frames, ignore_index=True)
    return template.iloc[0:0].copy()


def degree_counts(df: pd.DataFrame) -> Counter[str]:
    """Count how often each sequence appears as a heteromer endpoint."""
    counts: Counter[str] = Counter()
    counts.update(df["_cluster_a"])
    counts.update(df["_cluster_b"])
    return counts


def pair_same_species_counts(
    counts: Counter[str],
    blocked_pairs: set[tuple[str, str]],
) -> list[tuple[str, str]] | None:
    """Pair sequences from one species while matching the requested sequence degrees."""
    remaining_counts = Counter({cluster: count for cluster, count in counts.items() if count > 0})
    sampled_pairs: list[tuple[str, str]] = []

    # Greedily connect high-degree clusters first.  If no legal pair is possible,
    # this species group cannot be mirrored as negatives without changing degree.
    while remaining_counts:
        active_clusters = sorted(remaining_counts, key=lambda cluster: (-remaining_counts[cluster], cluster))
        chosen_pair: tuple[str, str] | None = None

        for cluster_a in active_clusters:
            for cluster_b in active_clusters:
                if cluster_a == cluster_b:
                    continue
                pair_key = canonical_pair(cluster_a, cluster_b)
                if pair_key not in blocked_pairs:
                    chosen_pair = (cluster_a, cluster_b)
                    break
            if chosen_pair is not None:
                break

        if chosen_pair is None:
            return None

        cluster_a, cluster_b = chosen_pair
        sampled_pairs.append((cluster_a, cluster_b))
        blocked_pairs.add(canonical_pair(cluster_a, cluster_b))

        remaining_counts[cluster_a] -= 1
        remaining_counts[cluster_b] -= 1
        if remaining_counts[cluster_a] == 0:
            del remaining_counts[cluster_a]
        if remaining_counts[cluster_b] == 0:
            del remaining_counts[cluster_b]

    return sampled_pairs


def pair_different_species_counts(
    counts_by_species: dict[str, Counter[str]],
    blocked_pairs: set[tuple[str, str]],
) -> list[tuple[str, str, str, str]] | None:
    """Pair sequences across species while matching sequence degrees."""
    remaining_by_species = {
        species: Counter({cluster: count for cluster, count in counts.items() if count > 0})
        for species, counts in counts_by_species.items()
    }
    remaining_by_species = {
        species: counts for species, counts in remaining_by_species.items() if sum(counts.values()) > 0
    }
    sampled_pairs: list[tuple[str, str, str, str]] = []

    # Greedily connect high-stub species first, then high-degree clusters within
    # those species.  The sampled pairs are always across different species.
    while remaining_by_species:
        species_order = sorted(
            remaining_by_species,
            key=lambda species: (-sum(remaining_by_species[species].values()), species),
        )
        chosen_pair: tuple[str, str, str, str] | None = None

        for species_a in species_order:
            for species_b in species_order:
                if species_a == species_b:
                    continue

                clusters_a = sorted(
                    remaining_by_species[species_a],
                    key=lambda cluster: (-remaining_by_species[species_a][cluster], cluster),
                )
                clusters_b = sorted(
                    remaining_by_species[species_b],
                    key=lambda cluster: (-remaining_by_species[species_b][cluster], cluster),
                )

                for cluster_a in clusters_a:
                    for cluster_b in clusters_b:
                        if cluster_a == cluster_b:
                            continue
                        pair_key = canonical_pair(cluster_a, cluster_b)
                        if pair_key not in blocked_pairs:
                            chosen_pair = (species_a, cluster_a, species_b, cluster_b)
                            break
                    if chosen_pair is not None:
                        break
                if chosen_pair is not None:
                    break
            if chosen_pair is not None:
                break

        if chosen_pair is None:
            return None

        species_a, cluster_a, species_b, cluster_b = chosen_pair
        sampled_pairs.append(chosen_pair)
        blocked_pairs.add(canonical_pair(cluster_a, cluster_b))

        remaining_by_species[species_a][cluster_a] -= 1
        remaining_by_species[species_b][cluster_b] -= 1
        if remaining_by_species[species_a][cluster_a] == 0:
            del remaining_by_species[species_a][cluster_a]
        if remaining_by_species[species_b][cluster_b] == 0:
            del remaining_by_species[species_b][cluster_b]
        if not remaining_by_species[species_a]:
            del remaining_by_species[species_a]
        if species_b in remaining_by_species and not remaining_by_species[species_b]:
            del remaining_by_species[species_b]

    return sampled_pairs


def build_negative_dataframe(
    sampled_pairs: list[tuple[str, str, str, str, str, str, str]],
    metadata: dict[str, dict[str, str]],
    split_value: str,
    original_columns: list[str],
) -> pd.DataFrame:
    """Convert sampled negative pairs back into the input table shape."""
    rows: list[dict[str, str]] = []
    safe_split = str(split_value).replace("/", "_").replace(" ", "_") or "blank"

    for row_number, (cluster_a, cluster_b, taxonomy_a, species_a, taxonomy_b, species_b, dimer_type) in enumerate(
        sampled_pairs
    ):
        row = {column: "" for column in original_columns}
        metadata_a = metadata.get(cluster_a, {})
        metadata_b = metadata.get(cluster_b, {})

        entity_a = metadata_a.get("entity", "")
        entity_b = metadata_b.get("entity", "")
        if not entity_a or not entity_b:
            raise ValueError(
                "Could not find representative entity names for sampled negative pair "
                f"{cluster_a},{cluster_b}"
            )

        cluster_100pct_a = metadata_a.get("cluster_100pct", "")
        cluster_100pct_b = metadata_b.get("cluster_100pct", "")
        if not cluster_100pct_a or not cluster_100pct_b:
            raise ValueError(
                "Could not find representative 100%-cluster ids for sampled negative pair "
                f"{cluster_a},{cluster_b}"
            )

        uniprot_a = metadata_a.get("uniprot", "")
        uniprot_b = metadata_b.get("uniprot", "")

        partition_a = metadata_a.get("partition", "")
        partition_b = metadata_b.get("partition", "")
        if partition_a and partition_a == partition_b:
            partition = partition_a
        elif partition_a or partition_b:
            partition = f"{partition_a}|{partition_b}"
        else:
            partition = ""

        row.update(
            {
                "assembly_id": f"negative_{safe_split}_{dimer_type}_{row_number}",
                "entity_pair": f"{entity_a},{entity_b}",
                "uniprot_pair": f"{uniprot_a}|{uniprot_b}" if uniprot_a or uniprot_b else "",
                "uniprot_1": uniprot_a,
                "uniprot_2": uniprot_b,
                "species_pair": f"{species_a}|{species_b}",
                "taxonomy_pair": f"{taxonomy_a}|{taxonomy_b}",
                "cluster_pair_100pct": f"{cluster_100pct_a},{cluster_100pct_b}",
                "new_cluster_pair": f"{cluster_a},{cluster_b}",
                "dimer_type": dimer_type,
                "species_1": species_a,
                "species_2": species_b,
                "taxonomy_1": taxonomy_a,
                "taxonomy_2": taxonomy_b,
                "partitions": partition,
                "split": split_value,
            }
        )
        rows.append(row)

    negative_df = pd.DataFrame(rows, columns=original_columns)
    if negative_df.empty:
        negative_df = pd.DataFrame(columns=original_columns)
    return add_helper_columns(negative_df, reset_order=True)


def print_stats(label: str, positives: pd.DataFrame, negatives: pd.DataFrame) -> None:
    """Print the balance-relevant counts for a positive/negative pair of tables."""
    pos_homo = int((positives["_species_relation"] == "homo").sum())
    pos_same = int((positives["_species_relation"] == "same_species").sum())
    pos_diff = int((positives["_species_relation"] == "different_species").sum())
    neg_homo = int((negatives["_species_relation"] == "homo").sum())
    neg_same = int((negatives["_species_relation"] == "same_species").sum())
    neg_diff = int((negatives["_species_relation"] == "different_species").sum())

    print(f"\n{label}")
    print(f"  positives: {len(positives)} total | {pos_homo} homomer | {pos_same} same-species heteromer | {pos_diff} different-species heteromer")
    print(f"  negatives: {len(negatives)} total | {neg_homo} homomer | {neg_same} same-species heteromer | {neg_diff} different-species heteromer")


def print_degree_balance(positives: pd.DataFrame, negatives: pd.DataFrame) -> None:
    """Print the largest sequence-degree mismatch between heteromer positives and negatives."""
    positive_counts = degree_counts(positives)
    negative_counts = degree_counts(negatives)
    all_clusters = set(positive_counts) | set(negative_counts)
    max_abs_difference = max(
        (abs(positive_counts.get(cluster, 0) - negative_counts.get(cluster, 0)) for cluster in all_clusters),
        default=0,
    )
    mismatched_clusters = sum(
        1 for cluster in all_clusters if positive_counts.get(cluster, 0) != negative_counts.get(cluster, 0)
    )
    print(f"  heteromer degree balance: {mismatched_clusters} mismatched sequences, max difference {max_abs_difference}")


def process_split(
    split_value: str,
    split_df: pd.DataFrame,
    original_columns: list[str],
    blocked_pairs: set[tuple[str, str]],
    keep_positive_homomers: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Balance positives and negatives for one split."""
    print(f"\n{'=' * 80}")
    print(f"Split: {split_value}")

    # Build one metadata table per sequence in this split.  Negative rows reuse
    # this metadata so sampled pairs still carry source entity, species,
    # taxonomy, and partition annotations.
    endpoint_a = split_df[
        ["_cluster_a", "_entity_a", "_cluster_100pct_a", "taxonomy_1", "species_1", "split", "partitions"]
    ].rename(
        columns={
            "_cluster_a": "cluster",
            "_entity_a": "entity",
            "_cluster_100pct_a": "cluster_100pct",
            "taxonomy_1": "taxonomy",
            "species_1": "species",
        }
    )
    endpoint_b = split_df[
        ["_cluster_b", "_entity_b", "_cluster_100pct_b", "taxonomy_2", "species_2", "split", "partitions"]
    ].rename(
        columns={
            "_cluster_b": "cluster",
            "_entity_b": "entity",
            "_cluster_100pct_b": "cluster_100pct",
            "taxonomy_2": "taxonomy",
            "species_2": "species",
        }
    )
    endpoint_a["uniprot"] = split_df["uniprot_1"] if "uniprot_1" in split_df.columns else ""
    endpoint_b["uniprot"] = split_df["uniprot_2"] if "uniprot_2" in split_df.columns else ""
    endpoints = pd.concat([endpoint_a, endpoint_b], ignore_index=True)
    metadata_df = endpoints.groupby("cluster", sort=False).agg(
        entity=("entity", most_common_value),
        cluster_100pct=("cluster_100pct", most_common_value),
        uniprot=("uniprot", most_common_value),
        taxonomy=("taxonomy", most_common_value),
        species=("species", most_common_value),
        split=("split", most_common_value),
        partition=("partitions", most_common_value),
    )
    metadata = metadata_df.to_dict("index")

    taxonomy_per_cluster = endpoints.groupby("cluster")["taxonomy"].agg(
        lambda values: {str(value).strip() for value in values if str(value).strip()}
    )
    multi_taxonomy_clusters = int((taxonomy_per_cluster.map(len) > 1).sum())
    if multi_taxonomy_clusters:
        print(f"  note: {multi_taxonomy_clusters} clusters have multiple taxonomy annotations; using most common value")

    entities_per_cluster = endpoints.groupby("cluster")["entity"].agg(
        lambda values: {str(value).strip() for value in values if str(value).strip()}
    )
    multi_entity_clusters = int((entities_per_cluster.map(len) > 1).sum())
    if multi_entity_clusters:
        print(f"  note: {multi_entity_clusters} clusters have multiple entity annotations; using most common value")

    empty_negatives = pd.DataFrame(columns=original_columns)
    empty_negatives = add_helper_columns(empty_negatives, reset_order=True)
    print_stats("Initial positives in split", split_df, empty_negatives)

    positive_homomers = split_df[split_df["_is_homo"]].sort_values("_row_order")

    # Same-species heteromers are checked before homomer negatives are drawn.
    # Otherwise a cluster that appears only in a same-species group removed below
    # can still be sampled as a negative homomer and become negative-only in the
    # final balanced output.
    heteromers = split_df[~split_df["_is_homo"]].copy()
    different_species_heteromers = heteromers[heteromers["_species_relation"] == "different_species"].copy()
    same_species_heteromers = heteromers[heteromers["_species_relation"] == "same_species"].copy()

    kept_same_species_groups: list[pd.DataFrame] = []
    removed_same_species_groups: list[pd.DataFrame] = []
    for species, species_group in same_species_heteromers.groupby("_species_key_a", sort=True):
        can_sample_group = pair_same_species_counts(degree_counts(species_group), set(blocked_pairs)) is not None
        if can_sample_group:
            kept_same_species_groups.append(species_group)
        else:
            removed_same_species_groups.append(species_group)

    kept_same_species_heteromers = concat_frames(kept_same_species_groups, split_df)
    removed_same_species_heteromers = concat_frames(removed_same_species_groups, split_df)
    kept_heteromers = concat_frames([kept_same_species_heteromers, different_species_heteromers], split_df)
    kept_heteromers = kept_heteromers.sort_values("_row_order").reset_index(drop=True)

    print("\nAfter removing same-species heteromer groups that cannot be mirrored")
    print(f"  removed same-species positive heteromers: {len(removed_same_species_heteromers)}")
    print_stats("  Heteromer positives kept before homomer sampling", kept_heteromers, empty_negatives)

    # Homomer negatives are self-pairs from the retained positive endpoint pool
    # that are not already positive self-interactors and not already sampled as
    # negatives in another split.
    retained_positive_clusters = set(kept_heteromers["_cluster_a"]) | set(kept_heteromers["_cluster_b"])
    if keep_positive_homomers:
        retained_positive_clusters.update(positive_homomers["_cluster_a"])
        retained_positive_clusters.update(positive_homomers["_cluster_b"])

    eligible_homomer_clusters = sorted(
        cluster for cluster in retained_positive_clusters if canonical_pair(cluster, cluster) not in blocked_pairs
    )
    homomer_negative_count = min(len(positive_homomers), len(eligible_homomer_clusters))

    sampled_homomer_pairs: list[tuple[str, str, str, str, str, str, str]] = []
    for cluster in eligible_homomer_clusters[:homomer_negative_count]:
        cluster_metadata = metadata.get(cluster, {})
        taxonomy = cluster_metadata.get("taxonomy", "")
        species = cluster_metadata.get("species", "")
        sampled_homomer_pairs.append((cluster, cluster, taxonomy, species, taxonomy, species, "homo"))

    negative_homomers = build_negative_dataframe(sampled_homomer_pairs, metadata, split_value, original_columns)
    print("\nAfter homomer negative sampling")
    print(f"  possible negative homomers after global blocking: {len(eligible_homomer_clusters)}")
    print_stats("  Homomer sampling result", positive_homomers, negative_homomers)

    # Either trim positive homomers to match the possible negative homomers, or
    # keep them all and remember which rows must later be filled by heteromers.
    if keep_positive_homomers:
        kept_homomers = positive_homomers.copy()
        overflow_homomers = positive_homomers.iloc[len(negative_homomers):].copy()
        removed_homomers = 0
        print("\nAfter keeping all positive homomers")
        print(f"  positive homomers kept: {len(kept_homomers)}")
        print(f"  homomer negatives missing and filled later with heteromers: {len(overflow_homomers)}")
        print_stats("  Homomer positives kept for weighting experiment", kept_homomers, negative_homomers)
    else:
        kept_homomers = positive_homomers.head(len(negative_homomers)).copy()
        overflow_homomers = positive_homomers.iloc[len(negative_homomers):].copy()
        removed_homomers = len(positive_homomers) - len(kept_homomers)
        print("\nAfter removing overflowing positive homomers")
        print(f"  removed positive homomers: {removed_homomers}")
        print_stats("  Homomer-balanced positives", kept_homomers, negative_homomers)

    positives = concat_frames([kept_homomers, kept_heteromers], split_df)
    positives = positives.sort_values("_row_order").reset_index(drop=True)

    print("\nAfter building positive set kept for negative sampling")
    print_stats("  Positive set kept for negative sampling", positives, negative_homomers)

    # Heteromer negatives are sampled against a split-local copy of the global
    # blocked pairs.  The copy receives homomers first, then each heteromer pair.
    local_blocked_pairs = set(blocked_pairs)
    local_blocked_pairs.update(negative_homomers["_pair_key"])

    sampled_heteromer_pairs: list[tuple[str, str, str, str, str, str, str]] = []

    # Same-species negatives preserve the sequence degree inside each species.
    for species, species_group in kept_same_species_heteromers.groupby("_species_key_a", sort=True):
        sampled_pairs = pair_same_species_counts(degree_counts(species_group), local_blocked_pairs)
        if sampled_pairs is None:
            raise RuntimeError(f"Internal error: same-species group {species} was feasible but sampling failed")

        first_row = species_group.iloc[0]
        if first_row["_species_key_a"] == species:
            taxonomy = first_row["taxonomy_1"]
            species_name = first_row["species_1"]
        else:
            taxonomy = first_row["taxonomy_2"]
            species_name = first_row["species_2"]

        for cluster_a, cluster_b in sampled_pairs:
            sampled_heteromer_pairs.append(
                (cluster_a, cluster_b, taxonomy, species_name, taxonomy, species_name, "hetero")
            )

    # Different-species negatives preserve each sequence degree and keep every
    # endpoint attached to its species group.
    different_species_counts: dict[str, Counter[str]] = defaultdict(Counter)
    species_labels: dict[str, tuple[str, str]] = {}
    for _, row in different_species_heteromers.iterrows():
        different_species_counts[row["_species_key_a"]][row["_cluster_a"]] += 1
        different_species_counts[row["_species_key_b"]][row["_cluster_b"]] += 1
        species_labels[row["_species_key_a"]] = (row["taxonomy_1"], row["species_1"])
        species_labels[row["_species_key_b"]] = (row["taxonomy_2"], row["species_2"])

    sampled_different_species_pairs = pair_different_species_counts(different_species_counts, local_blocked_pairs)
    if sampled_different_species_pairs is None:
        raise RuntimeError(
            f"Could not sample different-species heteromer negatives for split {split_value}; "
            "positive trimming for this case is not implemented because all previous datasets were feasible."
        )

    for species_a, cluster_a, species_b, cluster_b in sampled_different_species_pairs:
        taxonomy_a, species_name_a = species_labels.get(species_a, ("", ""))
        taxonomy_b, species_name_b = species_labels.get(species_b, ("", ""))
        sampled_heteromer_pairs.append(
            (cluster_a, cluster_b, taxonomy_a, species_name_a, taxonomy_b, species_name_b, "hetero")
        )

    # In keep-positive-homomers mode, there are fewer negative homomers than
    # positive homomers.  Fill the missing rows with extra heteromer negatives.
    if keep_positive_homomers:
        extra_needed = len(positives) - len(negative_homomers) - len(sampled_heteromer_pairs)
        extra_sampled = 0
        extra_same_species = 0
        extra_different_species = 0
        fallback_endpoint_count = 0

        if extra_needed > 0:
            cluster_species = {
                cluster: species_key(values.get("taxonomy", ""), values.get("species", ""))
                for cluster, values in metadata.items()
            }
            cluster_labels = {
                cluster: (values.get("taxonomy", ""), values.get("species", ""))
                for cluster, values in metadata.items()
            }
            positive_candidate_clusters = sorted(set(positives["_cluster_a"]) | set(positives["_cluster_b"]))

            # Build the current degree balance after the normal homomer and
            # heteromer negatives.  Filler heteromers should consume sequences
            # that are still underrepresented in the negative set.
            baseline_negative_heteromers = build_negative_dataframe(
                sampled_heteromer_pairs,
                metadata,
                split_value,
                original_columns,
            )
            current_negatives = concat_frames([negative_homomers, baseline_negative_heteromers], negative_homomers)
            positive_degree = degree_counts(positives)
            current_negative_degree = degree_counts(current_negatives)

            degree_deficits: Counter[str] = Counter()
            for cluster, positive_count in positive_degree.items():
                deficit = positive_count - current_negative_degree.get(cluster, 0)
                if deficit > 0:
                    degree_deficits[cluster] = deficit

            deficits_by_species: dict[str, Counter[str]] = defaultdict(Counter)
            for cluster, deficit in degree_deficits.items():
                deficits_by_species[cluster_species.get(cluster, "")][cluster] = deficit

            def add_extra_pair(cluster_a: str, cluster_b: str, forced_partner: bool = False) -> None:
                nonlocal extra_sampled, extra_same_species, extra_different_species, fallback_endpoint_count

                taxonomy_a, species_name_a = cluster_labels.get(cluster_a, ("", ""))
                taxonomy_b, species_name_b = cluster_labels.get(cluster_b, ("", ""))
                sampled_heteromer_pairs.append(
                    (cluster_a, cluster_b, taxonomy_a, species_name_a, taxonomy_b, species_name_b, "hetero")
                )
                local_blocked_pairs.add(canonical_pair(cluster_a, cluster_b))
                current_negative_degree[cluster_a] += 1
                current_negative_degree[cluster_b] += 1

                for cluster in [cluster_a, cluster_b]:
                    if degree_deficits.get(cluster, 0) > 0:
                        degree_deficits[cluster] -= 1
                        if degree_deficits[cluster] == 0:
                            del degree_deficits[cluster]

                    species = cluster_species.get(cluster, "")
                    if deficits_by_species[species].get(cluster, 0) > 0:
                        deficits_by_species[species][cluster] -= 1
                        if deficits_by_species[species][cluster] == 0:
                            del deficits_by_species[species][cluster]

                extra_sampled += 1
                if cluster_species.get(cluster_a) == cluster_species.get(cluster_b):
                    extra_same_species += 1
                else:
                    extra_different_species += 1
                if forced_partner:
                    fallback_endpoint_count += 1

            # First fill with same-species pairs where both endpoints still have
            # positive degree deficit.  This is the cleanest replacement for a
            # missing homomer negative.
            for species in sorted(deficits_by_species):
                while extra_sampled < extra_needed:
                    active_clusters = sorted(
                        deficits_by_species[species],
                        key=lambda cluster: (-deficits_by_species[species][cluster], cluster),
                    )
                    if len(active_clusters) < 2:
                        break

                    chosen_pair: tuple[str, str] | None = None
                    for index, cluster_a in enumerate(active_clusters):
                        for cluster_b in active_clusters[index + 1 :]:
                            if canonical_pair(cluster_a, cluster_b) not in local_blocked_pairs:
                                chosen_pair = (cluster_a, cluster_b)
                                break
                        if chosen_pair is not None:
                            break

                    if chosen_pair is None:
                        break

                    add_extra_pair(chosen_pair[0], chosen_pair[1])

            # If same-species deficit pairs are stuck, keep the anchor endpoint
            # deficit and choose the least overrepresented same-species partner.
            # This preserves the same-species count without creating the strong
            # repeated-partner artefact from the previous implementation.
            while extra_sampled < extra_needed and degree_deficits:
                active_clusters = sorted(degree_deficits, key=lambda cluster: (-degree_deficits[cluster], cluster))
                chosen_pair = None

                for cluster_a in active_clusters:
                    possible_partners = [
                        cluster
                        for cluster in positive_candidate_clusters
                        if cluster != cluster_a
                        and cluster_species.get(cluster) == cluster_species.get(cluster_a)
                        and canonical_pair(cluster_a, cluster) not in local_blocked_pairs
                    ]
                    possible_partners.sort(
                        key=lambda cluster: (
                            current_negative_degree.get(cluster, 0) - positive_degree.get(cluster, 0),
                            -positive_degree.get(cluster, 0),
                            cluster,
                        )
                    )
                    if possible_partners:
                        chosen_pair = (cluster_a, possible_partners[0])
                        break

                if chosen_pair is None:
                    break

                forced_partner = degree_deficits.get(chosen_pair[1], 0) <= 0
                add_extra_pair(chosen_pair[0], chosen_pair[1], forced_partner=forced_partner)

            # Last resort: pair remaining degree deficits across species.  This
            # is only used when no legal same-species heteromer can be found.
            while extra_sampled < extra_needed and degree_deficits:
                active_clusters = sorted(degree_deficits, key=lambda cluster: (-degree_deficits[cluster], cluster))
                chosen_pair = None

                for cluster_a in active_clusters:
                    for cluster_b in active_clusters:
                        if cluster_a == cluster_b:
                            continue
                        if cluster_species.get(cluster_a) == cluster_species.get(cluster_b):
                            continue
                        if canonical_pair(cluster_a, cluster_b) not in local_blocked_pairs:
                            chosen_pair = (cluster_a, cluster_b)
                            break
                    if chosen_pair is not None:
                        break

                if chosen_pair is None:
                    break

                add_extra_pair(chosen_pair[0], chosen_pair[1])

        print("\nAfter filling missing homomer negatives with extra heteromers")
        print(f"  extra heteromers needed: {extra_needed}")
        print(f"  extra heteromers sampled: {extra_sampled}")
        print(f"  extra same-species heteromers: {extra_same_species}")
        print(f"  extra different-species heteromers: {extra_different_species}")
        print(f"  forced same-species partner endpoints used: {fallback_endpoint_count}")

    negative_heteromers = build_negative_dataframe(sampled_heteromer_pairs, metadata, split_value, original_columns)
    negatives = concat_frames([negative_homomers, negative_heteromers], negative_homomers)
    negatives = negatives.reset_index(drop=True)

    print("\nFinal balanced split")
    print_stats("  Final split result", positives, negatives)
    if keep_positive_homomers:
        print("  heteromer degree balance includes filler heteromers replacing untrimmed homomers")
    print_degree_balance(kept_heteromers, negative_heteromers)
    print(f"  removed total positives in split: {len(split_df) - len(positives)}")

    return positives, negatives


def write_outputs(
    output_file: Path,
    positives: pd.DataFrame,
    negatives: pd.DataFrame,
    original_columns: list[str],
    positive_output: Path | None,
    negative_output: Path | None,
) -> None:
    """Write the combined table and optional positive-only / negative-only tables."""
    positives_out = positives[original_columns].copy()
    negatives_out = negatives[original_columns].copy()

    positives_out["label"] = "1"
    positives_out["sample_kind"] = "positive"
    negatives_out["label"] = "0"
    negatives_out["sample_kind"] = "negative"

    combined = pd.concat([positives_out, negatives_out], ignore_index=True)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_file, sep="\t", index=False)

    if positive_output is not None:
        positive_output.parent.mkdir(parents=True, exist_ok=True)
        positives_out.to_csv(positive_output, sep="\t", index=False)

    if negative_output is not None:
        negative_output.parent.mkdir(parents=True, exist_ok=True)
        negatives_out.to_csv(negative_output, sep="\t", index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Sample split-wise balanced negative interactions from a positive interaction TSV."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help=f"Input TSV. Default: {DEFAULT_INPUT}")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "Combined output TSV. By default, choose the canonical balanced_interactions_* "
            "filename from --input and --keep-positive-homomers."
        ),
    )
    parser.add_argument("--positive-output", type=Path, default=None, help="Optional TSV for kept positives only.")
    parser.add_argument("--negative-output", type=Path, default=None, help="Optional TSV for sampled negatives only.")
    parser.add_argument(
        "--keep-positive-homomers",
        action="store_true",
        help="Keep all positive homomers and fill missing negative rows with extra heteromer negatives.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_file = args.output or default_output_path(args.input, args.keep_positive_homomers)
    interactions, original_columns = load_and_prepare_interactions(args.input)
    split_order = list(interactions["split"].drop_duplicates())

    print(f"Loaded {len(interactions)} positive interactions from {args.input}")
    print(f"Found splits in input order: {', '.join(str(split) for split in split_order)}")

    # Block all known positives globally, then add sampled negatives after each
    # split so train/val/test never reuse the same negative interaction.
    blocked_pairs = set(interactions["_pair_key"])
    all_positive_splits: list[pd.DataFrame] = []
    all_negative_splits: list[pd.DataFrame] = []

    for split_value in split_order:
        split_df = interactions[interactions["split"] == split_value].copy()
        positives, negatives = process_split(
            str(split_value),
            split_df,
            original_columns,
            blocked_pairs,
            keep_positive_homomers=args.keep_positive_homomers,
        )
        all_positive_splits.append(positives)
        all_negative_splits.append(negatives)
        blocked_pairs.update(negatives["_pair_key"])

    balanced_positives = concat_frames(all_positive_splits, interactions)
    balanced_negatives = concat_frames(all_negative_splits, interactions)

    print(f"\n{'=' * 80}")
    print_stats("Combined output", balanced_positives, balanced_negatives)
    print_degree_balance(
        balanced_positives[~balanced_positives["_is_homo"]],
        balanced_negatives[~balanced_negatives["_is_homo"]],
    )

    write_outputs(
        output_file,
        balanced_positives,
        balanced_negatives,
        original_columns,
        args.positive_output,
        args.negative_output,
    )
    print(f"\nWrote combined balanced dataset to {output_file}")


if __name__ == "__main__":
    main()
