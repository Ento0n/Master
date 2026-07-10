# Dataset Creation Pipeline

This document describes the dataset creation pipeline transcribed from the handwritten workflow and cross-checked against the scripts in `data/` and the Slurm wrappers in `data/slurm/`.

The current pipeline is centered on `/nfs/scratch/pdb_dimers`. Most early scripts use this location as a hard-coded `DATA_DIR`; the later scripts expose command-line arguments for their main inputs and outputs.

## High-Level Flow

```text
RCSB biological assembly search/download
  -> sequence extraction
  -> sequence attribute filtering
  -> exact sequence reclustering
  -> duplicate interaction filtering
  -> all-vs-all sequence identity
  -> graph creation
  -> KaHIP partitioning
  -> partition/split assignment
  -> contact map generation
  -> suspicious contact filtering
  -> optional CD-HIT validation/test filtering
  -> negative sampling
  -> final positive/negative interaction TSVs
```

## Main Paths

| Path | Meaning |
| --- | --- |
| `/nfs/scratch/pdb_dimers` | Main dataset working directory. |
| `/nfs/scratch/pdb_dimers/assemblies` | Downloaded biological assembly mmCIF files. |
| `/nfs/scratch/pdb_dimers/KaHIP` | KaHIP graph mappings and partition outputs. |
| `/nfs/scratch/pdb_dimers/MMseqs2` | MMseqs2 output and sequence-identity graph files. |
| `/nfs/scratch/pdb_dimers/CD-HIT` | Per-split FASTA files and CD-HIT-2D outputs. |
| `/nfs/scratch/pdb_dimers/contact_maps` | Contact map `.npy` files and contact map summary. |

## Step 1: Download RCSB Protein Dimers

Script: `data/download_rcsb_dimers.py`

Typical Slurm wrapper:

```bash
sbatch data/slurm/download_rcsb_dimers.sh
```

Equivalent command:

```bash
python data/download_rcsb_dimers.py --outdir /nfs/scratch/pdb_dimers
```

Purpose:

- Queries the RCSB Search API for biological assemblies.
- Keeps protein-only assemblies with exactly two protein polymer instances.
- Fetches assembly, entity, UniProt, taxonomy, species, resolution, and 100% sequence-cluster metadata through the RCSB GraphQL API.
- Downloads the corresponding assembly mmCIF files unless disabled.

Important arguments:

| Argument | Default | Meaning |
| --- | --- | --- |
| `--outdir` | `rcsb_protein_dimers` | Output directory for TSVs and `assemblies/`. |
| `--batch-size` | `200` | GraphQL metadata batch size. |
| `--no-ligands` | off | Also require zero non-polymer instances in the assembly. |
| `--no-download` | off | Only write metadata TSVs, do not download coordinates. |
| `--limit` | `0` | Debug limit on candidate assemblies before metadata fetch. |

Main outputs:

| Output | Description |
| --- | --- |
| `all_candidate_assemblies.tsv` | Candidate positive dimer metadata. |
| `assemblies/*.cif.gz` | Downloaded biological assembly coordinate files. |

Core columns in `all_candidate_assemblies.tsv`:

- `assembly_id`, `pdb_id`, `assembly_number`
- `entity_pair`
- `uniprot_pair`, `uniprot_1`, `uniprot_2`
- `species_pair`, `species_1`, `species_2`
- `taxonomy_pair`, `taxonomy_1`, `taxonomy_2`
- `dimer_type`
- `cluster_pair_100pct`
- `resolution_best_angstrom`
- `modeled_polymer_monomer_count`
- `experimental_method`, `oligomeric_count`
- `download_url`, `local_filename`

Note: the script still mentions `selected_assemblies.tsv` in its module docstring and final log message, but the representative-selection write is currently commented out. The active output is `all_candidate_assemblies.tsv`.

## Step 2: Extract Entity Sequences

Script: `data/extract_sequences.py`

Command:

```bash
python data/extract_sequences.py
```

Purpose:

- Reads `all_candidate_assemblies.tsv`.
- Opens the corresponding mmCIF assembly files with `gemmi`.
- Extracts one sequence per unique `entity_name`.
- Reads `_pdbx_entity_remapping` and `_struct_asym` metadata to handle remapped entity IDs.
- Builds a binary residue mask per asym chain, where resolved residues are marked with `1` and missing residues with `0`.

Inputs:

- `/nfs/scratch/pdb_dimers/all_candidate_assemblies.tsv`
- `/nfs/scratch/pdb_dimers/assemblies/*.cif.gz`

Output:

- `/nfs/scratch/pdb_dimers/entity_sequences.tsv`

Output columns:

| Column | Meaning |
| --- | --- |
| `entity_name` | PDB/entity identifier such as `10BL_1`. |
| `cluster_id` | Original RCSB 100% sequence-identity cluster. |
| `sequence` | One-letter amino-acid sequence from the mmCIF entity. |
| `binary_mask` | Per-asym resolved-residue mask. |

Implementation caveat: the current file appears to contain an invalid nested quote in the line that builds `binary_mask`. If this script fails to parse, the f-string should be corrected before rerunning the extraction.

## Step 3: Filter by Sequence Attributes

Script: `data/filter_by_seq_att.py`

Command:

```bash
python data/filter_by_seq_att.py
```

Purpose:

- Removes interactions that contain sequences outside the accepted length range.
- Removes interactions with too many unknown residues.

Hard-coded thresholds:

| Filter | Threshold |
| --- | --- |
| Maximum length | `1250` residues |
| Minimum length | `50` residues |
| Maximum unknown residue fraction | `15%` `X` residues |

Inputs:

- `/nfs/scratch/pdb_dimers/all_candidate_assemblies.tsv`
- `/nfs/scratch/pdb_dimers/entity_sequences.tsv`

Output:

- `/nfs/scratch/pdb_dimers/filtered_interactions.tsv`

## Step 4: Recluster Exact Sequences at 100% Identity

Script: `data/cluster_by_100_seq_ident.py`

Command:

```bash
python data/cluster_by_100_seq_ident.py
```

Purpose:

- Restricts `entity_sequences.tsv` to entities still present after sequence-attribute filtering.
- Groups entities by exact sequence string.
- Assigns internal cluster IDs named `fix_<n>_100`.
- Adds `new_cluster_pair` to the interaction table.
- Detects inconsistent cases where one exact sequence maps to multiple original cluster IDs; interactions involving those PDB IDs are removed.

Inputs:

- `/nfs/scratch/pdb_dimers/entity_sequences.tsv`
- `/nfs/scratch/pdb_dimers/filtered_interactions.tsv`

Outputs:

| Output | Description |
| --- | --- |
| `filtered_interactions_with_clusters.tsv` | Interaction table with `new_cluster_pair`. |
| `unique_sequences.tsv` | One row per exact unique sequence with `new_cluster_id`. |

## Step 5: Remove Duplicate Interactions

Script: `data/filter_duplicate_interactions.py`

Command:

```bash
python data/filter_duplicate_interactions.py
```

Purpose:

- Treats `new_cluster_pair` as an unordered pair.
- Keeps only the first representative interaction per cluster pair.
- Removes redundant PDB assemblies that represent the same sequence-pair interaction.

Input:

- `/nfs/scratch/pdb_dimers/filtered_interactions_with_clusters.tsv`

Output:

- `/nfs/scratch/pdb_dimers/final_filtered_interactions.tsv`

## Step 6: Calculate All-vs-All Sequence Identity

Script: `data/calculate_all_vs_all_sequence_identity.py`

Command:

```bash
python data/calculate_all_vs_all_sequence_identity.py --local --global
```

Purpose:

- Computes pairwise sequence identity between every unique sequence used in `final_filtered_interactions.tsv`.
- Uses `edlib` in two modes:
  - local/HW mode, with the shorter sequence as query, for CD-HIT-style local identity.
  - global/NW mode for global identity.
- Requires at least one of `--local` or `--global`. Use both flags to reproduce the old calculation scope.
- Stores only the upper triangle of the symmetric all-vs-all result in compact `float32` NumPy arrays and streams the long-format TSV directly to disk.

Inputs:

- `/nfs/scratch/pdb_dimers/unique_sequences.tsv`
- `/nfs/scratch/pdb_dimers/final_filtered_interactions.tsv`

Outputs:

| Output | Description |
| --- | --- |
| `own_all_vs_all/local_all_vs_all_sequence_identity.triu.npy` | Compact upper-triangular local identities, written when `--local` is set. |
| `own_all_vs_all/global_all_vs_all_sequence_identity.triu.npy` | Compact upper-triangular global identities, written when `--global` is set. |
| `own_all_vs_all/all_vs_all_sequence_identity_cluster_order.tsv` | Cluster order for interpreting triangular arrays. |
| `own_all_vs_all/all_vs_all_sequence_identity_long_format.tsv` | Streamed upper-triangular pairwise identity table, unless `--no-long-format` is set. |

Compatibility note: pass `--write-square-tsv` to also write the old dense square matrix TSVs. This avoids pandas materialization but still creates very large files.

## Step 7: Create Graphs

Script: `data/create_graph.py`

The script supports two graph types.

### Interaction Graph

Command:

```bash
python data/create_graph.py --type interaction
```

Purpose:

- Creates an unweighted METIS graph where nodes are sequence clusters and edges are positive interactions.
- Skips self-connections.
- Writes an integer mapping from `new_cluster_id` to METIS node ID.
- Also writes a `.gexf` version for NetworkX visualization.

Inputs:

- `/nfs/scratch/pdb_dimers/final_filtered_interactions.tsv`

Outputs:

| Output | Description |
| --- | --- |
| `KaHIP/cluster_id_mapping.tsv` | Cluster ID to METIS integer node ID mapping. |
| `KaHIP/cluster.graph` | Unweighted METIS graph. |
| `KaHIP/cluster.graph.gexf` | NetworkX graph export. |

### Sequence-Identity Graph

Command:

```bash
python data/create_graph.py --type sequence_identity
```

Purpose:

- Creates a complete weighted graph where each node is an interaction pair from `final_filtered_interactions.tsv`.
- Edge weight is the maximum local sequence identity between any endpoint of pair A and any endpoint of pair B.
- Weight is written as an integer percentage for METIS/KaHIP.

Inputs:

- `/nfs/scratch/pdb_dimers/final_filtered_interactions.tsv`
- `/nfs/scratch/pdb_dimers/own_all_vs_all/local_all_vs_all_sequence_identity.triu.npy`
- `/nfs/scratch/pdb_dimers/own_all_vs_all/all_vs_all_sequence_identity_cluster_order.tsv`

Implementation notes:

- The triangular NumPy array is loaded with `mmap_mode="r"` and indexed through the cluster-order file.
- Diagonal lookups are treated as `1.0`, because a cluster compared to itself has 100% sequence identity.
- The graph is written row-by-row instead of materializing the complete adjacency list in memory.

Output:

- `/nfs/scratch/pdb_dimers/MMseqs2/sequence_identity.graph`

METIS header detail: the script writes `001`, indicating edge weights are present.

## Step 8: Partition the Sequence-Identity Graph with KaHIP

Slurm wrapper: `data/slurm/get_sequence_identity_partitions.sh`

Command used in the wrapper:

```bash
/nfs/home/students/a.spannagl/.local/kahip64/bin/kaffpa \
  /nfs/scratch/pdb_dimers/MMseqs2/sequence_identity.graph \
  --k=10 \
  --preconfiguration=strong \
  --imbalance=3 \
  --output_filename=/nfs/scratch/pdb_dimers/KaHIP/seq_ident_partitions_strong_output.txt
```

Purpose:

- Partitions the interaction-pair graph into 10 groups.
- The sequence-identity edge weights encourage similar interactions to stay together.
- The resulting partition IDs are later mapped to train, validation, and test splits.

Important KaHIP arguments:

| Argument | Meaning |
| --- | --- |
| `--k=10` | Produce 10 partitions. |
| `--preconfiguration=strong` | Use KaHIP's higher-quality preset. |
| `--imbalance=3` | Allow 3% partition imbalance. |
| `--output_filename` | Write one partition ID per graph node. |

Output:

- `/nfs/scratch/pdb_dimers/KaHIP/seq_ident_partitions_strong_output.txt`

Optional summary:

```bash
python data/create_partition_seq_ident_summary.py
```

This reads the partition file and `MMseqs2/sequence_identity.graph`, then writes:

- `/nfs/scratch/pdb_dimers/KaHIP/partition_sequence_identity_summary_strong.tsv`

The summary reports max/average sequence identity across partition pairs.

## Step 9: Attach Partitions and Assign Splits

This step is currently implemented manually in `notebooks/check_kahip_output.ipynb`.

The notebook loads:

```python
int_df = pd.read_csv("/nfs/scratch/pdb_dimers/final_filtered_interactions.tsv", sep="\t")
partitions = pd.read_csv(
    "/nfs/scratch/pdb_dimers/KaHIP/seq_ident_partitions_strong_output.txt",
    header=None,
)
```

Then attaches the partition IDs:

```python
partitions = partitions[0]
int_df["partitions"] = partitions
```

The observed split mapping is:

```python
int_df["split"] = int_df["partitions"].map({
    0: "val",
    1: "train",
    2: "train",
    3: "train",
    4: "train",
    5: "train",
    6: "train",
    7: "test",
    8: "train",
    9: "train",
})
```

Output:

- `/nfs/scratch/pdb_dimers/final_filtered_interactions_with_partitions.tsv`

Note: this notebook writes the TSV without `index=False`, so the file has a leading index column. Some downstream scripts explicitly tolerate this by renaming `Unnamed:` columns.

## Step 10: Generate Contact Maps

Script: `data/generate_contact_maps.py`

Typical Slurm wrapper:

```bash
sbatch data/slurm/generate_contact_maps.sh
```

Equivalent command:

```bash
python data/generate_contact_maps.py \
  --selection-tsv /nfs/scratch/pdb_dimers/final_filtered_interactions_with_partitions.tsv \
  --input-dir /nfs/scratch/pdb_dimers/assemblies \
  --output-dir /nfs/scratch/pdb_dimers/contact_maps \
  --force
```

Purpose:

- Builds one inter-chain interface contact map per selected positive dimer.
- Uses full sequence coordinates, including unresolved residues.
- Stores unknown residue pairs explicitly.

Contact map encoding:

| Value | Meaning |
| --- | --- |
| `-1` | Unknown pair because at least one residue is unresolved. |
| `0` | Known non-contact. |
| `1` | Known contact. |

Important arguments:

| Argument | Default | Meaning |
| --- | --- | --- |
| `--selection-tsv` | `/nfs/scratch/pdb_dimers/final_filtered_interactions_with_partitions.tsv` | Interaction rows to process. |
| `--input-dir` | `/nfs/scratch/pdb_dimers/assemblies` | Folder with `.cif` or `.cif.gz` files. |
| `--output-dir` | `/nfs/scratch/pdb_dimers/contact_maps` | Folder for `.npy` contact maps. |
| `--summary-tsv` | `<output-dir>/contact_map_summary.tsv` | Optional custom summary output path. |
| `--threshold` | `8.0` | Contact distance threshold in Angstrom. |
| `--map-type` | `interface` | Only interface maps are currently supported. |
| `--representative-atom` | `cbeta` | Use C-beta with C-alpha fallback, or force C-alpha. |
| `--chunk-size` | `1024` | Rows per SciPy distance calculation chunk. |
| `--allow-ambiguous-chain-selection` | off | Use first matching chains instead of failing when selection is ambiguous. |
| `--force` | off | Regenerate maps that already exist. |
| `--limit` | none | Process only the first N TSV rows for testing. |

Outputs:

| Output | Description |
| --- | --- |
| `contact_maps/*.npy` | One `int8` contact map per assembly. |
| `contact_maps/contact_map_summary.tsv` | Audit trail of chain selection, residues, missing residues, contacts, unknown pairs, and shape. |

## Step 11: Filter Suspicious Contact Maps

Script: `data/filter_suspicious_contacts.py`

Command:

```bash
python data/filter_suspicious_contacts.py
```

Purpose:

- Reads contact-map summary statistics.
- Computes contact-map sparsity as `contacts / (residues_1 * residues_2 - unknown_pairs)`.
- Removes structures with zero sparsity.
- Hard-codes additional problematic PDB IDs `8UG2` and `8UG0`.

Inputs:

- `/nfs/scratch/pdb_dimers/final_filtered_interactions_with_partitions.tsv`
- `/nfs/scratch/pdb_dimers/contact_maps/contact_map_summary.tsv`

Output:

- `/nfs/scratch/pdb_dimers/final_final_filtered_interactions_with_partitions.tsv`

Note: the current script writes with the default pandas index, so the output also contains a leading index column.

## Step 12: Optional FASTA and Similarity Diagnostics

### Whole-Dataset FASTA

Script: `data/create_fasta.py`

Command:

```bash
python data/create_fasta.py
```

Purpose:

- Reads `unique_sequences_with_partitioning.tsv`.
- Writes one FASTA record per sequence that has partition information.
- FASTA headers use `<partition>|<cluster>`.

Output:

- `/nfs/scratch/pdb_dimers/unique_sequences.fasta`

### MMseqs2 All-vs-All Diagnostic

Slurm wrapper: `data/slurm/run_mmseqs.sh`

Command:

```bash
mmseqs easy-search \
  /nfs/scratch/pdb_dimers/unique_sequences.fasta \
  /nfs/scratch/pdb_dimers/unique_sequences.fasta \
  /nfs/scratch/pdb_dimers/MMseqs2/hits.tsv \
  /nfs/scratch/pdb_dimers/MMseqs2/tmp \
  --format-output query,target,pident,nident,alnlen,qcov,tcov,evalue,bits \
  --min-seq-id 0.0 \
  -c 0.0 \
  --cov-mode 0 \
  -s 7.5 \
  --max-seqs 100000
```

Summary script:

```bash
python data/create_mmseqs2_summary.py
```

Output:

- `/nfs/scratch/pdb_dimers/MMseqs2/partition_identity_summary.tsv`

This diagnostic summarizes sequence identity between partitions based on MMseqs2 hits.

## Step 13: Optional CD-HIT-2D Split Filtering

The handwritten workflow contains a branch labeled "without CD-HIT"; the repository supports both with-CD-HIT and without-CD-HIT final datasets.

### Create Per-Split FASTA Files

Script: `data/create_fasta_per_split.py`

Command:

```bash
python data/create_fasta_per_split.py
```

Purpose:

- Reads positive interactions after contact filtering.
- Extracts unique cluster IDs for each split.
- Writes `train.fasta`, `val.fasta`, and `test.fasta` with headers `<index>_<split>|<new_cluster_id>`.

Inputs:

- `/nfs/scratch/pdb_dimers/final_final_filtered_interactions_with_partitions.tsv`
- `/nfs/scratch/pdb_dimers/unique_sequences.tsv`

Outputs:

| Output | Description |
| --- | --- |
| `CD-HIT/train.fasta` | Train split unique sequences. |
| `CD-HIT/val.fasta` | Validation split unique sequences. |
| `CD-HIT/test.fasta` | Test split unique sequences. |

Wrapper caveat: `data/slurm/create_fasta_per_split.sh` currently calls `create_fasta.py`, not `create_fasta_per_split.py`. Use the Python command above or fix the wrapper before submitting it.

### Run CD-HIT-2D

Train vs test:

```bash
cd-hit-2d \
  -i /nfs/scratch/pdb_dimers/CD-HIT/train.fasta \
  -i2 /nfs/scratch/pdb_dimers/CD-HIT/test.fasta \
  -o /nfs/scratch/pdb_dimers/CD-HIT/train_test.out \
  -n 2 \
  -c 0.4 \
  -s2 0 \
  -S2 999999
```

Train vs validation:

```bash
cd-hit-2d \
  -i /nfs/scratch/pdb_dimers/CD-HIT/train.fasta \
  -i2 /nfs/scratch/pdb_dimers/CD-HIT/val.fasta \
  -o /nfs/scratch/pdb_dimers/CD-HIT/train_val.out \
  -n 2 \
  -c 0.4 \
  -s2 0 \
  -S2 999999
```

Argument meaning:

| Argument | Meaning |
| --- | --- |
| `-i` | First database, here the train sequences. |
| `-i2` | Second database, here test or validation. |
| `-o` | Output FASTA-like keep-list. |
| `-n 2` | Word size appropriate for the low identity threshold. |
| `-c 0.4` | 40% identity threshold. |
| `-s2 0` | Do not require a minimum length ratio for the second input. |
| `-S2 999999` | Allow very large length differences for the second input. |

### Remove Similar Val/Test Sequences

Script: `data/remove_similar_sequences.py`

Command:

```bash
python data/remove_similar_sequences.py
```

Purpose:

- Reads CD-HIT output headers.
- Keeps validation rows only when both cluster IDs in `new_cluster_pair` are present in `train_val.out`.
- Keeps test rows only when both cluster IDs in `new_cluster_pair` are present in `train_test.out`.
- Leaves other splits unchanged.

Important arguments:

| Argument | Default | Meaning |
| --- | --- | --- |
| `--input` | `/nfs/scratch/pdb_dimers/final_final_filtered_interactions_with_partitions.tsv` | Interaction TSV to filter. |
| `--output` | `/nfs/scratch/pdb_dimers/final_final_filtered_interactions_with_partitions_cd_hit.tsv` | Filtered TSV. |
| `--val-clusters` | `/nfs/scratch/pdb_dimers/CD-HIT/train_val.out` | Allowed validation clusters. |
| `--test-clusters` | `/nfs/scratch/pdb_dimers/CD-HIT/train_test.out` | Allowed test clusters. |

Output:

- `/nfs/scratch/pdb_dimers/final_final_filtered_interactions_with_partitions_cd_hit.tsv`

## Step 14: Sample Negative Interactions

Script: `data/sample_negatives.py`

Without CD-HIT filtering:

```bash
python data/sample_negatives.py \
  --output /nfs/scratch/pdb_dimers/balanced_interactions_strict.tsv
```

With all positive homomers kept:

```bash
python data/sample_negatives.py \
  --output /nfs/scratch/pdb_dimers/balanced_interactions_keep_homomers.tsv \
  --keep-positive-homomers
```

After CD-HIT filtering:

```bash
python data/sample_negatives.py \
  --input /nfs/scratch/pdb_dimers/final_final_filtered_interactions_with_partitions_cd_hit.tsv \
  --output /nfs/scratch/pdb_dimers/balanced_interactions_strict_cd_hit.tsv
```

After CD-HIT filtering with all positive homomers kept:

```bash
python data/sample_negatives.py \
  --input /nfs/scratch/pdb_dimers/final_final_filtered_interactions_with_partitions_cd_hit.tsv \
  --output /nfs/scratch/pdb_dimers/balanced_interactions_keep_homomers_cd_hit.tsv \
  --keep-positive-homomers
```

Purpose:

- Starts from positive interactions.
- Samples negative interactions independently per split so train, validation, and test stay balanced.
- Globally blocks all known positive pairs.
- Also blocks already sampled negative pairs so different splits do not reuse the same negative interaction.

Default balancing strategy:

1. Sample as many negative homomers as possible from split-local sequences that are not known positive self-interactors.
2. In strict mode, trim overflowing positive homomers so positive and negative homomer counts match.
3. Keep all different-species heteromers and sample matching negative heteromers.
4. Keep a same-species heteromer group only when the negative sampler can reproduce the exact sequence degree inside that species group.

`--keep-positive-homomers` changes step 2: all positive homomers are retained, and missing negative rows are filled with extra heteromer negatives guided by positive-vs-negative degree deficits.

Important arguments:

| Argument | Default | Meaning |
| --- | --- | --- |
| `--input` | `/nfs/scratch/pdb_dimers/final_final_filtered_interactions_with_partitions.tsv` | Positive interaction TSV. |
| `--output` | `/nfs/scratch/pdb_dimers/balanced_positive_negative_interactions.tsv` | Combined positive/negative output. |
| `--positive-output` | none | Optional kept-positive-only TSV. |
| `--negative-output` | none | Optional sampled-negative-only TSV. |
| `--keep-positive-homomers` | off | Keep all positive homomers and fill with extra heteromer negatives. |

Final output columns:

- All input interaction columns.
- `label`: `1` for positives, `0` for negatives.
- `sample_kind`: `positive` or `negative`.

Observed produced outputs in the current scratch directory:

| Output | Description |
| --- | --- |
| `balanced_interactions_strict.tsv` | Strict positive/negative balancing without CD-HIT filtering. |
| `balanced_interactions_keep_homomers_cd_hit.tsv` | Keep-homomers variant after CD-HIT filtering. |
| `balanced_interactions_strict_cd_hit.tsv` | Strict balancing after CD-HIT filtering. |

## Optional Downstream Embeddings

Script: `data/generate_embeddings.py`

Command example:

```bash
python data/generate_embeddings.py \
  --input /nfs/scratch/pdb_dimers/unique_sequences.tsv \
  --output-dir /nfs/scratch/pdb_dimers/embeddings \
  --models all \
  --device auto
```

Purpose:

- Generates per-residue and mean embeddings for `unique_sequences.tsv`.
- Supports ESM2 `facebook/esm2_t33_650M_UR50D` and ProstT5 `Rostlab/ProstT5`.
- Saves tensors by `new_cluster_id`.

Important arguments:

| Argument | Default | Meaning |
| --- | --- | --- |
| `--models` | `all` | One or more of `all`, `esm2`, `prostt5`. |
| `--device` | `auto` | Torch device. |
| `--dtype` | `auto` | Model compute dtype. |
| `--storage-dtype` | `float16` | Saved tensor dtype. |
| `--batch-size` | `2` | Maximum chunks per forward pass. |
| `--max-batch-residues` | `2048` | Approximate residues per forward pass. |
| `--esm2-max-residues` | `1022` | ESM2 chunk size; `0` disables chunking. |
| `--prostt5-max-residues` | `1022` | ProstT5 chunk size; `0` disables chunking. |
| `--chunk-overlap` | `0` | Overlap between chunks; overlaps are averaged. |
| `--local-files-only` | off | Only use already cached Hugging Face models. |
| `--force` | off | Regenerate existing outputs. |
| `--limit` | none | Debug row limit. |
| `--start-index`, `--end-index` | none | Row range for array jobs. |

## Operational Notes

- Many scripts assume they are run in an environment with `pandas`, `gemmi`, `edlib`, `biopython`, `scipy`, `networkx`, and the external tools KaHIP, MMseqs2, and CD-HIT installed.
- Most scripts write into `/nfs/scratch/pdb_dimers`; rerunning a step can overwrite downstream inputs.
- Several scripts do not expose CLI arguments for thresholds or paths. Change the constants in the script or add CLI arguments before running the pipeline on a different directory.
- KaHIP partitioning and split assignment are order-sensitive: the partition output must have exactly one row per node in `sequence_identity.graph` and must match the row order of the interaction table used to create that graph.
- The CD-HIT branch is optional. For the "without CD-HIT" dataset, sample negatives directly from `final_final_filtered_interactions_with_partitions.tsv`.
- For the "with CD-HIT" dataset, run per-split FASTA creation, CD-HIT-2D, `remove_similar_sequences.py`, then negative sampling from the `_cd_hit.tsv` file.
