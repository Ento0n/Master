# Protein Interaction and Contact-Map Prediction

Research code for a master's thesis on predicting protein–protein interactions and inter-chain residue contacts. The project builds datasets from RCSB PDB biological assemblies, generates protein language model embeddings, and compares a one-hot baseline, an adapted D-SCRIPT model, and Transformer-based QUERY-PATCH models.

The repository covers data preparation, training, validation-based model and threshold selection, and result analysis. Datasets, pretrained embedding models, and trained checkpoints are stored separately and are not included in a clone.

## Models

| Model | Input and approach | Predictions |
| --- | --- | --- |
| Fully connected baseline | Concatenated, padded one-hot encodings of both protein sequences, followed by a multilayer perceptron. | Interaction probability |
| Adapted D-SCRIPT | Cached residue embeddings, a shared projection, a convolutional contact module, and an interaction pooling head. | Contact-map and interaction probability |
| QUERY-PATCH | Learned Transformer queries identify paired interface patches, which are combined with bilinear or convolutional residue compatibility. | Contact-map and interaction probability |

The adapted D-SCRIPT model is implemented locally in [scripts/models/dscript.py](scripts/models/dscript.py). QUERY-PATCH supports three named variants:

| Variant | Training arguments |
| --- | --- |
| QP-Additive | `--model query_patch --query-evidence-mode additive --compatibility-head bilinear` |
| QP-Bilin-IC | `--model query_patch --query-evidence-mode gated --compatibility-head bilinear` |
| QP-Conv-IC | `--model query_patch --query-evidence-mode gated --compatibility-head dscript` |

Additive mode sums query and compatibility evidence and uses the main interaction/contact losses. Gated mode uses query evidence to gate compatibility probabilities and adds branch supervision. Contact models support both the normal D-SCRIPT interaction head and maximum-contact pooling through `--int-mod-type normal|max`.

These flags select architectures; reproducing a particular experiment also requires its data, hyperparameters, seed, and checkpoint-selection settings.

## Repository layout

```text
scripts/
  data_pipeline/       Numbered dataset preparation and embedding scripts
  models/              Fully connected, adapted D-SCRIPT, and QUERY-PATCH architectures
  training_pipeline/   Lightning training, losses, metrics, and threshold selection
  analysis/            Dataset summaries and aggregate experiment reports
  figures/             Contact-map comparison figures
  misc/                Shared helpers, preparation, and benchmark utilities
docs/                  Dataset, metric, and model-history documentation
results/               Exported baseline, interaction, and contact summary CSVs
figures/               Exported contact-map comparisons and dataset plots
```

Training artifacts, W&B logs, and exploratory notebooks may also exist in the original workspace. They are not all tracked in Git.

## Setup

Run the commands below from the repository root. The project uses scripts directly and currently has no dependency lockfile or installable Python package.

```bash
git clone https://github.com/Ento0n/Master.git
cd Master

python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip

# Training, evaluation, and plotting
python -m pip install torch lightning torchmetrics numpy pandas wandb matplotlib

# Additional dependencies for dataset preparation and embedding generation
python -m pip install gemmi biopython scipy edlib transformers sentencepiece tqdm
```

Use a PyTorch build compatible with your machine and CUDA installation when running on a GPU. The installation commands above are unpinned; reproducing an existing experiment requires its original package versions as well as its data and training settings. Record the environment used for new runs alongside their outputs.

Preparing the full dataset also requires **KaHIP** (`kaffpa`) for graph partitioning and **CD-HIT** (`cd-hit-2d`) for the optional sequence-identity filtering. Some auxiliary workflows use MMseqs2. Downloading structures and uncached embedding model weights requires network access.

Check the training entrypoint without loading data or starting a run:

```bash
python -m scripts.training_pipeline.execute_pipeline --help
```

## Data preparation

The original dataset root is `/nfs/scratch/pdb_dimers`. Several preparation scripts hard-code this location; later scripts accept input/output paths. Adapt the paths in the Python and shell scripts under `scripts/data_pipeline/` before rebuilding the dataset elsewhere.

The current workflow is:

| Stage | Scripts in `scripts/data_pipeline/` | Purpose |
| --- | --- | --- |
| 1–2 | `01_download_rcsb_dimers.py`, `02_extract_sequences.py` | Download protein dimer assemblies and extract entity sequences and resolved-residue masks. |
| 3–5 | `03_filter_by_seq_att.py` through `05_filter_duplicate_interactions.py` | Filter sequence attributes, cluster identical sequences, and remove duplicate sequence-pair interactions. |
| 6–8 | `06_calculate_all_vs_all_sequence_identity.py`, `07_create_graph.py`, `08_KaHIP_partitioning.sh` | Compute sequence similarity, construct its graph, and partition with KaHIP. |
| 9 | `09_assign_splits.txt` | Records the manual partition-to-split assignment step. |
| 10–11 | `10_generate_contact_maps.py`, `11_filter_suspicious_contacts.py` | Generate structural contact labels and remove suspicious maps. |
| 12, optional | `12_1_create_fasta_per_split.py` through `12_4_remove_similar_sequences.py` | Filter validation/test sequences against training sequences with CD-HIT-2D. |
| 13 | `13_sample_negatives.py` | Sample negatives within each split while excluding known positives and duplicate negative pairs across splits. |
| 14 | `14_generate_embeddings.py` | Cache ESM2 or ProstT5 residue embeddings and mean embeddings. |

**Split assignment is a manual step.** Its reference notebook, `notebooks/check_kahip_output.ipynb`, is local and not tracked. A fresh clone therefore does not provide a fully automated rebuild. Supply the split-assigned table `dataset_iterations/05_KaHIP_partitions.tsv` before continuing with contact generation.

The resulting data layout is:

```text
pdb_dimers/
  assemblies/                         Biological assembly mmCIF files
  dataset_iterations/                 Intermediate interaction tables
  sequences/
    entity_sequences.tsv              Entity IDs and raw sequences
    unique_sequences.tsv              Exact-sequence cluster IDs and sequences
  own_all_vs_all/                      Pairwise sequence-identity arrays
  graphs/                             Sequence-identity graph and mappings
  KaHIP/                              Partition assignments
  CD-HIT/                             Split FASTA files and filtering outputs
  contact_maps/data/                   Structural contact-maps (.npy)
  final_datasets/                      Balanced interaction TSVs
  embeddings/esm2_t33_650M/            ESM2 tensors and manifest.tsv
  embeddings/prostt5/                  ProstT5 tensors and manifest.tsv
```

Negative sampling provides `strict` and `keep_homomers` variants, each with an optional `_cd_hit` suffix. Strict sampling balances positive and negative homomer counts by removing excess positive homomers. The keep-homomers mode retains positive homomers and fills the remaining negative count with heteromers. Sampled negatives exclude known positive pairs in the input; they are not experimentally confirmed non-interactions.

Once the preceding preparation steps are complete, generate the strict CD-HIT-filtered dataset and ESM2 embeddings as follows. Set `DATA_ROOT` to an **absolute path** to your prepared data and reuse it in the training examples below:

```bash
DATA_ROOT=/path/to/pdb_dimers

python scripts/data_pipeline/13_sample_negatives.py \
  --input "$DATA_ROOT/dataset_iterations/07_cdhit.tsv" \
  --output "$DATA_ROOT/final_datasets/balanced_interactions_strict_cd_hit.tsv"

python scripts/data_pipeline/14_generate_embeddings.py \
  --input "$DATA_ROOT/sequences/unique_sequences.tsv" \
  --output-dir "$DATA_ROOT/embeddings" \
  --models esm2 \
  --device cuda
```

An explicit negative-sampling output path is important: its inferred output otherwise goes beside the input table. Embedding generation supports overlapping chunks for long sequences and writes `manifest.tsv`, per-residue tensors, and mean tensors. Training uses the per-residue tensors.

### Training input formats

All interaction files are tab-separated tables. Provide nonempty `train` and `val` splits; include `test` for held-out evaluation. The current trainer registers early stopping against a validation metric, so a training-only table is insufficient for the documented workflow.

| Input | Required columns or contents |
| --- | --- |
| All interaction TSVs | `split` (`train`, `val`, `test`), `label` (`1` interacting, `0` sampled negative), `dimer_type` (`homo`, `hetero`) |
| Baseline interaction TSV | Additionally `entity_pair`, containing two comma-separated entity IDs |
| Baseline sequence TSV | `entity_name`, `sequence` |
| Contact-model interaction TSV | Additionally `new_cluster_pair` (two comma-separated cluster IDs), `pdb_id`, `assembly_number` |
| Embedding manifest TSV | `new_cluster_id`, `per_residue_path` |
| Embedding tensor | A raw PyTorch tensor of shape `[sequence_length, embedding_dim]` |
| Contact-map | An integer NumPy array of shape `[length_1, length_2]`, aligned with the pair's sequences |

Embedding paths are read as stored, rather than resolved relative to the manifest. Use absolute paths or paths valid from the working directory. The default embedding dimension is 1280 for ESM2; set `--embedding-dim` to the tensor width when using other embeddings.

Contact files are named `<lowercase_pdb_id>-assembly<assembly_number>.npy`, for example `10bl-assembly1.npy`. Labels are `1 = contact`, `0 = known non-contact`, and `-1 = unresolved`. Unknown cells and batch padding are excluded from contact losses and metrics. Every positive example must have a contact-map; synthetic negative pairs may leave `pdb_id` and `assembly_number` empty and have no structural map.

The contact generator uses an 8 Å distance cutoff by default, based on Cβ coordinates with Cα for glycine or missing Cβ atoms. The optional training flag `--negative-contact-maps` supplies all-zero targets for pairs without a real map. Those synthetic maps contribute to contact losses but are excluded from contact AUPRC.

## Training

The following examples assume `DATA_ROOT` has been set as above. They write to `runs/` and use offline W&B logging. The trainer always creates a W&B logger for the `master` project; for online tracking, authenticate with `wandb login` and omit `WANDB_MODE=offline`.

Default paths target the original workspace. Supply `--interactions`, `--output-dir`, and the model-specific sequence or embedding/contact paths explicitly, as shown below.

### Adapted D-SCRIPT model

```bash
WANDB_MODE=offline python -m scripts.training_pipeline.execute_pipeline \
  --model dscript \
  --interactions "$DATA_ROOT/final_datasets/balanced_interactions_strict_cd_hit.tsv" \
  --embedding-manifest "$DATA_ROOT/embeddings/esm2_t33_650M/manifest.tsv" \
  --contact-map-dir "$DATA_ROOT/contact_maps/data" \
  --loss-type balanced_bce \
  --interaction-loss-lambda 0.5 \
  --batch-size 8 \
  --max-epochs 100 \
  --seed 42 \
  --accelerator gpu --devices 1 \
  --output-dir ./runs --output-subdir dscript_esm2
```

For QUERY-PATCH, replace `--model dscript` with a variant's arguments from the model table and choose a different `--output-subdir`. Current QUERY-PATCH defaults are 64 learned queries, 4 attention heads, one decoder layer, gated evidence, and convolutional compatibility. `--projection-dim` must be divisible by `--query-heads`.

### One-hot baseline

```bash
WANDB_MODE=offline python -m scripts.training_pipeline.execute_pipeline \
  --model fully_connected \
  --interactions "$DATA_ROOT/final_datasets/balanced_interactions_strict_cd_hit.tsv" \
  --sequence-file "$DATA_ROOT/sequences/entity_sequences.tsv" \
  --selection-monitor val_loss \
  --selection-tiebreaker val_interaction_auprc \
  --max-epochs 100 \
  --seed 42 \
  --accelerator cpu --devices 1 \
  --output-dir ./runs --output-subdir one_hot
```

The baseline requires an explicit `--selection-tiebreaker` because the shared default is a contact metric that this model does not produce. It pads or truncates each sequence to `--max-length` (default 1250), using 20 amino-acid channels plus an unknown channel.

For a short trial, use a small input table retaining both training and validation rows, reduce `--max-epochs`, and set `--num-workers 0` if needed. `--limit-rows` takes the first rows of the TSV and can accidentally omit validation. CPU execution is supported for contact models with `--accelerator cpu --devices 1`; GPU memory requirements grow with the padded residue-pair map size, so reduce batch size for long proteins.

### Model selection and outputs

By default, adapted D-SCRIPT and QUERY-PATCH select checkpoints by **maximum validation interaction AUPRC** (`val_interaction_auprc`); the baseline selects by minimum validation loss (`val_loss`). The secondary metric is only reported from the selected epoch; it does not alter checkpoint selection. Early stopping defaults to 8 validation epochs without an improvement of at least 0.001.

For contact-focused checkpoint selection, use `--selection-monitor val_contact_auprc` together with `--selection-tiebreaker val_interaction_auprc`. Contact AUPRC is calculated separately for each structural map and averaged with equal weight per map. Training logs use `val_contact_auprc` and `test_contact_auprc`.

Training evaluates the best validation checkpoint on the test split and exports predictions. Use `--skip-test` during hyperparameter selection to suppress test evaluation and exports.

```text
runs/<output-subdir>/
  checkpoints/                      Best checkpoint and last.ckpt
  metrics/version_*/metrics.csv     Epoch-level training/validation/test metrics
  metrics/version_*/hparams.yaml    Logged model hyperparameters
  selected_validation_metrics.tsv  Selected epoch and validation scores
  loss_curve.png
  interaction_auprc_history.png
  contact_auprc_history.png         Contact models only
  test_predictions.tsv             Pair IDs, labels, probabilities, predictions
  contact_maps/                    Example true/predicted test contact-maps
```

Plots and test exports depend on the relevant data being present. Interaction predictions use probability ≥ 0.5; thresholded contact metrics use `--contact-threshold` (default 0.5). Metrics named `auprc` report average precision.

**Final contact results report AUPRC, precision, recall, and F1.** Each metric is calculated separately for every evaluable structural map and then averaged with equal weight per map. Unknown and padded cells are excluded, and synthetic negative maps do not enter these results. The final CSV columns are `test_auprc`, `test_precision`, `test_recall`, and `test_f1`; `validation_f1` records the score used to choose the contact cutoff. Aggregate CSVs report the mean and sample standard deviation across training seeds.

Training contact accuracy, precision, and recall also average the individual map scores equally. When `--negative-contact-maps` is enabled, those thresholded training metrics include synthetic negative maps; contact AUPRC and the final structural contact reports exclude them. Threshold-analysis summaries use `contact_auprc` and `accuracy`, `precision`, `recall`, and `f1` fields for each decision rule. Use [report_contact_thresholds.py](scripts/analysis/report_contact_thresholds.py) to generate the final result tables.

## Threshold selection and analysis

The threshold scripts evaluate **frozen adapted D-SCRIPT or QUERY-PATCH checkpoints**. They select a cutoff on validation data before applying it to test data:

| Script | Purpose |
| --- | --- |
| `scripts/training_pipeline/select_interaction_threshold.py` | Select an interaction threshold by validation accuracy. |
| `scripts/training_pipeline/select_contact_threshold.py` | Select a contact threshold by validation F1; cache known-cell scores and report contact metrics. |
| `scripts/analysis/report_seed_thresholds.py` | Verify and summarize the saved final interaction-threshold experiments across seeds. |
| `scripts/analysis/report_contact_thresholds.py` | Verify final contact-threshold experiments and export per-run metrics and summaries across seeds. |
| `scripts/analysis/report_baseline.py` | Summarize saved baseline results. |
| `scripts/analysis/summarize_dataset_splits.py` | Summarize interaction counts and sequence degrees by split. |

Contact threshold selection maximizes validation F1, averaged equally across contact-maps. The objective is named `f1` (`--selection-metric f1`). The final contact report expects adapted D-SCRIPT and QP-Conv-IC runs with seeds 42, 123, and 456 for each included embedding, evaluated on the same held-out cohort. These report scripts are tailored to the saved experiment comparisons.

For example, select a contact threshold for a completed adapted D-SCRIPT run:

```bash
python scripts/training_pipeline/select_contact_threshold.py \
  --model dscript \
  --checkpoint /path/to/best-checkpoint.ckpt \
  --interactions "$DATA_ROOT/final_datasets/balanced_interactions_strict_cd_hit.tsv" \
  --embedding-manifest "$DATA_ROOT/embeddings/esm2_t33_650M/manifest.tsv" \
  --contact-map-dir "$DATA_ROOT/contact_maps/data" \
  --test-predictions ./runs/dscript_esm2/test_predictions.tsv \
  --batch-size 8 \
  --device cuda \
  --output-dir ./runs/dscript_esm2_contact_threshold
```

Replace the checkpoint placeholder with the selected checkpoint, retain the original evaluation batch size, and choose a new output directory. Padding can affect model scores, so batch size matters for reproducing saved predictions. For historical checkpoints, use `--source-root` to load the matching archived implementation. These analyses require test predictions from a completed evaluation; runs made with `--skip-test` do not produce that file. Each analysis script exposes its arguments through `--help`.

## Further documentation

- [Dataset creation pipeline](docs/dataset_creation_pipeline.md): detailed preparation rationale and historical workflow.
- [Training metrics](docs/training_metrics.md): losses, masking conventions, and metric definitions.
- [QUERY-PATCH history](docs/query_patch_history.md): architecture changes and experimental stages.

The dataset guide predates the move from `data/` to `scripts/data_pipeline/`, and the QUERY-PATCH history records earlier experiment settings. Historical artifact links may refer to files available only in the original workspace. The training metrics guide documents the current metric names and averaging conventions; use the current source and CLI help for accepted argument names.
