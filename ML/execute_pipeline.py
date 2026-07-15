#!/usr/bin/env python
"""Train protein interaction models with PyTorch Lightning.

The default model is a D-SCRIPT-style architecture trained jointly on pair-level
interaction labels and available residue contact maps.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset

import lightning.pytorch as pl
from lightning.pytorch.callbacks import ModelCheckpoint
from lightning.pytorch.loggers import CSVLogger

from models import DScriptInteractionModel, FullyConnectedInteractionModel


# =============================================================================
# General settings, CLI, and shared helpers
# =============================================================================

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AMINO_ACID_TO_INDEX = {amino_acid: index for index, amino_acid in enumerate(AMINO_ACIDS)}
UNKNOWN_AMINO_ACID_INDEX = len(AMINO_ACIDS)
ONE_HOT_DEPTH = len(AMINO_ACIDS) + 1

# Dimer types are kept as small integer tensors so test-time metrics can be
# split by homodimer/heterodimer without passing strings through the DataLoader.
DIMER_TYPE_TO_INDEX = {"homo": 0, "hetero": 1}
INTERACTION_THRESHOLD = 0.5

# These defaults point at the current prepared dimer dataset. They can all be
# overridden from the command line, which is useful for small debug subsets or
# alternative embedding/contact-map directories.
DEFAULT_STRICT_FILE = Path("/nfs/scratch/pdb_dimers/balanced_interactions_strict.tsv")
DEFAULT_SEQUENCE_FILE = Path("/nfs/scratch/pdb_dimers/entity_sequences.tsv")
DEFAULT_EMBEDDING_MANIFEST = Path("/nfs/scratch/pdb_dimers/embeddings/esm2_t33_650M/manifest.tsv")
DEFAULT_CONTACT_MAP_DIR = Path("/nfs/scratch/pdb_dimers/contact_maps/data")
CONTACT_MAP_SAMPLE_COUNT = 5


def parse_args() -> argparse.Namespace:
    """Collect all command-line settings for data, model, and training."""
    parser = argparse.ArgumentParser(description="Train a PDB dimer interaction model.")
    parser.add_argument(
        "--model",
        choices=("dscript", "fully_connected"),
        help="Model family to train.",
    )
    parser.add_argument(
        "--interactions",
        type=Path,
        default=None,
        help="Interaction TSV. Defaults to the strict TSV in /nfs/scratch/pdb_dimers.",
    )
    parser.add_argument(
        "--sequence-file",
        type=Path,
        default=DEFAULT_SEQUENCE_FILE,
        help="Entity sequence TSV for the fully_connected baseline.",
    )
    parser.add_argument(
        "--embedding-manifest",
        type=Path,
        default=DEFAULT_EMBEDDING_MANIFEST,
        help="Embedding manifest with per_residue_path entries for D-SCRIPT.",
    )
    parser.add_argument(
        "--contact-map-dir",
        type=Path,
        default=DEFAULT_CONTACT_MAP_DIR,
        help="Directory containing <pdb_id>-assembly<assembly_number>.npy contact maps.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("ML/runs/dscript_interactions"),
        help="Where Lightning logs and checkpoints are written.",
    )
    parser.add_argument("--max-length", type=int, default=1250, help="One-hot baseline pad/truncate length.")
    parser.add_argument("--batch-size", type=int, default=16, help="Pairs per training batch.")
    parser.add_argument("--num-workers", type=int, default=4, help="DataLoader worker processes.")
    parser.add_argument("--max-epochs", type=int, default=10, help="Number of training epochs.")
    parser.add_argument("--learning-rate", type=float, default=1e-3, help="Adam learning rate.")
    parser.add_argument("--hidden-dim", type=int, default=512, help="Fully connected baseline hidden width.")
    parser.add_argument("--dropout", type=float, default=0.2, help="Fully connected baseline dropout.")
    parser.add_argument("--embedding-dim", type=int, default=1280, help="D-SCRIPT per-residue embedding size.")
    parser.add_argument("--projection-dim", type=int, default=100, help="D-SCRIPT projected residue size.")
    parser.add_argument("--contact-hidden-dim", type=int, default=50, help="D-SCRIPT contact module hidden channels.")
    parser.add_argument("--contact-width", type=int, default=7, help="Odd convolution width for contact prediction.")
    parser.add_argument("--projection-dropout", type=float, default=0.5, help="D-SCRIPT projection dropout.")
    parser.add_argument(
        "--interaction-loss-lambda",
        type=float,
        default=0.5,
        help="Weight for interaction BCE. Contact BCE weight is 1 - this value.",
    )
    parser.add_argument(
        "--contact-threshold",
        type=float,
        default=0.5,
        help="Threshold at which predicted contact probabilities count as contacts for metrics.",
    )
    parser.add_argument(
        "--loss-mode",
        choices=("auto", "basic", "dimer_type"),
        default="auto",
        help="'auto' uses dimer-type interaction weighting for keep_homomers and plain BCE otherwise.",
    )
    parser.add_argument("--accelerator", default="auto", help="Lightning accelerator, e.g. auto, cpu, gpu.")
    parser.add_argument("--devices", default="auto", help="Lightning devices, e.g. auto or 1.")
    parser.add_argument("--log-every-n-steps", type=int, default=50, help="How often Lightning logs progress.")
    parser.add_argument(
        "--limit-rows",
        type=int,
        default=None,
        help="Optional small-row limit for quick debugging. Leave unset for real training.",
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    """Fail early on invalid hyperparameters."""
    # ``interaction_loss_lambda`` interpolates between the two supervised tasks:
    # 1.0 means only pair-level interaction loss, 0.0 means only contact-map loss
    # whenever contact labels exist in the batch.
    if not 0.0 <= args.interaction_loss_lambda <= 1.0:
        raise ValueError("--interaction-loss-lambda must be between 0 and 1")

    # The threshold is used for metrics only. The contact-map loss uses the raw
    # logits directly because thresholding would make the loss non-differentiable.
    if not 0.0 <= args.contact_threshold <= 1.0:
        raise ValueError("--contact-threshold must be between 0 and 1")

    # The contact convolution is padded symmetrically. Requiring an odd kernel
    # keeps each output cell centered on the same residue pair.
    if args.contact_width <= 0 or args.contact_width % 2 == 0:
        raise ValueError("--contact-width must be a positive odd integer")


def split_pair(value: object, column_name: str) -> tuple[str, str]:
    """Split a comma-separated pair column and return stripped values."""
    # Both entity_pair and new_cluster_pair are stored as "left,right" strings.
    # Centralizing the parsing keeps the baseline and D-SCRIPT data paths
    # consistent and gives a useful error for malformed rows.
    pieces = str(value).split(",", maxsplit=1)
    if len(pieces) != 2 or not pieces[0].strip() or not pieces[1].strip():
        raise ValueError(f"Could not split {column_name} value into two entries: {value!r}")
    return pieces[0].strip(), pieces[1].strip()


def add_interaction_weights(examples: pd.DataFrame, loss_mode: str) -> pd.DataFrame:
    """Add per-example interaction BCE weights."""
    # Interaction weights affect only the pair-level BCE. They do not affect
    # contact-map supervision, because contact maps have their own per-cell BCE.
    train_df = examples[examples["split"] == "train"]
    if train_df.empty:
        raise ValueError("No training rows found. Expected split column to contain 'train'.")

    weight_map: dict[tuple[str, int], float] = {}
    if loss_mode != "basic":
        # In dimer_type mode, each (dimer_type, label) cell receives inverse
        # frequency weighting computed from the training split only.
        counts = train_df.groupby(["dimer_type", "label"]).size()
        total_examples = float(counts.sum())
        number_of_cells = float(len(counts))
        weight_map = {
            (dimer_type, int(label)): total_examples / (number_of_cells * float(count))
            for (dimer_type, label), count in counts.items()
        }

    examples = examples.copy()
    examples["loss_weight"] = [
        weight_map.get((row.dimer_type, int(row.label)), 1.0) for row in examples.itertuples(index=False)
    ]
    examples["loss_weight"] = examples["loss_weight"].astype("float32")
    return examples


def log_test_accuracy_by_dimer_type(
    module: pl.LightningModule,
    dimer_type_index: torch.Tensor,
    predictions: torch.Tensor,
    labels: torch.Tensor,
) -> None:
    """Log separate test accuracies for homo- and heterodimers."""
    for dimer_type_name, dimer_type_value in DIMER_TYPE_TO_INDEX.items():
        dimer_type_mask = dimer_type_index == dimer_type_value
        dimer_type_count = int(dimer_type_mask.sum().item())

        if dimer_type_count == 0:
            continue

        dimer_type_accuracy = (predictions[dimer_type_mask] == labels[dimer_type_mask]).float().mean()
        module.log(
            f"test_accuracy_{dimer_type_name}",
            dimer_type_accuracy,
            prog_bar=True,
            on_step=False,
            on_epoch=True,
            batch_size=dimer_type_count,
        )


def prepare_common_columns(examples: pd.DataFrame) -> pd.DataFrame:
    """Normalize label, split, and dimer-type columns."""
    # All model paths expect labels as numeric 0/1, split names in the same
    # train/val/test vocabulary, and lower-case dimer type names.
    examples = examples.copy()
    examples["label"] = pd.to_numeric(examples["label"], errors="raise").astype("int64")
    examples["split"] = examples["split"].astype(str).str.lower().replace({"validation": "val", "valid": "val"})
    examples["dimer_type"] = examples["dimer_type"].fillna("unknown").astype(str).str.lower()
    return examples


def print_split_summary(
    train_dataset: Dataset,
    val_dataset: Dataset,
    test_dataset: Dataset,
    loss_mode: str,
) -> None:
    """Print dataset sizes after setup."""
    print(
        "Prepared examples: "
        f"train={len(train_dataset)}, val={len(val_dataset)}, test={len(test_dataset)}, "
        f"loss_mode={loss_mode}"
    )


def save_loss_plot(metrics_path: Path, output_path: Path) -> None:
    """Plot epoch-level train/validation losses from Lightning's CSV metrics."""
    if not metrics_path.exists():
        print(f"Skipping loss plot because metrics file does not exist: {metrics_path}")
        return

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        print(f"Skipping loss plot because matplotlib is not installed: {exc}")
        return

    metrics = pd.read_csv(metrics_path)
    if "epoch" not in metrics.columns:
        print(f"Skipping loss plot because {metrics_path} has no epoch column")
        return

    # Lightning writes one sparse CSV row per logged metric. Grouping by epoch
    # reconstructs the train/validation loss curves that were logged with
    # on_epoch=True. If no validation split exists, the plot contains only train.
    loss_series: dict[str, pd.Series] = {}
    for metric_name, display_name in (("train_loss", "Train loss"), ("val_loss", "Validation loss")):
        if metric_name not in metrics.columns:
            continue

        metric_rows = metrics[["epoch", metric_name]].dropna()
        if metric_rows.empty:
            continue

        metric_rows["epoch"] = metric_rows["epoch"].astype(int)
        loss_series[display_name] = metric_rows.groupby("epoch")[metric_name].last().sort_index()

    if not loss_series:
        print(f"Skipping loss plot because no train_loss or val_loss rows were found in {metrics_path}")
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(8, 5))
    for display_name, values in loss_series.items():
        # Epochs are shown as 1-based labels because the saved checkpoints and
        # console progress are easier to compare to human-readable epoch counts.
        plt.plot(values.index + 1, values.values, marker="o", label=display_name)

    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.title("Training and Validation Loss")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=160)
    plt.close()
    print(f"Saved loss plot to {output_path}")


def move_tensor_batch_to_device(batch: dict[str, torch.Tensor], device: torch.device) -> dict[str, torch.Tensor]:
    """Move one DataLoader batch to the device that currently holds the model."""
    return {key: value.to(device) if isinstance(value, torch.Tensor) else value for key, value in batch.items()}


def safe_filename_part(value: object) -> str:
    """Convert row identifiers into a filesystem-safe filename fragment."""
    text = str(value)
    safe_text = "".join(character if character.isalnum() or character in ("-", "_", ".") else "_" for character in text)
    return safe_text.strip("_") or "unknown"


# =============================================================================
# Fully connected baseline
# =============================================================================


def one_hot_encode_sequence(sequence: str, max_length: int) -> torch.Tensor:
    """Convert one protein sequence to a padded one-hot tensor."""
    # The baseline has a fixed input size. Residues beyond max_length are
    # truncated and shorter sequences keep all-zero padding rows at the end.
    encoded = torch.zeros((max_length, ONE_HOT_DEPTH), dtype=torch.float32)
    trimmed_sequence = sequence.upper()[:max_length]

    for position, amino_acid in enumerate(trimmed_sequence):
        # Unknown or uncommon residue symbols share the final channel instead of
        # causing a failure during data loading.
        channel = AMINO_ACID_TO_INDEX.get(amino_acid, UNKNOWN_AMINO_ACID_INDEX)
        encoded[position, channel] = 1.0

    return encoded


class FullyConnectedDataset(Dataset):
    """Dataset that returns one flattened one-hot vector per protein pair."""

    def __init__(self, examples: pd.DataFrame, max_length: int) -> None:
        self.examples = examples.reset_index(drop=True)
        self.max_length = max_length

    def __len__(self) -> int:
        """Return the number of protein pairs in this split."""
        return len(self.examples)

    def __getitem__(self, index: int) -> dict[str, torch.Tensor]:
        """Encode both sequences and return tensors used by the Lightning module."""
        row = self.examples.iloc[index]

        # Each sequence becomes [max_length, ONE_HOT_DEPTH], then flattening and
        # concatenation produce one fixed-length vector for the protein pair.
        sequence_1 = one_hot_encode_sequence(row.sequence_1, self.max_length).flatten()
        sequence_2 = one_hot_encode_sequence(row.sequence_2, self.max_length).flatten()
        features = torch.cat((sequence_1, sequence_2), dim=0)

        return {
            "features": features,
            "label": torch.tensor(float(row.label), dtype=torch.float32),
            "loss_weight": torch.tensor(float(row.loss_weight), dtype=torch.float32),
            "dimer_type_index": torch.tensor(DIMER_TYPE_TO_INDEX.get(row.dimer_type, -1), dtype=torch.long),
        }


class FullyConnectedDataModule(pl.LightningDataModule):
    """LightningDataModule for the one-hot fully connected baseline."""

    def __init__(
        self,
        interaction_path: Path,
        sequence_path: Path,
        max_length: int,
        batch_size: int,
        num_workers: int,
        loss_mode: str,
        limit_rows: int | None,
    ) -> None:
        super().__init__()
        self.interaction_path = interaction_path
        self.sequence_path = sequence_path
        self.max_length = max_length
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.loss_mode = loss_mode
        self.limit_rows = limit_rows

        self.train_dataset: FullyConnectedDataset | None = None
        self.val_dataset: FullyConnectedDataset | None = None
        self.test_dataset: FullyConnectedDataset | None = None

    def setup(self, stage: str | None = None) -> None:
        """Load data, connect pairs to sequences, and split into datasets."""
        # Lightning may call setup more than once for different stages. The
        # guard avoids re-reading large TSV files once the datasets exist.
        if self.train_dataset is not None:
            return

        # The baseline needs sequence strings for both entities in each row.
        # D-SCRIPT uses embeddings instead and therefore follows a separate path.
        needed_columns = {"entity_pair", "split", "label", "dimer_type"}
        interaction_df = pd.read_csv(self.interaction_path, sep="\t", nrows=self.limit_rows)
        missing_columns = needed_columns.difference(interaction_df.columns)
        if missing_columns:
            raise ValueError(f"{self.interaction_path} is missing columns: {sorted(missing_columns)}")

        sequence_df = pd.read_csv(self.sequence_path, sep="\t", usecols=["entity_name", "sequence"])
        sequence_df = sequence_df.dropna(subset=["entity_name", "sequence"])
        sequence_df = sequence_df.drop_duplicates(subset=["entity_name"], keep="first")
        sequence_by_entity = dict(zip(sequence_df["entity_name"].astype(str), sequence_df["sequence"].astype(str)))

        # Convert the stored "entity_1,entity_2" string into two lookup keys and
        # attach the corresponding amino-acid sequences.
        examples = interaction_df[["entity_pair", "split", "label", "dimer_type"]].copy()
        examples[["entity_1", "entity_2"]] = examples["entity_pair"].apply(
            lambda value: pd.Series(split_pair(value, "entity_pair"))
        )
        examples["sequence_1"] = examples["entity_1"].map(sequence_by_entity)
        examples["sequence_2"] = examples["entity_2"].map(sequence_by_entity)

        missing_entity_1 = set(examples.loc[examples["sequence_1"].isna(), "entity_1"])
        missing_entity_2 = set(examples.loc[examples["sequence_2"].isna(), "entity_2"])
        missing_entities = sorted(missing_entity_1.union(missing_entity_2))
        if missing_entities:
            preview = ", ".join(missing_entities[:10])
            raise ValueError(f"Missing sequences for {len(missing_entities)} entities. First examples: {preview}")

        # From here on, the baseline and D-SCRIPT paths share the same label,
        # split, dimer-type, and interaction-weight conventions.
        examples = prepare_common_columns(examples)
        examples = add_interaction_weights(examples, self.loss_mode)

        self.train_dataset = FullyConnectedDataset(examples[examples["split"] == "train"], self.max_length)
        self.val_dataset = FullyConnectedDataset(examples[examples["split"] == "val"], self.max_length)
        self.test_dataset = FullyConnectedDataset(examples[examples["split"] == "test"], self.max_length)
        print_split_summary(self.train_dataset, self.val_dataset, self.test_dataset, self.loss_mode)

    def train_dataloader(self) -> DataLoader:
        """Return shuffled training batches."""
        return DataLoader(self.train_dataset, batch_size=self.batch_size, shuffle=True, num_workers=self.num_workers)

    def val_dataloader(self) -> DataLoader | None:
        """Return validation batches when a validation split exists."""
        if self.val_dataset is None or len(self.val_dataset) == 0:
            return None
        return DataLoader(self.val_dataset, batch_size=self.batch_size, shuffle=False, num_workers=self.num_workers)

    def test_dataloader(self) -> DataLoader | None:
        """Return test batches when a test split exists."""
        if self.test_dataset is None or len(self.test_dataset) == 0:
            return None
        return DataLoader(self.test_dataset, batch_size=self.batch_size, shuffle=False, num_workers=self.num_workers)


class FullyConnectedLightningModule(pl.LightningModule):
    """Lightning wrapper around the fully connected interaction model."""

    def __init__(self, input_dim: int, hidden_dim: int, dropout: float, learning_rate: float) -> None:
        super().__init__()
        self.save_hyperparameters()
        self.model = FullyConnectedInteractionModel(input_dim=input_dim, hidden_dim=hidden_dim, dropout=dropout)
        self.learning_rate = learning_rate

    def forward(self, features: torch.Tensor) -> torch.Tensor:
        """Run a batch through the model and return logits."""
        return self.model(features)

    def shared_step(self, batch: dict[str, torch.Tensor], stage: str) -> torch.Tensor:
        """Compute weighted BCE loss and simple unweighted accuracy for one batch."""
        labels = batch["label"]
        weights = batch["loss_weight"]
        logits = self(batch["features"])

        # The baseline has only one supervised objective: predict whether the
        # protein pair interacts. Any dimer-type balancing is applied here.
        per_example_loss = F.binary_cross_entropy_with_logits(logits, labels, reduction="none")
        loss = (per_example_loss * weights).mean()

        # Accuracy is deliberately unweighted so it remains easy to interpret.
        predictions = (torch.sigmoid(logits) >= INTERACTION_THRESHOLD).float()
        accuracy = (predictions == labels).float().mean()

        self.log(f"{stage}_loss", loss, prog_bar=True, on_step=False, on_epoch=True, batch_size=labels.size(0))
        self.log(f"{stage}_accuracy", accuracy, prog_bar=True, on_step=False, on_epoch=True, batch_size=labels.size(0))

        if stage == "test":
            log_test_accuracy_by_dimer_type(self, batch["dimer_type_index"], predictions, labels)

        return loss

    def training_step(self, batch: dict[str, torch.Tensor], _batch_idx: int) -> torch.Tensor:
        """Train on one batch."""
        return self.shared_step(batch, "train")

    def validation_step(self, batch: dict[str, torch.Tensor], _batch_idx: int) -> torch.Tensor:
        """Validate on one batch."""
        return self.shared_step(batch, "val")

    def test_step(self, batch: dict[str, torch.Tensor], _batch_idx: int) -> torch.Tensor:
        """Evaluate on one test batch."""
        return self.shared_step(batch, "test")

    def configure_optimizers(self) -> torch.optim.Optimizer:
        """Use Adam as a simple default optimizer."""
        return torch.optim.Adam(self.parameters(), lr=self.learning_rate)


# =============================================================================
# D-SCRIPT model with real contact-map supervision
# =============================================================================


class DScriptDataset(Dataset):
    """Dataset that loads per-residue embeddings and optional true contact maps."""

    def __init__(self, examples: pd.DataFrame) -> None:
        self.examples = examples.reset_index(drop=True)

    def __len__(self) -> int:
        """Return the number of protein pairs in this split."""
        return len(self.examples)

    def __getitem__(self, index: int) -> dict[str, object]:
        """Load embedding tensors and a true contact map when available."""
        row = self.examples.iloc[index]

        # Embeddings are stored as one tensor per cluster/entity:
        # [residues, embedding_dim]. They are kept variable-length here and are
        # padded later by collate_dscript_batch.
        embedding_1 = self.load_embedding(Path(row.embedding_path_1))

        # Homomer pairs can point to the same embedding file for both proteins.
        # Reusing the already loaded tensor avoids one redundant disk read.
        if row.embedding_path_1 == row.embedding_path_2:
            embedding_2 = embedding_1
        else:
            embedding_2 = self.load_embedding(Path(row.embedding_path_2))

        # Positive examples should have a real contact map. Negative examples
        # usually do not; None is allowed and becomes an all -1 padded target.
        # When a map exists, its -1 cells later mask unresolved residue pairs out
        # of both contact supervision and interaction-score pooling.
        contact_map = self.load_contact_map(row.contact_map_path, embedding_1.shape[0], embedding_2.shape[0])

        return {
            "embedding_1": embedding_1,
            "embedding_2": embedding_2,
            "contact_map": contact_map,
            "label": torch.tensor(float(row.label), dtype=torch.float32),
            "loss_weight": torch.tensor(float(row.loss_weight), dtype=torch.float32),
            "dimer_type_index": torch.tensor(DIMER_TYPE_TO_INDEX.get(row.dimer_type, -1), dtype=torch.long),
        }

    @staticmethod
    def load_embedding(path: Path) -> torch.Tensor:
        """Load one per-residue embedding tensor from disk."""
        tensor = torch.load(path, map_location="cpu")
        if not isinstance(tensor, torch.Tensor) or tensor.ndim != 2:
            raise ValueError(f"Expected a [residues, embedding_dim] tensor at {path}")

        # The model and losses expect floating-point tensors. Loading on CPU
        # keeps workers independent of the training device; Lightning moves the
        # final batch to GPU/CPU as configured.
        return tensor.to(dtype=torch.float32)

    @staticmethod
    def load_contact_map(path_value: object, residues_1: int, residues_2: int) -> torch.Tensor | None:
        """Load one true contact map, preserving -1 as unknown."""
        if pd.isna(path_value) or not str(path_value):
            return None

        path = Path(str(path_value))
        contact_map = np.load(path)

        # This check is important: the supervised contact loss assumes that
        # contact_map[i, j] describes residue i from protein 1 and residue j from
        # protein 2. A mismatch would silently train on shifted/wrong labels.
        if contact_map.shape != (residues_1, residues_2):
            raise ValueError(
                f"Contact map shape mismatch for {path}: got {contact_map.shape}, "
                f"expected {(residues_1, residues_2)}"
            )

        # Values are 0 = known non-contact, 1 = known contact, -1 = unknown.
        # Keeping -1 intact lets the loss mask out unknown cells later.
        return torch.as_tensor(contact_map.astype(np.float32, copy=False))


def collate_dscript_batch(items: list[dict[str, object]]) -> dict[str, torch.Tensor]:
    """Pad variable-length embeddings and contact maps within one D-SCRIPT batch."""
    batch_size = len(items)
    embedding_dim = items[0]["embedding_1"].shape[1]

    # Each batch is padded only to the largest protein lengths in that batch.
    # This avoids a global max-length tensor and keeps contact maps much smaller.
    max_residues_1 = max(item["embedding_1"].shape[0] for item in items)
    max_residues_2 = max(item["embedding_2"].shape[0] for item in items)

    # Padded embedding rows are zeros. The boolean masks below tell the model
    # which residue positions are real and which are just padding.
    embeddings_1 = torch.zeros((batch_size, max_residues_1, embedding_dim), dtype=torch.float32)
    embeddings_2 = torch.zeros((batch_size, max_residues_2, embedding_dim), dtype=torch.float32)
    mask_1 = torch.zeros((batch_size, max_residues_1), dtype=torch.bool)
    mask_2 = torch.zeros((batch_size, max_residues_2), dtype=torch.bool)

    # Contact-map padding is -1, the same code used by the real data for
    # "unknown". This means padded cells and genuinely unknown cells are both
    # ignored by the contact-map loss.
    contact_maps = torch.full((batch_size, max_residues_1, max_residues_2), -1.0, dtype=torch.float32)
    interaction_pair_masks = torch.zeros((batch_size, max_residues_1, max_residues_2), dtype=torch.bool)

    for batch_index, item in enumerate(items):
        embedding_1 = item["embedding_1"]
        embedding_2 = item["embedding_2"]
        residues_1 = embedding_1.shape[0]
        residues_2 = embedding_2.shape[0]

        if embedding_1.shape[1] != embedding_dim or embedding_2.shape[1] != embedding_dim:
            raise ValueError("All embeddings in a D-SCRIPT batch must have the same embedding dimension")

        # Copy each variable-length tensor into the top-left valid region of the
        # padded batch tensor, then mark the copied residue positions as real.
        embeddings_1[batch_index, :residues_1] = embedding_1
        embeddings_2[batch_index, :residues_2] = embedding_2
        mask_1[batch_index, :residues_1] = True
        mask_2[batch_index, :residues_2] = True

        contact_map = item["contact_map"]
        if contact_map is not None:
            # Real contact maps are copied into the same top-left valid region.
            # Remaining cells stay -1 and are ignored by the loss. For examples
            # with structural supervision, interaction pooling uses only known
            # residue pairs, so unresolved -1 cells cannot become unconstrained
            # evidence for the pair-level interaction prediction.
            contact_maps[batch_index, :residues_1, :residues_2] = contact_map
            interaction_pair_masks[batch_index, :residues_1, :residues_2] = contact_map >= 0
        else:
            # Examples without any contact map, usually synthetic/non-structural
            # negatives, do not identify unresolved residues. Keep the previous
            # behavior and allow all real residue pairs to contribute to the
            # interaction objective for those examples.
            interaction_pair_masks[batch_index, :residues_1, :residues_2] = True

    return {
        "embedding_1": embeddings_1,
        "embedding_2": embeddings_2,
        "mask_1": mask_1,
        "mask_2": mask_2,
        "contact_map": contact_maps,
        "interaction_pair_mask": interaction_pair_masks,
        "label": torch.stack([item["label"] for item in items]),
        "loss_weight": torch.stack([item["loss_weight"] for item in items]),
        "dimer_type_index": torch.stack([item["dimer_type_index"] for item in items]),
    }


class DScriptDataModule(pl.LightningDataModule):
    """LightningDataModule for D-SCRIPT over precomputed per-residue embeddings."""

    def __init__(
        self,
        interaction_path: Path,
        embedding_manifest: Path,
        contact_map_dir: Path,
        batch_size: int,
        num_workers: int,
        loss_mode: str,
        limit_rows: int | None,
    ) -> None:
        super().__init__()
        self.interaction_path = interaction_path
        self.embedding_manifest = embedding_manifest
        self.contact_map_dir = contact_map_dir
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.loss_mode = loss_mode
        self.limit_rows = limit_rows

        self.train_dataset: DScriptDataset | None = None
        self.val_dataset: DScriptDataset | None = None
        self.test_dataset: DScriptDataset | None = None

    def setup(self, stage: str | None = None) -> None:
        """Load interaction rows and connect them to embedding/contact-map files."""
        # As in the baseline DataModule, setup can be called multiple times by
        # Lightning. Once datasets exist, no additional disk scan is needed.
        if self.train_dataset is not None:
            return

        # D-SCRIPT does not need raw sequences here. It needs cluster ids for
        # embedding lookup and PDB/assembly ids for contact-map lookup.
        needed_columns = {"new_cluster_pair", "split", "label", "dimer_type", "pdb_id", "assembly_number"}
        interaction_df = pd.read_csv(self.interaction_path, sep="\t", nrows=self.limit_rows)
        missing_columns = needed_columns.difference(interaction_df.columns)
        if missing_columns:
            raise ValueError(f"{self.interaction_path} is missing columns: {sorted(missing_columns)}")

        embedding_df = pd.read_csv(self.embedding_manifest, sep="\t")
        embedding_by_cluster = self.build_embedding_lookup(embedding_df)

        # Split the stored "cluster_1,cluster_2" key and attach one per-residue
        # embedding path to each side of the pair.
        examples = interaction_df[
            ["new_cluster_pair", "split", "label", "dimer_type", "pdb_id", "assembly_number"]
        ].copy()
        examples[["cluster_1", "cluster_2"]] = examples["new_cluster_pair"].apply(
            lambda value: pd.Series(split_pair(value, "new_cluster_pair"))
        )
        examples["embedding_path_1"] = examples["cluster_1"].map(embedding_by_cluster)
        examples["embedding_path_2"] = examples["cluster_2"].map(embedding_by_cluster)

        missing_cluster_1 = set(examples.loc[examples["embedding_path_1"].isna(), "cluster_1"])
        missing_cluster_2 = set(examples.loc[examples["embedding_path_2"].isna(), "cluster_2"])
        missing_clusters = sorted(missing_cluster_1.union(missing_cluster_2))
        if missing_clusters:
            preview = ", ".join(missing_clusters[:10])
            raise ValueError(f"Missing embeddings for {len(missing_clusters)} clusters. First examples: {preview}")

        # Contact maps exist for real PDB assemblies. Negative examples can be
        # synthetic/non-structural pairs and therefore may not have pdb_id or a
        # contact map. Those rows simply receive no contact supervision.
        examples["contact_map_path"] = [
            self.resolve_contact_map_path(row.pdb_id, row.assembly_number) for row in examples.itertuples(index=False)
        ]
        examples = prepare_common_columns(examples)

        # Positive structural examples are expected to have contact labels. If a
        # positive map is missing, failing here is better than silently training
        # the contact objective on an incomplete positive set.
        missing_positive_maps = examples[(examples["label"] == 1) & examples["contact_map_path"].isna()]
        if not missing_positive_maps.empty:
            preview = missing_positive_maps[["pdb_id", "assembly_number"]].head(10).to_dict("records")
            raise ValueError(f"Missing contact maps for {len(missing_positive_maps)} positive examples: {preview}")

        # The interaction objective still uses the same class/dimer-type weights
        # as the baseline. Contact-map loss is handled inside the LightningModule.
        examples = add_interaction_weights(examples, self.loss_mode)

        self.train_dataset = DScriptDataset(examples[examples["split"] == "train"])
        self.val_dataset = DScriptDataset(examples[examples["split"] == "val"])
        self.test_dataset = DScriptDataset(examples[examples["split"] == "test"])

        print_split_summary(self.train_dataset, self.val_dataset, self.test_dataset, self.loss_mode)
        print(f"Contact-supervised examples: {examples['contact_map_path'].notna().sum()}/{len(examples)}")

    @staticmethod
    def build_embedding_lookup(embedding_df: pd.DataFrame) -> dict[str, str]:
        """Return a cluster-id to per-residue embedding path mapping."""
        needed_columns = {"new_cluster_id", "per_residue_path"}
        missing_columns = needed_columns.difference(embedding_df.columns)
        if missing_columns:
            raise ValueError(f"Embedding manifest is missing columns: {sorted(missing_columns)}")

        # The manifest may contain additional metadata columns. Only the stable
        # cluster id and per-residue tensor path are needed for training.
        embedding_df = embedding_df.dropna(subset=["new_cluster_id", "per_residue_path"])
        return dict(zip(embedding_df["new_cluster_id"].astype(str), embedding_df["per_residue_path"].astype(str)))

    def resolve_contact_map_path(self, pdb_id: object, assembly_number: object) -> str | float:
        """Resolve a row's real contact map path, or NaN when none exists."""
        if pd.isna(pdb_id) or pd.isna(assembly_number):
            return np.nan

        # Contact maps are saved as e.g. "10bl-assembly1.npy". Lower-casing the
        # PDB id makes the lookup independent of TSV capitalization.
        path = self.contact_map_dir / f"{str(pdb_id).lower()}-assembly{int(assembly_number)}.npy"
        return str(path) if path.exists() else np.nan

    def train_dataloader(self) -> DataLoader:
        """Return shuffled D-SCRIPT training batches."""
        return DataLoader(
            self.train_dataset,
            batch_size=self.batch_size,
            shuffle=True,
            num_workers=self.num_workers,
            collate_fn=collate_dscript_batch,
        )

    def val_dataloader(self) -> DataLoader | None:
        """Return validation batches when a validation split exists."""
        if self.val_dataset is None or len(self.val_dataset) == 0:
            return None
        return DataLoader(
            self.val_dataset,
            batch_size=self.batch_size,
            shuffle=False,
            num_workers=self.num_workers,
            collate_fn=collate_dscript_batch,
        )

    def test_dataloader(self) -> DataLoader | None:
        """Return test batches when a test split exists."""
        if self.test_dataset is None or len(self.test_dataset) == 0:
            return None
        return DataLoader(
            self.test_dataset,
            batch_size=self.batch_size,
            shuffle=False,
            num_workers=self.num_workers,
            collate_fn=collate_dscript_batch,
        )


class DScriptLightningModule(pl.LightningModule):
    """Lightning wrapper for supervised D-SCRIPT interaction and contact-map learning."""

    def __init__(
        self,
        embedding_dim: int,
        projection_dim: int,
        contact_hidden_dim: int,
        contact_width: int,
        projection_dropout: float,
        learning_rate: float,
        interaction_loss_lambda: float,
        contact_threshold: float,
    ) -> None:
        super().__init__()
        self.save_hyperparameters()
        self.model = DScriptInteractionModel(
            embedding_dim=embedding_dim,
            projection_dim=projection_dim,
            contact_hidden_dim=contact_hidden_dim,
            contact_width=contact_width,
            projection_dropout=projection_dropout,
        )
        self.learning_rate = learning_rate
        self.interaction_loss_lambda = interaction_loss_lambda
        self.contact_threshold = contact_threshold

    def forward(
        self,
        embedding_1: torch.Tensor,
        embedding_2: torch.Tensor,
        mask_1: torch.Tensor,
        mask_2: torch.Tensor,
        interaction_pair_mask: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """Run a batch through D-SCRIPT and return interaction logits."""
        return self.model(embedding_1, embedding_2, mask_1, mask_2, interaction_pair_mask=interaction_pair_mask)

    def shared_step(self, batch: dict[str, torch.Tensor], stage: str) -> torch.Tensor:
        """Compute lambda-weighted interaction and contact-map losses."""
        labels = batch["label"]
        weights = batch["loss_weight"]

        # One D-SCRIPT forward pass gives all outputs needed for both tasks:
        # - interaction_logits: one raw logit per protein pair,
        # - contact_probabilities: sigmoid(contact_logits), values in [0, 1],
        # - contact_logits: raw per-cell logits for stable contact BCE.
        interaction_logits, contact_probabilities, contact_logits = self.model(
            batch["embedding_1"],
            batch["embedding_2"],
            batch["mask_1"],
            batch["mask_2"],
            interaction_pair_mask=batch["interaction_pair_mask"],
            return_contact_map=True,
            return_contact_logits=True,
        )

        # Pair-level interaction loss is computed for every row in the batch,
        # including rows without a contact map. The optional weights come from
        # dimer-type/label balancing in add_interaction_weights.
        per_example_interaction_loss = F.binary_cross_entropy_with_logits(
            interaction_logits,
            labels,
            reduction="none",
        )
        interaction_loss = (per_example_interaction_loss * weights).mean()
        contact_loss, contact_metrics = self.compute_contact_loss_and_metrics(
            contact_logits,
            contact_probabilities,
            batch,
        )

        # Lambda controls the tradeoff between the two tasks:
        #   lambda = 1.0 -> only interaction prediction,
        #   lambda = 0.0 -> only contact-map prediction when maps are present.
        # If a batch has no known contact cells, falling back to interaction loss
        # avoids multiplying by a zero contact loss and wasting the batch.
        if contact_metrics["known_contacts"] > 0:
            loss = self.interaction_loss_lambda * interaction_loss + (1.0 - self.interaction_loss_lambda) * contact_loss
        else:
            loss = interaction_loss

        # Interaction classification uses the same probability cutoff as the
        # prediction export, so reported accuracy and saved labels agree.
        interaction_predictions = (torch.sigmoid(interaction_logits) >= INTERACTION_THRESHOLD).float()
        interaction_accuracy = (interaction_predictions == labels).float().mean()

        batch_size = labels.size(0)
        self.log(f"{stage}_loss", loss, prog_bar=True, on_step=False, on_epoch=True, batch_size=batch_size)
        self.log(
            f"{stage}_interaction_loss",
            interaction_loss,
            prog_bar=False,
            on_step=False,
            on_epoch=True,
            batch_size=batch_size,
        )
        self.log(
            f"{stage}_contact_loss",
            contact_loss,
            prog_bar=False,
            on_step=False,
            on_epoch=True,
            batch_size=batch_size,
        )
        self.log(
            f"{stage}_accuracy",
            interaction_accuracy,
            prog_bar=True,
            on_step=False,
            on_epoch=True,
            batch_size=batch_size,
        )

        if contact_metrics["known_contacts"] > 0:
            self.log_contact_metrics(stage, contact_metrics, batch_size)

        if stage == "test":
            log_test_accuracy_by_dimer_type(self, batch["dimer_type_index"], interaction_predictions, labels)

        return loss

    def compute_contact_loss_and_metrics(
        self,
        contact_logits: torch.Tensor,
        contact_probabilities: torch.Tensor,
        batch: dict[str, torch.Tensor],
    ) -> tuple[torch.Tensor, dict[str, torch.Tensor | int]]:
        """Compute contact BCE on known map cells and thresholded contact metrics."""
        true_contact_map = batch["contact_map"]

        # Contact targets use:
        #   1  = known contact,
        #   0  = known non-contact,
        #  -1  = unknown or padding.
        # Only cells with target >= 0 are supervised or counted in metrics.
        known_mask = true_contact_map >= 0
        known_contacts = int(known_mask.sum().item())
        if known_contacts == 0:
            zero = contact_logits.new_tensor(0.0)
            return zero, {"known_contacts": 0}

        # Clamp maps -1 -> 0 before BCE, then multiply by known_mask so those
        # clamped unknown cells contribute exactly zero to the final loss.
        contact_targets = true_contact_map.clamp_min(0.0)
        per_cell_loss = F.binary_cross_entropy_with_logits(
            contact_logits,
            contact_targets,
            reduction="none",
        )
        known_mask_float = known_mask.to(dtype=per_cell_loss.dtype)

        # Average contact BCE within each supervised example first, then average
        # examples. This prevents large proteins/maps from dominating the contact
        # objective merely because they have more residue pairs.
        per_example_known = known_mask_float.flatten(start_dim=1).sum(dim=1)
        per_example_loss = (per_cell_loss * known_mask_float).flatten(start_dim=1).sum(dim=1)
        per_example_loss = per_example_loss / per_example_known.clamp_min(1.0)
        contact_loss = per_example_loss[per_example_known > 0].mean()

        # The threshold is used only for interpretable contact-map metrics. The
        # differentiable loss above still trains on raw logits.
        targets = true_contact_map[known_mask]
        predicted_contacts = contact_probabilities[known_mask] >= self.contact_threshold
        true_contacts = targets.bool()
        contact_accuracy = (predicted_contacts == true_contacts).float().mean()

        # clamp_min(1.0) avoids NaN precision/recall in batches with no predicted
        # positives or no actual positives. In those edge cases the numerator is
        # zero, so the metric becomes zero.
        true_positives = (predicted_contacts & true_contacts).sum().float()
        predicted_positive = predicted_contacts.sum().float()
        actual_positive = true_contacts.sum().float()
        precision = true_positives / predicted_positive.clamp_min(1.0)
        recall = true_positives / actual_positive.clamp_min(1.0)

        return contact_loss, {
            "known_contacts": known_contacts,
            "contact_accuracy": contact_accuracy,
            "contact_precision": precision,
            "contact_recall": recall,
        }

    def log_contact_metrics(
        self,
        stage: str,
        metrics: dict[str, torch.Tensor | int],
        batch_size: int,
    ) -> None:
        """Log contact-map metrics computed with the configured threshold."""
        self.log(
            f"{stage}_contact_accuracy",
            metrics["contact_accuracy"],
            prog_bar=False,
            on_step=False,
            on_epoch=True,
            batch_size=batch_size,
        )
        self.log(
            f"{stage}_contact_precision",
            metrics["contact_precision"],
            prog_bar=False,
            on_step=False,
            on_epoch=True,
            batch_size=batch_size,
        )
        self.log(
            f"{stage}_contact_recall",
            metrics["contact_recall"],
            prog_bar=False,
            on_step=False,
            on_epoch=True,
            batch_size=batch_size,
        )

    def training_step(self, batch: dict[str, torch.Tensor], _batch_idx: int) -> torch.Tensor:
        """Train on one batch."""
        return self.shared_step(batch, "train")

    def validation_step(self, batch: dict[str, torch.Tensor], _batch_idx: int) -> torch.Tensor:
        """Validate on one batch."""
        return self.shared_step(batch, "val")

    def test_step(self, batch: dict[str, torch.Tensor], _batch_idx: int) -> torch.Tensor:
        """Evaluate on one test batch."""
        return self.shared_step(batch, "test")

    def configure_optimizers(self) -> torch.optim.Optimizer:
        """Use Adam as a simple default optimizer."""
        return torch.optim.Adam(self.parameters(), lr=self.learning_rate)


# =============================================================================
# Training orchestration and model selection
# =============================================================================


def build_data_module_and_model(
    args: argparse.Namespace,
    interaction_path: Path,
    loss_mode: str,
) -> tuple[pl.LightningDataModule, pl.LightningModule]:
    """Create the requested data module and Lightning module."""
    # The fully connected baseline and D-SCRIPT need different tensors, so both
    # the DataModule and LightningModule are selected together here.
    if args.model == "fully_connected":
        input_dim = 2 * args.max_length * ONE_HOT_DEPTH
        data_module = FullyConnectedDataModule(
            interaction_path=interaction_path,
            sequence_path=args.sequence_file,
            max_length=args.max_length,
            batch_size=args.batch_size,
            num_workers=args.num_workers,
            loss_mode=loss_mode,
            limit_rows=args.limit_rows,
        )
        model = FullyConnectedLightningModule(
            input_dim=input_dim,
            hidden_dim=args.hidden_dim,
            dropout=args.dropout,
            learning_rate=args.learning_rate,
        )
        return data_module, model

    # D-SCRIPT is the default path. It reads per-residue embeddings and, where
    # available, real contact maps from /nfs/scratch/pdb_dimers/contact_maps.
    data_module = DScriptDataModule(
        interaction_path=interaction_path,
        embedding_manifest=args.embedding_manifest,
        contact_map_dir=args.contact_map_dir,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        loss_mode=loss_mode,
        limit_rows=args.limit_rows,
    )
    model = DScriptLightningModule(
        embedding_dim=args.embedding_dim,
        projection_dim=args.projection_dim,
        contact_hidden_dim=args.contact_hidden_dim,
        contact_width=args.contact_width,
        projection_dropout=args.projection_dropout,
        learning_rate=args.learning_rate,
        interaction_loss_lambda=args.interaction_loss_lambda,
        contact_threshold=args.contact_threshold,
    )
    return data_module, model


def predict_interaction_logits(
    model: pl.LightningModule,
    batch: dict[str, torch.Tensor],
) -> torch.Tensor:
    """Return one raw interaction logit per pair for either supported model."""
    # Prediction export intentionally mirrors the training/evaluation forward
    # paths. The fully connected baseline receives one fixed vector per pair,
    # while D-SCRIPT receives padded residue embeddings plus boolean masks that
    # identify real residues versus padding.
    if "features" in batch:
        return model(batch["features"])

    return model(
        batch["embedding_1"],
        batch["embedding_2"],
        batch["mask_1"],
        batch["mask_2"],
        batch.get("interaction_pair_mask"),
    )


def make_prediction_dataloader(
    data_module: pl.LightningDataModule,
    dataset: Dataset,
) -> DataLoader:
    """Return an ordered DataLoader for prediction export."""
    kwargs = {
        "batch_size": data_module.batch_size,
        "shuffle": False,
        "num_workers": data_module.num_workers,
    }

    # D-SCRIPT examples have variable residue lengths. The same collate function
    # used during training pads embeddings/maps and creates masks, so prediction
    # sees the exact tensor shapes and padding semantics used for evaluation.
    if isinstance(dataset, DScriptDataset):
        kwargs["collate_fn"] = collate_dscript_batch

    return DataLoader(dataset, **kwargs)


def collect_split_predictions(
    split_name: str,
    dataset: Dataset,
    data_module: pl.LightningDataModule,
    model: pl.LightningModule,
    device: torch.device,
) -> pd.DataFrame:
    """Score one split and return metadata with true/predicted labels."""
    probabilities: list[float] = []
    predicted_labels: list[int] = []

    loader = make_prediction_dataloader(data_module, dataset)
    with torch.inference_mode():
        for batch in loader:
            batch = move_tensor_batch_to_device(batch, device)
            logits = predict_interaction_logits(model, batch)

            # Interaction probabilities are sigmoid(logit) in [0, 1]. The
            # predicted label uses the same threshold as the training metrics:
            # probability >= 0.5 means predicted interaction label 1.
            batch_probabilities = torch.sigmoid(logits).detach().cpu()
            batch_predictions = (batch_probabilities >= INTERACTION_THRESHOLD).to(torch.long)
            probabilities.extend(batch_probabilities.tolist())
            predicted_labels.extend(batch_predictions.tolist())

    metadata = dataset.examples.reset_index(drop=True).copy()
    if len(metadata) != len(probabilities):
        raise RuntimeError(
            f"Prediction count mismatch for {split_name}: "
            f"{len(probabilities)} predictions for {len(metadata)} rows"
        )

    if "split" not in metadata.columns:
        metadata.insert(0, "split", split_name)

    metadata["label"] = pd.to_numeric(metadata["label"], errors="raise").astype("int64")
    metadata["predicted_label"] = predicted_labels
    metadata["predicted_probability"] = probabilities

    # Keep identifiers and labels near the front, while excluding heavy/internal
    # training columns such as raw sequences, embedding paths, and loss weights.
    preferred_columns = [
        "split",
        "entity_pair",
        "new_cluster_pair",
        "entity_1",
        "entity_2",
        "cluster_1",
        "cluster_2",
        "pdb_id",
        "assembly_number",
        "dimer_type",
        "label",
        "predicted_label",
        "predicted_probability",
    ]
    excluded_columns = {
        "sequence_1",
        "sequence_2",
        "embedding_path_1",
        "embedding_path_2",
        "contact_map_path",
        "loss_weight",
    }
    output_columns = [column for column in preferred_columns if column in metadata.columns]
    output_columns.extend(
        column for column in metadata.columns if column not in output_columns and column not in excluded_columns
    )
    return metadata[output_columns]


def save_prediction_table(
    data_module: pl.LightningDataModule,
    model: pl.LightningModule,
    output_path: Path,
) -> None:
    """Save labels, predicted labels, and probabilities for test-set rows."""
    model.eval()
    device = next(model.parameters()).device

    # The prediction table is intentionally restricted to held-out test rows.
    # Train/validation labels are still used for fitting and model selection,
    # but they are not exported here to keep the final list focused on the
    # unbiased evaluation split.
    if data_module.test_dataset is None or len(data_module.test_dataset) == 0:
        print("Skipping prediction export because no test rows are available")
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    predictions = collect_split_predictions("test", data_module.test_dataset, data_module, model, device)
    predictions.to_csv(output_path, sep="\t", index=False)
    print(f"Saved test predictions for {len(predictions)} examples to {output_path}")


def save_contact_map_png(
    matrix: np.ndarray,
    output_path: Path,
    title: str,
) -> None:
    """Save one residue-by-residue contact matrix as a categorical PNG."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.patches import Patch

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Contact maps are categorical, not continuous:
    #   -1 = unknown/unresolved residue pair or padding,
    #    0 = known non-contact,
    #    1 = known contact.
    # Predicted maps are thresholded before this function, so their cells are 0/1.
    label_colors = {
        -1: "#bdbdbd",
        0: "#f7f7f7",
        1: "#d73027",
    }
    color_map = ListedColormap(
        [label_colors[-1], label_colors[0], label_colors[1]],
        name="contact_map_labels",
    )
    color_norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], color_map.N)

    figure, axis = plt.subplots(figsize=(6, 5))
    axis.imshow(
        matrix,
        aspect="auto",
        interpolation="nearest",
        cmap=color_map,
        norm=color_norm,
    )
    axis.set_xlabel("Protein 2 residue index")
    axis.set_ylabel("Protein 1 residue index")
    axis.set_title(title)

    # The same legend is used for true and predicted maps so exported PNGs can be
    # compared directly. Predicted maps may not contain -1, but true maps can.
    legend_handles = [
        Patch(
            facecolor=label_colors[1],
            edgecolor="black",
            linewidth=0.25,
            label="1 contact",
        ),
        Patch(
            facecolor=label_colors[0],
            edgecolor="black",
            linewidth=0.25,
            label="0 no contact",
        ),
        Patch(
            facecolor=label_colors[-1],
            edgecolor="black",
            linewidth=0.25,
            label="-1 unknown",
        ),
    ]
    axis.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.14),
        ncol=3,
        frameon=False,
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(figure)


def contact_map_sample_stem(row: pd.Series, sample_number: int, test_index: int) -> str:
    """Build a readable filename stem for one exported test contact-map sample."""
    # PDB/assembly ids are the clearest structural identifiers for D-SCRIPT
    # contact maps. If they are unavailable, fall back to cluster ids and the
    # original test-set row index.
    if "pdb_id" in row and pd.notna(row.pdb_id):
        assembly = f"_assembly{int(row.assembly_number)}" if pd.notna(row.get("assembly_number")) else ""
        identifier = f"{row.pdb_id}{assembly}"
    elif "new_cluster_pair" in row and pd.notna(row.new_cluster_pair):
        identifier = row.new_cluster_pair
    else:
        identifier = f"test_index_{test_index}"

    return f"sample_{sample_number:02d}_test_index_{test_index}_{safe_filename_part(identifier)}"


def save_test_contact_map_samples(
    data_module: pl.LightningDataModule,
    model: pl.LightningModule,
    output_dir: Path,
    sample_count: int = CONTACT_MAP_SAMPLE_COUNT,
) -> None:
    """Save predicted and true contact maps for random test examples with labels."""
    if not isinstance(data_module.test_dataset, DScriptDataset) or not isinstance(model, DScriptLightningModule):
        print("Skipping contact-map PNG export because contact maps are only available for the D-SCRIPT model")
        return

    test_dataset = data_module.test_dataset
    examples = test_dataset.examples
    candidates = examples.index[examples["contact_map_path"].notna()].to_numpy()
    if len(candidates) == 0:
        print("Skipping contact-map PNG export because no test rows have true contact maps")
        return

    try:
        import matplotlib  # noqa: F401
    except ImportError as exc:
        print(f"Skipping contact-map PNG export because matplotlib is not installed: {exc}")
        return

    # Draw without replacement from test rows that have a true contact map. This
    # avoids negative/synthetic rows where the contact target is entirely unknown
    # and therefore cannot be shown as a ground-truth map.
    sample_size = min(sample_count, len(candidates))
    sampled_indices = np.random.default_rng().choice(candidates, size=sample_size, replace=False)

    model.eval()
    device = next(model.parameters()).device
    for sample_number, test_index in enumerate(sampled_indices, start=1):
        item = test_dataset[int(test_index)]
        batch = collate_dscript_batch([item])
        batch = move_tensor_batch_to_device(batch, device)

        with torch.inference_mode():
            _, predicted_contact_map = model.model(
                batch["embedding_1"],
                batch["embedding_2"],
                batch["mask_1"],
                batch["mask_2"],
                interaction_pair_mask=batch["interaction_pair_mask"],
                return_contact_map=True,
            )

        # With a one-example batch, the collated tensors have no padding beyond
        # the current pair's true lengths. The exported predicted map uses the
        # same probability threshold as contact metrics, giving categorical 0/1
        # cells. True maps keep the dataset convention: -1 unknown, 0 non-contact,
        # and 1 contact.
        predicted_map = (
            predicted_contact_map[0].detach().cpu().numpy() >= model.contact_threshold
        ).astype(np.int8)
        true_map = batch["contact_map"][0].detach().cpu().numpy().astype(np.int8, copy=False)
        row = examples.iloc[int(test_index)]
        stem = contact_map_sample_stem(row, sample_number, int(test_index))

        save_contact_map_png(
            predicted_map,
            output_dir / f"{stem}_predicted_contact_map.png",
            f"Predicted Contact Map (threshold >= {model.contact_threshold:g})",
        )
        save_contact_map_png(
            true_map,
            output_dir / f"{stem}_true_contact_map.png",
            "True Contact Map",
        )

    print(f"Saved predicted and true contact-map PNGs for {sample_size} test examples to {output_dir}")


def main() -> None:
    """Wire the data module, model, trainer, and final test run together."""
    args = parse_args()
    validate_args(args)

    # Set pytorch lightning seed
    pl.seed_everything(42)

    # If no interaction file is passed, train on the strict balanced split. The
    # old keep_homomers file name is still detected for the dimer-type weighting
    # mode so existing baseline runs keep their previous behavior.
    interaction_path = args.interactions or DEFAULT_STRICT_FILE
    loss_mode = args.loss_mode
    if loss_mode == "auto":
        loss_mode = "dimer_type" if "keep_homomers" in interaction_path.name else "basic"

    data_module, model = build_data_module_and_model(args, interaction_path, loss_mode)

    # Explicit setup before Trainer construction lets us decide whether a
    # validation set exists and therefore which checkpointing strategy to use.
    data_module.setup("fit")

    has_validation = data_module.val_dataset is not None and len(data_module.val_dataset) > 0
    if has_validation:
        # With validation data, keep the checkpoint with the lowest validation
        # loss and also save "last" for resuming/debugging.
        checkpoint_callback = ModelCheckpoint(
            dirpath=args.output_dir / "checkpoints",
            filename=f"{args.model}" + "-{epoch:02d}-{val_loss:.4f}",
            monitor="val_loss",
            mode="min",
            save_top_k=1,
            save_last=True,
        )
    else:
        # Without validation data, there is no model-selection metric. Save all
        # epoch checkpoints plus "last" so training output is still recoverable.
        checkpoint_callback = ModelCheckpoint(
            dirpath=args.output_dir / "checkpoints",
            filename=f"{args.model}" + "-{epoch:02d}",
            save_top_k=-1,
            save_last=True,
        )

    csv_logger = CSVLogger(save_dir=args.output_dir, name="metrics")
    trainer = pl.Trainer(
        max_epochs=args.max_epochs,
        accelerator=args.accelerator,
        devices=int(args.devices) if args.devices.isdigit() else args.devices,
        default_root_dir=args.output_dir,
        callbacks=[checkpoint_callback],
        logger=csv_logger,
        log_every_n_steps=args.log_every_n_steps,
    )

    trainer.fit(model, datamodule=data_module)

    if data_module.test_dataset and len(data_module.test_dataset) > 0:
        # If validation existed, test the best validation checkpoint. Otherwise
        # test the current final weights.
        trainer.test(model, datamodule=data_module, ckpt_path="best" if has_validation else None)

    prediction_model = model
    if has_validation and checkpoint_callback.best_model_path:
        # The prediction table should describe the same selected model as the
        # validation-based checkpoint/test result, not necessarily the final
        # epoch weights after training.
        prediction_model = type(model).load_from_checkpoint(checkpoint_callback.best_model_path)
        prediction_model.to(trainer.strategy.root_device)

    csv_logger.save()
    save_loss_plot(Path(csv_logger.log_dir) / "metrics.csv", args.output_dir / "loss_curve.png")
    save_prediction_table(data_module, prediction_model, args.output_dir / "test_predictions.tsv")
    save_test_contact_map_samples(data_module, prediction_model, args.output_dir / "contact_maps")


if __name__ == "__main__":
    main()
