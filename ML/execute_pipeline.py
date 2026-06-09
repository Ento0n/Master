#!/usr/bin/env python
"""Train a simple one-hot protein interaction classifier with PyTorch Lightning.

Example:
    python ML/execute_pipeline.py \
        --interactions /nfs/scratch/pdb_dimers/balanced_interactions_keep_homomers.tsv \
        --sequence-file /nfs/scratch/pdb_dimers/entity_sequences.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Optional, Tuple

import pandas as pd

import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset

import lightning.pytorch as pl
from lightning.pytorch.callbacks import ModelCheckpoint

from models import FullyConnectedInteractionModel


AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AMINO_ACID_TO_INDEX = {amino_acid: index for index, amino_acid in enumerate(AMINO_ACIDS)}
UNKNOWN_AMINO_ACID_INDEX = len(AMINO_ACIDS)
ONE_HOT_DEPTH = len(AMINO_ACIDS) + 1
DIMER_TYPE_TO_INDEX = {"homo": 0, "hetero": 1}

DEFAULT_STRICT_FILE = Path("/nfs/scratch/pdb_dimers/balanced_interactions_strict.tsv")
DEFAULT_SEQUENCE_FILE = Path("/nfs/scratch/pdb_dimers/entity_sequences.tsv")


def parse_args() -> argparse.Namespace:
    """Collect all command-line settings for data, model, and training."""
    parser = argparse.ArgumentParser(description="Train a one-hot PDB dimer interaction model.")
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
        help="Entity sequence TSV with columns entity_name and sequence.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("ML/runs/one_hot_interactions"),
        help="Where Lightning logs and checkpoints are written.",
    )
    parser.add_argument("--max-length", type=int, default=1250, help="Pad or truncate sequences to this length.")
    parser.add_argument("--batch-size", type=int, default=16, help="Pairs per training batch.")
    parser.add_argument("--num-workers", type=int, default=4, help="DataLoader worker processes.")
    parser.add_argument("--max-epochs", type=int, default=10, help="Number of training epochs.")
    parser.add_argument("--learning-rate", type=float, default=1e-3, help="Adam learning rate.")
    parser.add_argument("--hidden-dim", type=int, default=512, help="Width of the first hidden layer.")
    parser.add_argument("--dropout", type=float, default=0.2, help="Dropout used between Linear layers.")
    parser.add_argument(
        "--loss-mode",
        choices=("auto", "basic", "dimer_type"),
        default="auto",
        help="'auto' uses dimer-type weighting for keep_homomers and plain BCE otherwise.",
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


def read_interaction_table(path: Path, limit_rows: Optional[int]) -> pd.DataFrame:
    """Read and validate the interaction TSV."""
    needed_columns = {"entity_pair", "split", "label", "dimer_type"}
    interaction_df = pd.read_csv(path, sep="\t", nrows=limit_rows)
    missing_columns = needed_columns.difference(interaction_df.columns)
    if missing_columns:
        raise ValueError(f"{path} is missing columns: {sorted(missing_columns)}")
    return interaction_df


def read_sequence_table(path: Path) -> Dict[str, str]:
    """Read sequences and return a simple entity_name -> sequence dictionary."""
    sequence_df = pd.read_csv(path, sep="\t", usecols=["entity_name", "sequence"])
    sequence_df = sequence_df.dropna(subset=["entity_name", "sequence"])
    sequence_df = sequence_df.drop_duplicates(subset=["entity_name"], keep="first")
    return dict(zip(sequence_df["entity_name"].astype(str), sequence_df["sequence"].astype(str)))


def add_entities_and_sequences(interaction_df: pd.DataFrame, sequence_by_entity: Dict[str, str]) -> pd.DataFrame:
    """Split entity_pair into two entity IDs and attach both protein sequences."""
    examples = interaction_df[["entity_pair", "split", "label", "dimer_type"]].copy()
    entity_columns = examples["entity_pair"].astype(str).str.split(",", n=1, expand=True)

    if entity_columns.shape[1] != 2:
        bad_examples = examples["entity_pair"].head(5)
        raise ValueError(f"Could not split entity_pair values into two entities: {bad_examples.tolist()}")

    missing_second_entity = entity_columns[1].isna()
    if missing_second_entity.any():
        bad_examples = examples.loc[missing_second_entity, "entity_pair"].head(5)
        raise ValueError(f"Could not split some entity_pair values into two entities: {bad_examples.tolist()}")

    examples["entity_1"] = entity_columns[0].str.strip()
    examples["entity_2"] = entity_columns[1].str.strip()
    examples["sequence_1"] = examples["entity_1"].map(sequence_by_entity)
    examples["sequence_2"] = examples["entity_2"].map(sequence_by_entity)

    missing_entity_1 = set(examples.loc[examples["sequence_1"].isna(), "entity_1"])
    missing_entity_2 = set(examples.loc[examples["sequence_2"].isna(), "entity_2"])
    missing_entities = sorted(missing_entity_1.union(missing_entity_2))
    if missing_entities:
        preview = ", ".join(missing_entities[:10])
        raise ValueError(f"Missing sequences for {len(missing_entities)} entities. First examples: {preview}")

    examples["label"] = pd.to_numeric(examples["label"], errors="raise").astype("int64")
    examples["split"] = examples["split"].astype(str).str.lower().replace({"validation": "val", "valid": "val"})
    examples["dimer_type"] = examples["dimer_type"].fillna("unknown").astype(str).str.lower()
    return examples


def choose_loss_mode(requested_mode: str, interaction_path: Path) -> str:
    """Choose plain BCE or dimer-type weighted BCE."""
    if requested_mode != "auto":
        return requested_mode
    if "keep_homomers" in interaction_path.name:
        return "dimer_type"
    return "basic"


def make_weight_map(train_df: pd.DataFrame, loss_mode: str) -> Dict[Tuple[str, int], float]:
    """Build per-sample BCE weights from the training split.

    For keep_homomers, each (dimer_type, label) cell receives the same total
    weight. This directly handles cases like many positive homomers but fewer
    negative homomers.
    """
    if loss_mode == "basic":
        return {}

    counts = train_df.groupby(["dimer_type", "label"]).size()
    total_examples = float(counts.sum())
    number_of_cells = float(len(counts))
    return {
        (dimer_type, int(label)): total_examples / (number_of_cells * float(count))
        for (dimer_type, label), count in counts.items()
    }


def add_loss_weights(examples: pd.DataFrame, weight_map: Dict[Tuple[str, int], float]) -> pd.DataFrame:
    """Attach the sample weight used by binary cross entropy."""
    examples = examples.copy()
    examples["loss_weight"] = [
        weight_map.get((row.dimer_type, int(row.label)), 1.0) for row in examples.itertuples(index=False)
    ]
    examples["loss_weight"] = examples["loss_weight"].astype("float32")
    return examples


def one_hot_encode_sequence(sequence: str, max_length: int) -> torch.Tensor:
    """Convert one protein sequence to a padded one-hot tensor.

    The 20 common amino acids get their own channels, unknown letters like X
    use the final channel, and padding stays all zeros.
    """
    encoded = torch.zeros((max_length, ONE_HOT_DEPTH), dtype=torch.float32)
    trimmed_sequence = sequence.upper()[:max_length]

    for position, amino_acid in enumerate(trimmed_sequence):
        channel = AMINO_ACID_TO_INDEX.get(amino_acid, UNKNOWN_AMINO_ACID_INDEX)
        encoded[position, channel] = 1.0

    return encoded


class InteractionDataset(Dataset):
    """Dataset that returns one flattened one-hot vector per protein pair."""

    def __init__(self, examples: pd.DataFrame, max_length: int) -> None:
        self.examples = examples.reset_index(drop=True)
        self.max_length = max_length

    def __len__(self) -> int:
        """Return the number of protein pairs in this split."""
        return len(self.examples)

    def __getitem__(self, index: int) -> Dict[str, torch.Tensor]:
        """Encode both sequences and return tensors used by the Lightning module."""
        row = self.examples.iloc[index]
        sequence_1 = one_hot_encode_sequence(row.sequence_1, self.max_length).flatten()
        sequence_2 = one_hot_encode_sequence(row.sequence_2, self.max_length).flatten()
        features = torch.cat((sequence_1, sequence_2), dim=0)

        return {
            "features": features,
            "label": torch.tensor(float(row.label), dtype=torch.float32),
            "loss_weight": torch.tensor(float(row.loss_weight), dtype=torch.float32),
            "dimer_type_index": torch.tensor(DIMER_TYPE_TO_INDEX.get(row.dimer_type, -1), dtype=torch.long),
        }


class InteractionDataModule(pl.LightningDataModule):
    """LightningDataModule that reads TSV files and creates train/val/test loaders."""

    def __init__(
        self,
        interaction_path: Path,
        sequence_path: Path,
        max_length: int,
        batch_size: int,
        num_workers: int,
        loss_mode: str,
        limit_rows: Optional[int],
    ) -> None:
        super().__init__()
        self.interaction_path = interaction_path
        self.sequence_path = sequence_path
        self.max_length = max_length
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.loss_mode = loss_mode
        self.limit_rows = limit_rows

        self.train_dataset: Optional[InteractionDataset] = None
        self.val_dataset: Optional[InteractionDataset] = None
        self.test_dataset: Optional[InteractionDataset] = None

    def setup(self, stage: Optional[str] = None) -> None:
        """Load data, connect pairs to sequences, and split into datasets."""
        if self.train_dataset is not None:
            return

        interaction_df = read_interaction_table(self.interaction_path, self.limit_rows)
        sequence_by_entity = read_sequence_table(self.sequence_path)
        examples = add_entities_and_sequences(interaction_df, sequence_by_entity)

        train_df = examples[examples["split"] == "train"].copy()
        if train_df.empty:
            raise ValueError("No training rows found. Expected split column to contain 'train'.")

        weight_map = make_weight_map(train_df, self.loss_mode)
        examples = add_loss_weights(examples, weight_map)

        self.train_dataset = InteractionDataset(examples[examples["split"] == "train"], self.max_length)
        self.val_dataset = InteractionDataset(examples[examples["split"] == "val"], self.max_length)
        self.test_dataset = InteractionDataset(examples[examples["split"] == "test"], self.max_length)

        print(
            "Prepared examples: "
            f"train={len(self.train_dataset)}, val={len(self.val_dataset)}, test={len(self.test_dataset)}, "
            f"loss_mode={self.loss_mode}"
        )

    def train_dataloader(self) -> DataLoader:
        """Return shuffled training batches."""
        return DataLoader(
            self.train_dataset,
            batch_size=self.batch_size,
            shuffle=True,
            num_workers=self.num_workers,
        )

    def val_dataloader(self) -> Optional[DataLoader]:
        """Return validation batches when a validation split exists."""
        if self.val_dataset is None or len(self.val_dataset) == 0:
            return None
        return DataLoader(
            self.val_dataset,
            batch_size=self.batch_size,
            shuffle=False,
            num_workers=self.num_workers,
        )

    def test_dataloader(self) -> Optional[DataLoader]:
        """Return test batches when a test split exists."""
        if self.test_dataset is None or len(self.test_dataset) == 0:
            return None
        return DataLoader(
            self.test_dataset,
            batch_size=self.batch_size,
            shuffle=False,
            num_workers=self.num_workers,
        )


class InteractionLightningModule(pl.LightningModule):
    """Lightning wrapper around the fully connected interaction model."""

    def __init__(self, input_dim: int, hidden_dim: int, dropout: float, learning_rate: float) -> None:
        super().__init__()
        self.save_hyperparameters()
        self.model = FullyConnectedInteractionModel(input_dim=input_dim, hidden_dim=hidden_dim, dropout=dropout)
        self.learning_rate = learning_rate

    def forward(self, features: torch.Tensor) -> torch.Tensor:
        """Run a batch through the model and return logits."""
        return self.model(features)

    def shared_step(self, batch: Dict[str, torch.Tensor], stage: str) -> torch.Tensor:
        """Compute weighted BCE loss and simple unweighted accuracy for one batch."""
        labels = batch["label"]
        weights = batch["loss_weight"]
        logits = self(batch["features"])

        per_example_loss = F.binary_cross_entropy_with_logits(logits, labels, reduction="none")
        loss = (per_example_loss * weights).mean()

        probabilities = torch.sigmoid(logits)
        predictions = (probabilities >= 0.5).float()
        accuracy = (predictions == labels).float().mean()

        self.log(f"{stage}_loss", loss, prog_bar=True, on_step=False, on_epoch=True, batch_size=labels.size(0))
        self.log(
            f"{stage}_accuracy",
            accuracy,
            prog_bar=True,
            on_step=False,
            on_epoch=True,
            batch_size=labels.size(0),
        )

        if stage == "test":
            self.log_test_accuracy_by_dimer_type(batch["dimer_type_index"], predictions, labels)

        return loss

    def log_test_accuracy_by_dimer_type(
        self,
        dimer_type_index: torch.Tensor,
        predictions: torch.Tensor,
        labels: torch.Tensor,
    ) -> None:
        """Log separate test accuracies for homo- and heterodimers.

        Each metric is averaged across only the examples from that dimer class,
        so the final test output can show whether the model behaves differently
        on homomers versus heteromers.
        """
        for dimer_type_name, dimer_type_value in DIMER_TYPE_TO_INDEX.items():
            dimer_type_mask = dimer_type_index == dimer_type_value
            dimer_type_count = int(dimer_type_mask.sum().item())

            if dimer_type_count == 0:
                continue

            dimer_type_accuracy = (predictions[dimer_type_mask] == labels[dimer_type_mask]).float().mean()
            self.log(
                f"test_accuracy_{dimer_type_name}",
                dimer_type_accuracy,
                prog_bar=True,
                on_step=False,
                on_epoch=True,
                batch_size=dimer_type_count,
            )

    def training_step(self, batch: Dict[str, torch.Tensor], batch_idx: int) -> torch.Tensor:
        """Train on one batch."""
        return self.shared_step(batch, "train")

    def validation_step(self, batch: Dict[str, torch.Tensor], batch_idx: int) -> torch.Tensor:
        """Validate on one batch."""
        return self.shared_step(batch, "val")

    def test_step(self, batch: Dict[str, torch.Tensor], batch_idx: int) -> torch.Tensor:
        """Evaluate on one test batch."""
        return self.shared_step(batch, "test")

    def configure_optimizers(self) -> torch.optim.Optimizer:
        """Use Adam as a simple default optimizer."""
        return torch.optim.Adam(self.parameters(), lr=self.learning_rate)


def parse_devices(devices: str):
    """Let Lightning accept either 'auto' or a plain integer device count."""
    if devices.isdigit():
        return int(devices)
    return devices


def main() -> None:
    """Wire the data module, model, trainer, and final test run together."""
    args = parse_args()
    interaction_path = args.interactions or DEFAULT_STRICT_FILE
    loss_mode = choose_loss_mode(args.loss_mode, interaction_path)
    input_dim = 2 * args.max_length * ONE_HOT_DEPTH

    data_module = InteractionDataModule(
        interaction_path=interaction_path,
        sequence_path=args.sequence_file,
        max_length=args.max_length,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        loss_mode=loss_mode,
        limit_rows=args.limit_rows,
    )
    data_module.setup("fit")

    has_validation = data_module.val_dataset is not None and len(data_module.val_dataset) > 0
    if has_validation:
        checkpoint_callback = ModelCheckpoint(
            dirpath=args.output_dir / "checkpoints",
            filename="interaction-{epoch:02d}-{val_loss:.4f}",
            monitor="val_loss",
            mode="min",
            save_top_k=1,
            save_last=True,
        )
    else:
        checkpoint_callback = ModelCheckpoint(
            dirpath=args.output_dir / "checkpoints",
            filename="interaction-{epoch:02d}",
            save_top_k=-1,
            save_last=True,
        )

    model = InteractionLightningModule(
        input_dim=input_dim,
        hidden_dim=args.hidden_dim,
        dropout=args.dropout,
        learning_rate=args.learning_rate,
    )

    trainer = pl.Trainer(
        max_epochs=args.max_epochs,
        accelerator=args.accelerator,
        devices=parse_devices(args.devices),
        default_root_dir=args.output_dir,
        callbacks=[checkpoint_callback],
        log_every_n_steps=args.log_every_n_steps,
    )

    trainer.fit(model, datamodule=data_module)

    if data_module.test_dataset and len(data_module.test_dataset) > 0:
        trainer.test(model, datamodule=data_module, ckpt_path="best" if has_validation else None)


if __name__ == "__main__":
    main()
