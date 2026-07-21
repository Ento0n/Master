import argparse
import json
import os
import re
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

import pandas as pd


DEFAULT_INPUT = "/nfs/scratch/pdb_dimers/sequences/unique_sequences.tsv"
DEFAULT_OUTPUT_DIR = "/nfs/scratch/pdb_dimers/embeddings"

MODEL_SPECS = {
    "esm2": {
        "folder": "esm2_t33_650M",
        "hf_model": "facebook/esm2_t33_650M_UR50D",
        "default_max_residues": 1022,
        "description": "ESM2 t33 650M UR50D",
    },
    "prostt5": {
        "folder": "prostt5",
        "hf_model": "Rostlab/ProstT5",
        "default_max_residues": 1022,
        "description": "ProstT5",
    },
}

CANONICAL_WITH_X = set("ACDEFGHIKLMNPQRSTVWYX")


@dataclass(frozen=True)
class SequenceRecord:
    row_index: int
    cluster_id: str
    sequence: str
    output_name: str


@dataclass(frozen=True)
class ChunkRecord:
    record_index: int
    chunk_index: int
    start: int
    end: int
    sequence: str

    @property
    def length(self) -> int:
        return self.end - self.start


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate per-residue and mean embeddings for sequences in "
            "unique_sequences.tsv. Output tensors are named by new_cluster_id."
        )
    )
    parser.add_argument("--input", default=DEFAULT_INPUT, help="Input TSV with sequence and new_cluster_id columns.")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR, help="Root directory for embedding outputs.")
    parser.add_argument(
        "--models",
        nargs="+",
        choices=["all", "esm2", "prostt5"],
        default=["all"],
        help="Models to run. Default: all.",
    )
    parser.add_argument(
        "--device",
        default="auto",
        help="Torch device. Use 'auto', 'cpu', 'cuda', or e.g. 'cuda:0'. Default: auto.",
    )
    parser.add_argument(
        "--dtype",
        choices=["auto", "float32", "float16", "bfloat16"],
        default="auto",
        help="Model compute dtype. Default: float16 on CUDA, float32 on CPU.",
    )
    parser.add_argument(
        "--storage-dtype",
        choices=["float32", "float16", "bfloat16"],
        default="float16",
        help="Dtype used for saved embedding tensors. Default: float16.",
    )
    parser.add_argument("--batch-size", type=int, default=2, help="Maximum number of chunks per forward pass.")
    parser.add_argument(
        "--max-batch-residues",
        type=int,
        default=2048,
        help="Approximate maximum residues per forward pass across all chunks.",
    )
    parser.add_argument(
        "--esm2-max-residues",
        type=int,
        default=MODEL_SPECS["esm2"]["default_max_residues"],
        help="Maximum residues per ESM2 chunk. Use 0 to disable chunking.",
    )
    parser.add_argument(
        "--prostt5-max-residues",
        type=int,
        default=MODEL_SPECS["prostt5"]["default_max_residues"],
        help="Maximum residues per ProstT5 chunk. Use 0 to disable chunking.",
    )
    parser.add_argument(
        "--chunk-overlap",
        type=int,
        default=0,
        help="Residue overlap between chunks. Overlapping residues are averaged. Default: 0.",
    )
    parser.add_argument(
        "--prostt5-prefix",
        default="<AA2fold>",
        help="Task prefix prepended to ProstT5 amino-acid inputs. Use '' to disable.",
    )
    parser.add_argument("--cache-dir", default=None, help="Optional Hugging Face cache directory.")
    parser.add_argument(
        "--local-files-only",
        action="store_true",
        help="Only load models already present in the local Hugging Face cache.",
    )
    parser.add_argument("--force", action="store_true", help="Regenerate outputs even if both tensor files exist.")
    parser.add_argument("--limit", type=int, default=None, help="Optional limit for testing.")
    parser.add_argument(
        "--start-index",
        type=int,
        default=None,
        help="Optional 0-based inclusive row start for array jobs after loading the TSV.",
    )
    parser.add_argument(
        "--end-index",
        type=int,
        default=None,
        help="Optional 0-based exclusive row end for array jobs after loading the TSV.",
    )
    return parser.parse_args()


def timestamp() -> str:
    return datetime.now().isoformat(timespec="seconds")


def resolve_models(model_args: list[str]) -> list[str]:
    if "all" in model_args:
        return ["esm2", "prostt5"]
    return model_args


def safe_output_name(cluster_id: str) -> str:
    name = str(cluster_id).strip()
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
    safe = safe.strip("._")
    if not safe:
        raise ValueError(f"Could not create a safe output filename for cluster id {cluster_id!r}")
    return safe


def load_records(args: argparse.Namespace) -> list[SequenceRecord]:
    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input TSV does not exist: {input_path}")

    df = pd.read_csv(input_path, sep="\t", usecols=["sequence", "new_cluster_id"])
    df = df.dropna(subset=["sequence", "new_cluster_id"]).reset_index(drop=True)

    if args.start_index is not None or args.end_index is not None:
        start = args.start_index if args.start_index is not None else 0
        end = args.end_index if args.end_index is not None else len(df)
        df = df.iloc[start:end].reset_index(drop=True)

    if args.limit is not None:
        df = df.head(args.limit).reset_index(drop=True)

    if df["new_cluster_id"].duplicated().any():
        duplicates = df.loc[df["new_cluster_id"].duplicated(), "new_cluster_id"].head(10).tolist()
        raise ValueError(f"new_cluster_id must be unique. Example duplicates: {duplicates}")

    records = []
    for row_index, row in df.iterrows():
        sequence = str(row["sequence"]).strip().upper()
        cluster_id = str(row["new_cluster_id"]).strip()
        if not sequence:
            continue
        records.append(
            SequenceRecord(
                row_index=row_index,
                cluster_id=cluster_id,
                sequence=sequence,
                output_name=safe_output_name(cluster_id),
            )
        )
    return records


def model_dirs(output_dir: Path, model_key: str) -> tuple[Path, Path, Path]:
    model_dir = output_dir / MODEL_SPECS[model_key]["folder"]
    per_residue_dir = model_dir / "per_residue"
    mean_dir = model_dir / "mean"
    return model_dir, per_residue_dir, mean_dir


def output_paths(output_dir: Path, model_key: str, record: SequenceRecord) -> tuple[Path, Path]:
    _, per_residue_dir, mean_dir = model_dirs(output_dir, model_key)
    return per_residue_dir / f"{record.output_name}.pt", mean_dir / f"{record.output_name}.pt"


def ensure_dirs(output_dir: Path, model_key: str) -> None:
    model_dir, per_residue_dir, mean_dir = model_dirs(output_dir, model_key)
    model_dir.mkdir(parents=True, exist_ok=True)
    per_residue_dir.mkdir(parents=True, exist_ok=True)
    mean_dir.mkdir(parents=True, exist_ok=True)


def sanitize_for_model(sequence: str, model_key: str, tokenizer=None) -> tuple[str, int]:
    if model_key == "prostt5":
        sanitized = "".join(residue if residue in CANONICAL_WITH_X else "X" for residue in sequence)
        return sanitized, sum(a != b for a, b in zip(sequence, sanitized))

    if model_key == "esm2" and tokenizer is not None:
        vocab = tokenizer.get_vocab()
        allowed = {residue for residue in set(sequence) if residue in vocab}
        allowed.add("X")
        sanitized = "".join(residue if residue in allowed else "X" for residue in sequence)
        return sanitized, sum(a != b for a, b in zip(sequence, sanitized))

    sanitized = "".join(residue if residue in CANONICAL_WITH_X else "X" for residue in sequence)
    return sanitized, sum(a != b for a, b in zip(sequence, sanitized))


def max_residues_for_model(args: argparse.Namespace, model_key: str, longest_sequence: int) -> int:
    value = args.esm2_max_residues if model_key == "esm2" else args.prostt5_max_residues
    if value <= 0:
        return longest_sequence
    return value


def make_chunks(record_index: int, sequence: str, max_residues: int, overlap: int) -> list[ChunkRecord]:
    if max_residues <= 0:
        raise ValueError("max_residues must be positive after argument resolution")
    if overlap < 0:
        raise ValueError("--chunk-overlap must be >= 0")
    if overlap >= max_residues:
        raise ValueError("--chunk-overlap must be smaller than the maximum chunk size")

    chunks = []
    start = 0
    chunk_index = 0
    step = max_residues - overlap
    while start < len(sequence):
        end = min(start + max_residues, len(sequence))
        chunks.append(
            ChunkRecord(
                record_index=record_index,
                chunk_index=chunk_index,
                start=start,
                end=end,
                sequence=sequence[start:end],
            )
        )
        if end == len(sequence):
            break
        start += step
        chunk_index += 1
    return chunks


def iter_batches(
    chunks: Iterable[ChunkRecord],
    batch_size: int,
    max_batch_residues: int,
) -> Iterable[list[ChunkRecord]]:
    batch = []
    residues = 0
    for chunk in chunks:
        if batch and (len(batch) >= batch_size or residues + chunk.length > max_batch_residues):
            yield batch
            batch = []
            residues = 0
        batch.append(chunk)
        residues += chunk.length
    if batch:
        yield batch


def resolve_device(torch, device_arg: str):
    if device_arg == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    device = torch.device(device_arg)
    if device.type == "cuda" and not torch.cuda.is_available():
        cuda_visible_devices = os.environ.get("CUDA_VISIBLE_DEVICES", "<unset>")
        raise RuntimeError(
            "CUDA was requested, but this PyTorch installation cannot use CUDA. "
            f"torch={torch.__version__}, torch.version.cuda={torch.version.cuda}, "
            f"CUDA_VISIBLE_DEVICES={cuda_visible_devices}. Install a CUDA-enabled "
            "PyTorch build in the environment used by the Slurm job, or run with --device cpu."
        )
    return device


def resolve_dtype(torch, dtype_arg: str, device) -> object:
    if dtype_arg == "auto":
        return torch.float16 if device.type == "cuda" else torch.float32
    return {
        "float32": torch.float32,
        "float16": torch.float16,
        "bfloat16": torch.bfloat16,
    }[dtype_arg]


def load_model_and_tokenizer(model_key: str, args: argparse.Namespace):
    try:
        import torch
        from transformers import AutoModel, AutoTokenizer, T5EncoderModel
    except ImportError as exc:
        raise ImportError(
            "Embedding generation requires torch and transformers. "
            "Install them in the environment used for the Slurm job."
        ) from exc

    spec = MODEL_SPECS[model_key]
    device = resolve_device(torch, args.device)
    dtype = resolve_dtype(torch, args.dtype, device)
    tokenizer_kwargs = {
        "cache_dir": args.cache_dir,
        "local_files_only": args.local_files_only,
    }
    if model_key == "prostt5":
        tokenizer_kwargs["do_lower_case"] = False
    model_kwargs = {
        "cache_dir": args.cache_dir,
        "local_files_only": args.local_files_only,
    }
    if dtype is not torch.float32:
        model_kwargs["dtype"] = dtype

    print(f"[{timestamp()}] Loading {spec['description']} from {spec['hf_model']} on {device}", flush=True)
    tokenizer = AutoTokenizer.from_pretrained(spec["hf_model"], **tokenizer_kwargs)
    tokenizer.padding_side = "right"

    model_class = T5EncoderModel if model_key == "prostt5" else AutoModel
    model = model_class.from_pretrained(spec["hf_model"], **model_kwargs)
    model.to(device)
    model.eval()

    if device.type == "cuda":
        torch.backends.cuda.matmul.allow_tf32 = True

    return torch, tokenizer, model, device


def format_model_inputs(model_key: str, chunks: list[ChunkRecord], args: argparse.Namespace) -> tuple[list[str], int]:
    if model_key == "esm2":
        return [chunk.sequence for chunk in chunks], 1

    prefix = args.prostt5_prefix.strip()
    spaced_sequences = [" ".join(chunk.sequence) for chunk in chunks]
    if not prefix:
        return spaced_sequences, 0
    return [f"{prefix} {sequence}" for sequence in spaced_sequences], None


def prostt5_prefix_token_count(tokenizer, prefix: str) -> int:
    prefix = prefix.strip()
    if not prefix:
        return 0
    encoded = tokenizer(prefix, add_special_tokens=False)
    return len(encoded["input_ids"])


def embed_batch(model_key: str, chunks: list[ChunkRecord], tokenizer, model, device, args: argparse.Namespace, torch):
    texts, residue_start = format_model_inputs(model_key, chunks, args)
    if model_key == "prostt5" and residue_start is None:
        residue_start = prostt5_prefix_token_count(tokenizer, args.prostt5_prefix)

    inputs = tokenizer(
        texts,
        add_special_tokens=True,
        padding=True,
        truncation=False,
        return_tensors="pt",
    )
    inputs = {key: value.to(device) for key, value in inputs.items()}

    with torch.no_grad():
        outputs = model(**inputs)

    hidden = outputs.last_hidden_state
    embeddings = []
    for batch_index, chunk in enumerate(chunks):
        embedding = hidden[batch_index, residue_start : residue_start + chunk.length].detach().cpu()
        if embedding.shape[0] != chunk.length:
            token_count = int(inputs["attention_mask"][batch_index].sum().item())
            raise RuntimeError(
                f"Unexpected tokenization for {model_key}: chunk length {chunk.length}, "
                f"residue_start {residue_start}, non-padding tokens {token_count}, "
                f"extracted residues {embedding.shape[0]}"
            )
        embeddings.append(embedding)
    return embeddings


def cast_tensor(tensor, dtype_name: str, torch):
    dtype = {
        "float32": torch.float32,
        "float16": torch.float16,
        "bfloat16": torch.bfloat16,
    }[dtype_name]
    return tensor.to(dtype=dtype)


def atomic_torch_save(tensor, path: Path, torch) -> None:
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    torch.save(tensor, tmp_path)
    os.replace(tmp_path, path)


class SequenceAccumulator:
    def __init__(self, length: int, hidden_size: int, torch):
        self.torch = torch
        self.sum = torch.zeros((length, hidden_size), dtype=torch.float32)
        self.count = torch.zeros((length, 1), dtype=torch.float32)

    def add(self, start: int, end: int, embedding) -> None:
        self.sum[start:end] += embedding.to(dtype=self.torch.float32)
        self.count[start:end] += 1.0

    def finalize(self):
        if bool((self.count == 0).any()):
            missing = self.torch.where(self.count.squeeze(1) == 0)[0][:10].tolist()
            raise RuntimeError(f"Missing residue embeddings at positions: {missing}")
        return self.sum / self.count


def tqdm_or_plain(iterable, **kwargs):
    try:
        from tqdm import tqdm

        return tqdm(iterable, **kwargs)
    except ImportError:
        return iterable


def write_model_metadata(model_dir: Path, model_key: str, args: argparse.Namespace, records: list[SequenceRecord]) -> None:
    metadata = {
        "created_at": timestamp(),
        "model_key": model_key,
        "model_description": MODEL_SPECS[model_key]["description"],
        "hf_model": MODEL_SPECS[model_key]["hf_model"],
        "input_tsv": str(Path(args.input)),
        "output_layout": {
            "per_residue": "per_residue/<new_cluster_id>.pt",
            "mean": "mean/<new_cluster_id>.pt",
        },
        "tensor_contents": {
            "per_residue": "[sequence_length, hidden_size] torch tensor, special tokens removed",
            "mean": "[hidden_size] torch tensor, mean over per-residue tensor",
        },
        "storage_dtype": args.storage_dtype,
        "chunk_overlap": args.chunk_overlap,
        "num_records_in_run": len(records),
    }
    if model_key == "prostt5":
        metadata["prostt5_prefix"] = args.prostt5_prefix
    with open(model_dir / "metadata.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")


def write_manifest(output_dir: Path, model_key: str, records: list[SequenceRecord], args: argparse.Namespace) -> None:
    rows = []
    for record in records:
        per_residue_path, mean_path = output_paths(output_dir, model_key, record)
        rows.append(
            {
                "new_cluster_id": record.cluster_id,
                "sequence_length": len(record.sequence),
                "per_residue_path": str(per_residue_path),
                "mean_path": str(mean_path),
                "per_residue_exists": per_residue_path.exists(),
                "mean_exists": mean_path.exists(),
            }
        )
    model_dir, _, _ = model_dirs(output_dir, model_key)
    pd.DataFrame(rows).to_csv(model_dir / "manifest.tsv", sep="\t", index=False)


def run_model(model_key: str, records: list[SequenceRecord], args: argparse.Namespace) -> None:
    output_dir = Path(args.output_dir)
    ensure_dirs(output_dir, model_key)
    model_dir, _, _ = model_dirs(output_dir, model_key)

    pending = []
    for record in records:
        per_residue_path, mean_path = output_paths(output_dir, model_key, record)
        if not args.force and per_residue_path.exists() and mean_path.exists():
            continue
        pending.append(record)

    print(
        f"[{timestamp()}] {MODEL_SPECS[model_key]['description']}: "
        f"{len(pending)} records pending, {len(records) - len(pending)} already complete",
        flush=True,
    )
    write_model_metadata(model_dir, model_key, args, records)
    if not pending:
        write_manifest(output_dir, model_key, records, args)
        return

    torch, tokenizer, model, device = load_model_and_tokenizer(model_key, args)

    sanitized_records = []
    replacements = 0
    for record in pending:
        sequence, changed = sanitize_for_model(record.sequence, model_key, tokenizer=tokenizer)
        replacements += changed
        sanitized_records.append(
            SequenceRecord(
                row_index=record.row_index,
                cluster_id=record.cluster_id,
                sequence=sequence,
                output_name=record.output_name,
            )
        )
    if replacements:
        print(
            f"[{timestamp()}] {MODEL_SPECS[model_key]['description']}: "
            f"replaced {replacements} unsupported residue symbols with X",
            flush=True,
        )

    longest_sequence = max(len(record.sequence) for record in sanitized_records)
    max_residues = max_residues_for_model(args, model_key, longest_sequence)

    all_chunks = []
    remaining_chunks = {}
    for record_index, record in enumerate(sanitized_records):
        chunks = make_chunks(record_index, record.sequence, max_residues, args.chunk_overlap)
        remaining_chunks[record_index] = len(chunks)
        all_chunks.extend(chunks)

    print(
        f"[{timestamp()}] {MODEL_SPECS[model_key]['description']}: "
        f"processing {len(all_chunks)} chunks with max {max_residues} residues/chunk",
        flush=True,
    )

    accumulators = {}
    processed_records = 0
    batches = iter_batches(all_chunks, args.batch_size, args.max_batch_residues)
    progress = tqdm_or_plain(batches, total=None, desc=MODEL_SPECS[model_key]["folder"])
    for batch in progress:
        batch_embeddings = embed_batch(model_key, batch, tokenizer, model, device, args, torch)
        for chunk, embedding in zip(batch, batch_embeddings):
            record = sanitized_records[chunk.record_index]
            if chunk.record_index not in accumulators:
                hidden_size = embedding.shape[1]
                accumulators[chunk.record_index] = SequenceAccumulator(len(record.sequence), hidden_size, torch)

            accumulator = accumulators[chunk.record_index]
            accumulator.add(chunk.start, chunk.end, embedding)
            remaining_chunks[chunk.record_index] -= 1

            if remaining_chunks[chunk.record_index] == 0:
                per_residue = accumulator.finalize()
                mean = per_residue.mean(dim=0)

                original_record = pending[chunk.record_index]
                per_residue_path, mean_path = output_paths(output_dir, model_key, original_record)
                atomic_torch_save(cast_tensor(per_residue, args.storage_dtype, torch), per_residue_path, torch)
                atomic_torch_save(cast_tensor(mean, args.storage_dtype, torch), mean_path, torch)

                del accumulators[chunk.record_index]
                processed_records += 1
                if processed_records % 100 == 0:
                    print(
                        f"[{timestamp()}] {MODEL_SPECS[model_key]['description']}: "
                        f"saved {processed_records}/{len(pending)} records",
                        flush=True,
                    )

    write_manifest(output_dir, model_key, records, args)
    print(
        f"[{timestamp()}] {MODEL_SPECS[model_key]['description']}: "
        f"finished {processed_records} generated records",
        flush=True,
    )


def validate_args(args: argparse.Namespace) -> None:
    if args.batch_size <= 0:
        raise ValueError("--batch-size must be > 0")
    if args.max_batch_residues <= 0:
        raise ValueError("--max-batch-residues must be > 0")
    if args.limit is not None and args.limit <= 0:
        raise ValueError("--limit must be > 0 when provided")
    if args.start_index is not None and args.start_index < 0:
        raise ValueError("--start-index must be >= 0")
    if args.end_index is not None and args.end_index < 0:
        raise ValueError("--end-index must be >= 0")
    if (
        args.start_index is not None
        and args.end_index is not None
        and args.start_index >= args.end_index
    ):
        raise ValueError("--start-index must be smaller than --end-index")


def main() -> None:
    args = parse_args()
    validate_args(args)
    records = load_records(args)
    models = resolve_models(args.models)

    print(f"[{timestamp()}] Loaded {len(records)} sequences from {args.input}", flush=True)
    print(f"[{timestamp()}] Output root: {args.output_dir}", flush=True)
    for model_key in models:
        run_model(model_key, records, args)
    print(f"[{timestamp()}] All requested embedding generation finished", flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[{timestamp()}] ERROR: {exc}", file=sys.stderr, flush=True)
        raise
