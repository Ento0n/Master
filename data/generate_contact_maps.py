import argparse
import csv
import gzip
import sys
from collections import OrderedDict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from Bio.PDB import MMCIFParser
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
except ImportError as exc:
    raise SystemExit("Biopython is required. Install it with: pip install biopython") from exc

try:
    from scipy.spatial.distance import cdist
except ImportError as exc:
    raise SystemExit("SciPy is required. Install it with: pip install scipy") from exc


DEFAULT_SELECTION_TSV = "/nfs/scratch/pdb_dimers/final_filtered_interactions_with_partitions.tsv"
DEFAULT_INPUT_DIR = "/nfs/scratch/pdb_dimers/assemblies"
DEFAULT_OUTPUT_DIR = "/nfs/scratch/pdb_dimers/contact_maps"
DEFAULT_THRESHOLD_ANGSTROM = 8.0


@dataclass(frozen=True)
class InteractionRecord:
    row_index: int
    pdb_id: str
    assembly_number: str
    entity_names: tuple[str, str]
    entity_ids: tuple[str, str]
    structure_path: Path
    output_path: Path


@dataclass(frozen=True)
class ResidueRecord:
    chain_id: str
    seq_id: str
    seq_index: int
    resname: str
    coordinate: np.ndarray


@dataclass(frozen=True)
class ChainSelection:
    requested_entity_ids: tuple[str, str]
    remapped_entity_ids: tuple[str, str]
    chain_ids: tuple[str, str]
    sequence_lengths: tuple[int, int]
    first_chain_records: list[ResidueRecord]
    second_chain_records: list[ResidueRecord]


@dataclass(frozen=True)
class ContactMapStats:
    contacts: int
    known_pairs: int
    unknown_pairs: int


# Parse command-line options for metadata input, structure input, output, and contact definition.
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate controlled full-sequence inter-chain contact maps from selected PDB dimer assembly rows. "
            "Output .npy arrays are int8 with -1 for unknown residue pairs, 0 for known non-contacts, "
            "and 1 for known contacts."
        )
    )
    parser.add_argument(
        "--selection-tsv",
        default=DEFAULT_SELECTION_TSV,
        help="TSV with pdb_id, assembly_number, entity_pair, and optionally local_filename.",
    )
    parser.add_argument("--input-dir", default=DEFAULT_INPUT_DIR, help="Folder containing .cif or .cif.gz files.")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR, help="Folder where .npy contact maps are saved.")
    parser.add_argument(
        "--summary-tsv",
        default=None,
        help="Optional mapping summary path. Default: <output-dir>/contact_map_summary.tsv.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=DEFAULT_THRESHOLD_ANGSTROM,
        help="Residue contact distance threshold in Angstrom. Default: 8.0.",
    )
    parser.add_argument(
        "--map-type",
        choices=["interface"],
        default="interface",
        help="Inter-chain interface map only. Rows are sequence 1 residues and columns are sequence 2 residues.",
    )
    parser.add_argument(
        "--representative-atom",
        choices=["cbeta", "ca"],
        default="cbeta",
        help="Use C-beta atoms with C-alpha fallback, or always use C-alpha. Default: cbeta.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=1024,
        help="Number of row residues processed per SciPy distance calculation. Default: 1024.",
    )
    parser.add_argument(
        "--allow-ambiguous-chain-selection",
        action="store_true",
        help="If an entity maps to too many coordinate chains, use the first matching chains instead of failing.",
    )
    parser.add_argument("--force", action="store_true", help="Regenerate maps that already exist.")
    parser.add_argument("--limit", type=int, default=None, help="Optional number of TSV rows to process.")
    return parser.parse_args()


# Open plain text and gzip-compressed mmCIF files through the same text interface.
def open_text_file(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("rt")


# Create a Biopython mmCIF parser that uses label asym IDs, matching _struct_asym and _atom_site.label_asym_id.
def make_parser() -> MMCIFParser:
    try:
        return MMCIFParser(QUIET=True, auth_chains=False, auth_residues=False)
    except TypeError as exc:
        raise SystemExit(
            "This script needs a Biopython version whose MMCIFParser supports "
            "auth_chains=False and auth_residues=False."
        ) from exc


# Parse one mmCIF file into both Biopython coordinates and a category dictionary.
def parse_structure_and_metadata(path: Path):
    parser = make_parser()
    with open_text_file(path) as handle:
        structure = parser.get_structure(path.stem, handle)
    with open_text_file(path) as handle:
        mmcif_dict = MMCIF2Dict(handle)
    return structure, mmcif_dict


# Keep the assembly basename identical and replace only the structure suffix with .npy.
def output_path_for_structure(structure_path: Path, output_dir: Path) -> Path:
    output_name = structure_path.name
    for suffix in (".cif.gz", ".cif"):
        if output_name.endswith(suffix):
            output_name = output_name[: -len(suffix)]
            break
    return output_dir / f"{output_name}.npy"


# Derive the expected assembly filename when the selection TSV has no local_filename column.
def filename_from_row(row: pd.Series) -> str:
    if "local_filename" in row and pd.notna(row["local_filename"]):
        return str(row["local_filename"]).strip()

    pdb_id = str(row["pdb_id"]).strip().lower()
    assembly_number = normalize_assembly_number(row["assembly_number"])
    return f"{pdb_id}-assembly{assembly_number}.cif.gz"


# Normalize values read by pandas so assembly 1 is always represented as "1", not "1.0".
def normalize_assembly_number(value) -> str:
    assembly_number = str(value).strip()
    if assembly_number.endswith(".0"):
        assembly_number = assembly_number[:-2]
    return assembly_number


# Check that a selected row points to the assembly file named by its pdb_id and assembly_number.
def validate_structure_filename(pdb_id: str, assembly_number: str, filename: str) -> None:
    expected_stem = f"{pdb_id.lower()}-assembly{assembly_number}"
    actual_stem = filename
    for suffix in (".cif.gz", ".cif"):
        if actual_stem.endswith(suffix):
            actual_stem = actual_stem[: -len(suffix)]
            break

    if actual_stem.lower() != expected_stem:
        raise ValueError(
            f"Selection row points to {filename!r}, but pdb_id/assembly_number expect {expected_stem!r}."
        )


# Extract the numeric entity ID from names such as 10BL_1.
def entity_id_from_name(entity_name: str) -> str:
    entity_name = entity_name.strip()
    if "_" not in entity_name:
        raise ValueError(f"Could not parse entity name {entity_name!r}; expected a value like 10BL_1.")
    return entity_name.rsplit("_", 1)[1].strip()


# Load selected interaction rows; each row controls one output contact map.
def load_interaction_records(args: argparse.Namespace) -> list[InteractionRecord]:
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    df = pd.read_csv(args.selection_tsv, sep="\t")

    required_columns = {"pdb_id", "assembly_number", "entity_pair"}
    missing_columns = required_columns.difference(df.columns)
    if missing_columns:
        raise ValueError(f"Selection TSV is missing required columns: {sorted(missing_columns)}")

    if args.limit is not None:
        df = df.head(args.limit).reset_index(drop=True)

    records = []
    seen_outputs = set()
    for row_index, row in df.iterrows():
        entity_names = tuple(part.strip() for part in str(row["entity_pair"]).split(","))
        if len(entity_names) != 2:
            raise ValueError(f"Row {row_index} has invalid entity_pair: {row['entity_pair']!r}")

        pdb_id = str(row["pdb_id"]).strip()
        assembly_number = normalize_assembly_number(row["assembly_number"])
        filename = filename_from_row(row)
        validate_structure_filename(pdb_id, assembly_number, filename)
        structure_path = input_dir / filename
        output_path = output_path_for_structure(structure_path, output_dir)
        if output_path in seen_outputs:
            raise ValueError(f"Duplicate output path from selection TSV: {output_path}")
        seen_outputs.add(output_path)

        records.append(
            InteractionRecord(
                row_index=row_index,
                pdb_id=pdb_id,
                assembly_number=assembly_number,
                entity_names=(entity_names[0], entity_names[1]),
                entity_ids=(entity_id_from_name(entity_names[0]), entity_id_from_name(entity_names[1])),
                structure_path=structure_path,
                output_path=output_path,
            )
        )

    return records


# Normalize MMCIF2Dict values because single values and loop values may be represented differently.
def as_list(value) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item) for item in value]
    return [str(value)]


# Apply _pdbx_entity_remapping in the same direction as extract_sequences.py.
def remap_entity_id(entity_id: str, mmcif_dict: dict) -> str:
    remapped_ids = as_list(mmcif_dict.get("_pdbx_entity_remapping.entity_id"))
    original_ids = as_list(mmcif_dict.get("_pdbx_entity_remapping.orig_entity_id"))

    for remapped_id, original_id in zip(remapped_ids, original_ids):
        if original_id == entity_id and remapped_id != original_id:
            return remapped_id

    return entity_id


# Read the structure-asym mapping that connects entity IDs to concrete coordinate chain IDs.
def struct_asym_to_entity_id(mmcif_dict: dict) -> OrderedDict[str, str]:
    asym_ids = as_list(mmcif_dict.get("_struct_asym.id"))
    entity_ids = as_list(mmcif_dict.get("_struct_asym.entity_id"))

    mapping = OrderedDict()
    for asym_id, entity_id in zip(asym_ids, entity_ids):
        mapping[asym_id] = entity_id
    return mapping


# Read full polymer sequence lengths so unresolved residues can stay in their original sequence positions.
def entity_sequence_lengths(mmcif_dict: dict) -> dict[str, int]:
    entity_ids = as_list(mmcif_dict.get("_entity_poly_seq.entity_id"))
    seq_nums = as_list(mmcif_dict.get("_entity_poly_seq.num"))

    lengths = {}
    for entity_id, seq_num in zip(entity_ids, seq_nums):
        try:
            sequence_position = int(seq_num)
        except ValueError:
            continue
        lengths[entity_id] = max(lengths.get(entity_id, 0), sequence_position)

    return lengths


# Convert a Biopython residue ID into a compact sequence-position string.
def residue_seq_id(residue) -> str:
    _, residue_number, insertion_code = residue.id
    insertion_code = insertion_code.strip()
    if insertion_code:
        return f"{residue_number}{insertion_code}"
    return str(residue_number)


def residue_seq_index(residue) -> int | None:
    _, residue_number, _ = residue.id
    if residue_number <= 0:
        return None
    return residue_number - 1


# Select the residue coordinate used for contact calculation.
def representative_coordinate(residue, representative_atom: str) -> np.ndarray | None:
    if representative_atom == "ca":
        atom_name = "CA"
    else:
        atom_name = "CB" if residue.has_id("CB") else "CA"

    if not residue.has_id(atom_name):
        return None

    return residue[atom_name].get_coord().astype(np.float32)


# Convert one selected Biopython chain into residue records with one coordinate per residue.
def residue_records_for_chain(chain, representative_atom: str) -> list[ResidueRecord]:
    records = []
    chain_id = str(chain.id).strip() or "_"
    seen_seq_indices = set()

    for residue in chain:
        coordinate = representative_coordinate(residue, representative_atom)
        if coordinate is None:
            continue
        seq_index = residue_seq_index(residue)
        if seq_index is None or seq_index in seen_seq_indices:
            continue
        seen_seq_indices.add(seq_index)

        records.append(
            ResidueRecord(
                chain_id=chain_id,
                seq_id=residue_seq_id(residue),
                seq_index=seq_index,
                resname=residue.get_resname(),
                coordinate=coordinate,
            )
        )

    return records


# Extract coordinate-bearing chains from the first model, preserving the assembly file chain order.
def residue_records_by_chain(structure, representative_atom: str) -> OrderedDict[str, list[ResidueRecord]]:
    model = next(structure.get_models())
    chains = OrderedDict()

    for chain in model:
        records = residue_records_for_chain(chain, representative_atom)
        if records:
            chains[str(chain.id).strip() or "_"] = records

    return chains


# Select the concrete coordinate chains for a TSV row's requested entity pair.
def select_chains_for_record(
    record: InteractionRecord,
    mmcif_dict: dict,
    structure,
    representative_atom: str,
    allow_ambiguous: bool,
) -> ChainSelection:
    asym_to_entity = struct_asym_to_entity_id(mmcif_dict)
    sequence_lengths_by_entity = entity_sequence_lengths(mmcif_dict)
    chain_records = residue_records_by_chain(structure, representative_atom)
    remapped_entities = tuple(remap_entity_id(entity_id, mmcif_dict) for entity_id in record.entity_ids)

    def sequence_length_for_entity(entity_id: str) -> int:
        if entity_id not in sequence_lengths_by_entity:
            raise ValueError(f"Entity {entity_id} has no full polymer sequence in _entity_poly_seq.")
        return sequence_lengths_by_entity[entity_id]

    sequence_lengths = tuple(sequence_length_for_entity(entity_id) for entity_id in remapped_entities)

    def coordinate_chains_for_entity(entity_id: str) -> list[str]:
        candidate_chains = [asym_id for asym_id, mapped_entity_id in asym_to_entity.items() if mapped_entity_id == entity_id]
        return [chain_id for chain_id in candidate_chains if chain_id in chain_records]

    def records_within_sequence(chain_id: str, sequence_length: int) -> list[ResidueRecord]:
        records = chain_records[chain_id]
        in_bounds = [residue_record for residue_record in records if residue_record.seq_index < sequence_length]
        skipped = len(records) - len(in_bounds)
        if skipped:
            print(
                f"Warning: {record.structure_path.name} chain {chain_id} has {skipped} coordinate residues "
                f"outside the {sequence_length}-residue entity sequence; skipping them.",
                file=sys.stderr,
            )
        if not in_bounds:
            raise ValueError(f"Selected chain {chain_id} has no coordinate residues within its entity sequence.")
        return in_bounds

    if remapped_entities[0] == remapped_entities[1]:
        candidate_chains = coordinate_chains_for_entity(remapped_entities[0])
        if len(candidate_chains) < 2:
            raise ValueError(
                f"Entity {remapped_entities[0]} maps to {len(candidate_chains)} coordinate chains; "
                "two distinct chains are required for a homomer contact map."
            )
        if len(candidate_chains) > 2:
            if not allow_ambiguous:
                raise ValueError(
                    f"Entity {remapped_entities[0]} maps to {len(candidate_chains)} coordinate chains "
                    f"({candidate_chains}); chain-pair selection is ambiguous."
                )
            print(
                f"Warning: {record.structure_path.name} entity {remapped_entities[0]} has "
                f"{len(candidate_chains)} coordinate chains; using {candidate_chains[0]} and {candidate_chains[1]}",
                file=sys.stderr,
            )
        selected_chain_ids = (candidate_chains[0], candidate_chains[1])
    else:
        first_candidates = coordinate_chains_for_entity(remapped_entities[0])
        second_candidates = coordinate_chains_for_entity(remapped_entities[1])
        if not first_candidates:
            raise ValueError(f"Entity {remapped_entities[0]} has no coordinate chain with a representative atom.")
        if not second_candidates:
            raise ValueError(f"Entity {remapped_entities[1]} has no coordinate chain with a representative atom.")
        if len(first_candidates) > 1:
            if not allow_ambiguous:
                raise ValueError(
                    f"Entity {remapped_entities[0]} maps to {len(first_candidates)} coordinate chains "
                    f"({first_candidates}); chain selection is ambiguous."
                )
            print(
                f"Warning: {record.structure_path.name} entity {remapped_entities[0]} has "
                f"{len(first_candidates)} coordinate chains; using {first_candidates[0]}",
                file=sys.stderr,
            )
        if len(second_candidates) > 1:
            if not allow_ambiguous:
                raise ValueError(
                    f"Entity {remapped_entities[1]} maps to {len(second_candidates)} coordinate chains "
                    f"({second_candidates}); chain selection is ambiguous."
                )
            print(
                f"Warning: {record.structure_path.name} entity {remapped_entities[1]} has "
                f"{len(second_candidates)} coordinate chains; using {second_candidates[0]}",
                file=sys.stderr,
            )
        selected_chain_ids = (first_candidates[0], second_candidates[0])

    first_records = records_within_sequence(selected_chain_ids[0], sequence_lengths[0])
    second_records = records_within_sequence(selected_chain_ids[1], sequence_lengths[1])

    return ChainSelection(
        requested_entity_ids=record.entity_ids,
        remapped_entity_ids=(remapped_entities[0], remapped_entities[1]),
        chain_ids=selected_chain_ids,
        sequence_lengths=(sequence_lengths[0], sequence_lengths[1]),
        first_chain_records=first_records,
        second_chain_records=second_records,
    )


# Stack residue coordinates into the matrix format expected by SciPy cdist.
def coordinates_from_records(records: list[ResidueRecord]) -> np.ndarray:
    return np.array([record.coordinate for record in records], dtype=np.float32)


# Compute a padded int8 interface map. -1 means unknown because at least one residue is unresolved.
def build_contact_map(
    row_records: list[ResidueRecord],
    column_records: list[ResidueRecord],
    row_sequence_length: int,
    column_sequence_length: int,
    threshold: float,
    chunk_size: int,
) -> np.ndarray:
    contact_map = np.full((row_sequence_length, column_sequence_length), -1, dtype=np.int8)
    row_coordinates = coordinates_from_records(row_records)
    column_coordinates = coordinates_from_records(column_records)
    row_indices = np.array([record.seq_index for record in row_records], dtype=np.intp)
    column_indices = np.array([record.seq_index for record in column_records], dtype=np.intp)

    for start in range(0, len(row_coordinates), chunk_size):
        end = min(start + chunk_size, len(row_coordinates))
        distances = cdist(row_coordinates[start:end], column_coordinates, metric="euclidean")
        contact_map[row_indices[start:end, None], column_indices[None, :]] = (distances <= threshold).astype(np.int8)

    return contact_map


# Build a rectangular selected-chain interface map in full sequence coordinates.
def make_interface_map(selection: ChainSelection, threshold: float, chunk_size: int) -> np.ndarray:
    return build_contact_map(
        selection.first_chain_records,
        selection.second_chain_records,
        selection.sequence_lengths[0],
        selection.sequence_lengths[1],
        threshold,
        chunk_size,
    )


# Process one selected interaction row and write its inter-chain contact map as a NumPy int8 array.
def process_record(record: InteractionRecord, args: argparse.Namespace) -> tuple[tuple[int, int], ChainSelection, ContactMapStats]:
    if not record.structure_path.exists():
        raise FileNotFoundError(f"Assembly file does not exist: {record.structure_path}")

    structure, mmcif_dict = parse_structure_and_metadata(record.structure_path)
    selection = select_chains_for_record(
        record,
        mmcif_dict,
        structure,
        args.representative_atom,
        args.allow_ambiguous_chain_selection,
    )

    contact_map = make_interface_map(selection, args.threshold, args.chunk_size)
    stats = ContactMapStats(
        contacts=int((contact_map == 1).sum()),
        known_pairs=int((contact_map != -1).sum()),
        unknown_pairs=int((contact_map == -1).sum()),
    )

    record.output_path.parent.mkdir(parents=True, exist_ok=True)
    np.save(record.output_path, contact_map)
    return contact_map.shape, selection, stats


# Write a compact audit trail of which TSV entity pair became which concrete coordinate-chain pair.
def write_summary(summary_path: Path, summary_rows: list[dict]) -> None:
    if not summary_rows:
        return

    summary_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "row_index",
        "pdb_id",
        "assembly_number",
        "structure_file",
        "output_file",
        "requested_entity_name_1",
        "requested_entity_name_2",
        "requested_entity_1",
        "requested_entity_2",
        "remapped_entity_1",
        "remapped_entity_2",
        "chain_1",
        "chain_2",
        "residues_1",
        "residues_2",
        "resolved_residues_1",
        "resolved_residues_2",
        "missing_residues_1",
        "missing_residues_2",
        "known_pairs",
        "unknown_pairs",
        "contacts",
        "shape",
    ]

    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)


# Run the batch job over selected interaction rows and report failures without stopping the whole run.
def main() -> None:
    args = parse_args()

    if args.chunk_size <= 0:
        raise ValueError("--chunk-size must be greater than zero.")

    output_dir = Path(args.output_dir)
    summary_path = Path(args.summary_tsv) if args.summary_tsv else output_dir / "contact_map_summary.tsv"

    print(f"Starting contact map generation at {datetime.now().isoformat(timespec='seconds')}")
    print(f"Selection TSV: {args.selection_tsv}")
    print(f"Input directory: {args.input_dir}")
    print(f"Output directory: {output_dir}")

    records = load_interaction_records(args)
    if not records:
        raise ValueError("No interaction rows found in the selection TSV.")

    failures = []
    summary_rows = []
    for index, record in enumerate(records, start=1):
        if record.output_path.exists() and not args.force:
            print(f"[{index}/{len(records)}] Skipping existing map: {record.output_path.name}")
            continue

        try:
            shape, selection, stats = process_record(record, args)
            summary_rows.append(
                {
                    "row_index": record.row_index,
                    "pdb_id": record.pdb_id,
                    "assembly_number": record.assembly_number,
                    "structure_file": record.structure_path.name,
                    "output_file": record.output_path.name,
                    "requested_entity_name_1": record.entity_names[0],
                    "requested_entity_name_2": record.entity_names[1],
                    "requested_entity_1": selection.requested_entity_ids[0],
                    "requested_entity_2": selection.requested_entity_ids[1],
                    "remapped_entity_1": selection.remapped_entity_ids[0],
                    "remapped_entity_2": selection.remapped_entity_ids[1],
                    "chain_1": selection.chain_ids[0],
                    "chain_2": selection.chain_ids[1],
                    "residues_1": selection.sequence_lengths[0],
                    "residues_2": selection.sequence_lengths[1],
                    "resolved_residues_1": len(selection.first_chain_records),
                    "resolved_residues_2": len(selection.second_chain_records),
                    "missing_residues_1": selection.sequence_lengths[0] - len(selection.first_chain_records),
                    "missing_residues_2": selection.sequence_lengths[1] - len(selection.second_chain_records),
                    "known_pairs": stats.known_pairs,
                    "unknown_pairs": stats.unknown_pairs,
                    "contacts": stats.contacts,
                    "shape": f"{shape[0]}x{shape[1]}",
                }
            )
            print(
                f"[{index}/{len(records)}] Saved {record.output_path.name} "
                f"from chains {selection.chain_ids[0]} x {selection.chain_ids[1]} "
                f"with shape {shape[0]} x {shape[1]}, {stats.contacts} contacts, "
                f"{stats.unknown_pairs} unknown pairs"
            )
        except Exception as exc:
            failures.append((record.structure_path.name, str(exc)))
            print(f"[{index}/{len(records)}] Failed {record.structure_path.name}: {exc}", file=sys.stderr)

    write_summary(summary_path, summary_rows)
    if summary_rows:
        print(f"Mapping summary saved to: {summary_path}")

    if failures:
        print("\nFailed structures:", file=sys.stderr)
        for name, message in failures:
            print(f"- {name}: {message}", file=sys.stderr)
        sys.exit(1)

    print(f"Finished contact map generation at {datetime.now().isoformat(timespec='seconds')}")


if __name__ == "__main__":
    main()
