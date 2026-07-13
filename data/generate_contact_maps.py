"""Generate full-sequence inter-chain contact maps for selected PDB dimers.

Each input TSV row describes one positive dimer assembly and one entity pair.
The script resolves that pair to concrete coordinate chains in the assembly
mmCIF file, extracts one representative atom per resolved residue, and writes a
rectangular NumPy contact map.

Map axes and labels:
- rows are the full polymer sequence positions of entity 1,
- columns are the full polymer sequence positions of entity 2,
- values are int8 labels: -1 unknown, 0 known non-contact, 1 known contact.

Unknown cells are kept explicitly. The map is initialized to -1, and only
residue pairs where both residues have usable coordinates are overwritten with
0 or 1. This preserves unresolved residues in their original sequence positions
so later training code can mask unknown labels instead of silently shortening
the proteins.
"""

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


# ---------------------------------------------------------------------------
# Configuration and compact data containers
# ---------------------------------------------------------------------------

DEFAULT_SELECTION_TSV = "/nfs/scratch/pdb_dimers/final_filtered_interactions_with_partitions.tsv"
DEFAULT_INPUT_DIR = "/nfs/scratch/pdb_dimers/assemblies"
DEFAULT_OUTPUT_DIR = "/nfs/scratch/pdb_dimers/contact_maps"
DEFAULT_THRESHOLD_ANGSTROM = 8.0
STRUCTURE_SUFFIXES = (".cif.gz", ".cif")
CONTACT_MAP_DATA_SUBDIR = "data"

# Contact maps are int8 arrays with this complete label set.
CONTACT_UNKNOWN = -1
CONTACT_ABSENT = 0
CONTACT_PRESENT = 1


@dataclass(frozen=True)
class InteractionRecord:
    """One selected TSV row after resolving its structure and output paths."""

    row_index: int
    pdb_id: str
    assembly_number: str
    entity_names: tuple[str, str]
    entity_ids: tuple[str, str]
    structure_path: Path
    output_path: Path


@dataclass(frozen=True)
class ResidueRecord:
    """One resolved residue with the coordinate used for distance checks.

    ``seq_index`` is zero-based and indexes the full polymer sequence axis in
    the output map. ``coordinate`` has shape [3], dtype float32, in Angstrom.
    """

    chain_id: str
    seq_id: str
    seq_index: int
    resname: str
    coordinate: np.ndarray


@dataclass(frozen=True)
class ChainSelection:
    """The concrete coordinate chains selected for one requested entity pair."""

    requested_entity_ids: tuple[str, str]
    remapped_entity_ids: tuple[str, str]
    chain_ids: tuple[str, str]
    sequence_lengths: tuple[int, int]
    first_chain_records: list[ResidueRecord]
    second_chain_records: list[ResidueRecord]


@dataclass(frozen=True)
class ContactMapStats:
    """Counts used by the contact-map audit summary."""

    contacts: int
    known_pairs: int
    unknown_pairs: int


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse paths and contact-definition options for the batch run."""

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
    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help="Folder for contact-map outputs. .npy maps are saved in <output-dir>/data/.",
    )
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


# ---------------------------------------------------------------------------
# mmCIF file loading
# ---------------------------------------------------------------------------

def open_text_file(path: Path):
    """Open plain text and gzip-compressed mmCIF files through one interface."""

    if path.name.endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("rt")


def make_parser() -> MMCIFParser:
    """Create a Biopython parser that exposes label asym and label sequence IDs.

    The selection and _struct_asym metadata use label asym IDs. Using label
    residue numbering also lets ``residue.id`` align with _entity_poly_seq
    sequence positions, which are the axes of the output contact map.
    """

    try:
        return MMCIFParser(QUIET=True, auth_chains=False, auth_residues=False)
    except TypeError as exc:
        raise SystemExit(
            "This script needs a Biopython version whose MMCIFParser supports "
            "auth_chains=False and auth_residues=False."
        ) from exc


def parse_structure_and_metadata(path: Path):
    """Return Biopython coordinates plus the raw mmCIF category dictionary.

    Biopython consumes file handles during parsing, so the structure parser and
    MMCIF2Dict each get their own handle. The structure object is used for atom
    coordinates; the dictionary is used for entity remapping, chain-to-entity
    mapping, and full polymer sequence lengths.
    """

    parser = make_parser()
    with open_text_file(path) as handle:
        structure = parser.get_structure(path.stem, handle)
    with open_text_file(path) as handle:
        mmcif_dict = MMCIF2Dict(handle)
    return structure, mmcif_dict


# ---------------------------------------------------------------------------
# Selection TSV parsing
# ---------------------------------------------------------------------------

def strip_structure_suffix(filename: str) -> str:
    """Remove a supported structure suffix while leaving the basename unchanged."""

    for suffix in STRUCTURE_SUFFIXES:
        if filename.endswith(suffix):
            return filename[: -len(suffix)]
    return filename


def output_path_for_structure(structure_path: Path, output_dir: Path) -> Path:
    """Use the assembly basename under the output directory's data subfolder."""

    return output_dir / CONTACT_MAP_DATA_SUBDIR / f"{strip_structure_suffix(structure_path.name)}.npy"


def normalize_assembly_number(value) -> str:
    """Normalize pandas values so assembly 1 is written as ``"1"``, not ``"1.0"``."""

    assembly_number = str(value).strip()
    if assembly_number.endswith(".0"):
        assembly_number = assembly_number[:-2]
    return assembly_number


def filename_from_row(row: pd.Series) -> str:
    """Read the structure filename from the row or derive the default name."""

    if "local_filename" in row and pd.notna(row["local_filename"]):
        return str(row["local_filename"]).strip()

    pdb_id = str(row["pdb_id"]).strip().lower()
    assembly_number = normalize_assembly_number(row["assembly_number"])
    return f"{pdb_id}-assembly{assembly_number}.cif.gz"


def validate_structure_filename(pdb_id: str, assembly_number: str, filename: str) -> None:
    """Fail early when a TSV row points to a different PDB assembly file."""

    expected_stem = f"{pdb_id.lower()}-assembly{assembly_number}"
    actual_stem = strip_structure_suffix(filename)

    if actual_stem.lower() != expected_stem:
        raise ValueError(
            f"Selection row points to {filename!r}, but pdb_id/assembly_number expect {expected_stem!r}."
        )


def entity_id_from_name(entity_name: str) -> str:
    """Extract the numeric entity ID from names such as ``10BL_1``."""

    entity_name = entity_name.strip()
    if "_" not in entity_name:
        raise ValueError(f"Could not parse entity name {entity_name!r}; expected a value like 10BL_1.")
    return entity_name.rsplit("_", 1)[1].strip()


def entity_pair_from_row(row_index: int, row: pd.Series) -> tuple[tuple[str, str], tuple[str, str]]:
    """Parse the TSV entity-pair column into display names and numeric IDs."""

    entity_names = tuple(part.strip() for part in str(row["entity_pair"]).split(","))
    if len(entity_names) != 2:
        raise ValueError(f"Row {row_index} has invalid entity_pair: {row['entity_pair']!r}")

    entity_ids = (entity_id_from_name(entity_names[0]), entity_id_from_name(entity_names[1]))
    return (entity_names[0], entity_names[1]), entity_ids


def load_interaction_records(args: argparse.Namespace) -> list[InteractionRecord]:
    """Load selected interactions and resolve deterministic output paths.

    Each output map is keyed by the structure filename, so duplicate output
    paths in the selection TSV are treated as an error instead of overwriting a
    previous row. The optional ``--limit`` is applied after reading the TSV and
    before duplicate checks, making small dry runs deterministic.
    """

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
        entity_names, entity_ids = entity_pair_from_row(row_index, row)
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
                entity_names=entity_names,
                entity_ids=entity_ids,
                structure_path=structure_path,
                output_path=output_path,
            )
        )

    return records


# ---------------------------------------------------------------------------
# mmCIF metadata helpers
# ---------------------------------------------------------------------------

def as_list(value) -> list[str]:
    """Normalize MMCIF2Dict values, which may be absent, scalar, or list-like."""

    if value is None:
        return []
    if isinstance(value, list):
        return [str(item) for item in value]
    return [str(value)]


def remap_entity_id(entity_id: str, mmcif_dict: dict) -> str:
    """Apply _pdbx_entity_remapping in the same direction as extract_sequences.py.

    Some assemblies remap original entity IDs before coordinates are assigned to
    asym chains. The TSV stores the requested/original IDs, while _struct_asym
    can refer to the remapped IDs. If no remapping row exists, the ID is already
    usable as-is.
    """

    remapped_ids = as_list(mmcif_dict.get("_pdbx_entity_remapping.entity_id"))
    original_ids = as_list(mmcif_dict.get("_pdbx_entity_remapping.orig_entity_id"))

    for remapped_id, original_id in zip(remapped_ids, original_ids):
        if original_id == entity_id and remapped_id != original_id:
            return remapped_id

    return entity_id


def struct_asym_to_entity_id(mmcif_dict: dict) -> OrderedDict[str, str]:
    """Read the label-asym-chain to entity-ID mapping from _struct_asym."""

    asym_ids = as_list(mmcif_dict.get("_struct_asym.id"))
    entity_ids = as_list(mmcif_dict.get("_struct_asym.entity_id"))

    mapping = OrderedDict()
    for asym_id, entity_id in zip(asym_ids, entity_ids):
        mapping[asym_id] = entity_id
    return mapping


def entity_sequence_lengths(mmcif_dict: dict) -> dict[str, int]:
    """Return full polymer sequence lengths by entity ID.

    Lengths come from _entity_poly_seq instead of observed coordinates, because
    unresolved residues must remain present as unknown rows or columns in the
    output map. Sequence positions are 1-based in mmCIF and become zero-based
    only when indexing arrays.
    """

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


# ---------------------------------------------------------------------------
# Coordinate extraction and chain selection
# ---------------------------------------------------------------------------

def residue_seq_id(residue) -> str:
    """Return a compact sequence-position string for audit/debug output."""

    _, residue_number, insertion_code = residue.id
    insertion_code = insertion_code.strip()
    if insertion_code:
        return f"{residue_number}{insertion_code}"
    return str(residue_number)


def residue_seq_index(residue) -> int | None:
    """Convert Biopython's label sequence ID into a zero-based array index."""

    _, residue_number, _ = residue.id
    if residue_number <= 0:
        return None
    return residue_number - 1


def representative_coordinate(residue, representative_atom: str) -> np.ndarray | None:
    """Select the residue coordinate used for inter-chain distance checks.

    ``representative_atom="ca"`` always uses C-alpha. ``"cbeta"`` uses C-beta
    when present and falls back to C-alpha, which keeps glycine and residues
    missing a C-beta from becoming unknown solely because of the atom choice.
    Returned coordinates have shape [3] and dtype float32.
    """

    if representative_atom == "ca":
        atom_name = "CA"
    else:
        atom_name = "CB" if residue.has_id("CB") else "CA"

    if not residue.has_id(atom_name):
        return None

    return residue[atom_name].get_coord().astype(np.float32)


def residue_records_for_chain(chain, representative_atom: str) -> list[ResidueRecord]:
    """Extract one coordinate-bearing record per sequence position in a chain.

    Residues without the selected representative atom remain unresolved and are
    left as unknown in the final map. If multiple residues map to the same
    zero-based sequence index, the first coordinate in chain order is kept so
    each map cell has exactly one row or column coordinate.
    """

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


def residue_records_by_chain(structure, representative_atom: str) -> OrderedDict[str, list[ResidueRecord]]:
    """Extract coordinate-bearing residues from the first model in assembly order."""

    model = next(structure.get_models())
    chains = OrderedDict()

    for chain in model:
        records = residue_records_for_chain(chain, representative_atom)
        if records:
            chains[str(chain.id).strip() or "_"] = records

    return chains


def coordinate_chains_by_entity(
    mmcif_dict: dict,
    chain_records: OrderedDict[str, list[ResidueRecord]],
) -> OrderedDict[str, list[str]]:
    """Group coordinate-bearing label asym IDs by entity, preserving file order."""

    chains_by_entity = OrderedDict()
    for asym_id, entity_id in struct_asym_to_entity_id(mmcif_dict).items():
        if asym_id in chain_records:
            chains_by_entity.setdefault(entity_id, []).append(asym_id)
    return chains_by_entity


def choose_chain_for_entity(
    record: InteractionRecord,
    entity_id: str,
    candidate_chains: list[str],
    allow_ambiguous: bool,
) -> str:
    """Select one coordinate chain for a heteromer entity.

    Multiple coordinate chains for one requested entity are ambiguous because
    the TSV specifies entities, not chain IDs. By default the script fails so
    the user can inspect the assembly. With ``--allow-ambiguous-chain-selection``
    it keeps the first chain in mmCIF order and records a warning.
    """

    if not candidate_chains:
        raise ValueError(f"Entity {entity_id} has no coordinate chain with a representative atom.")

    if len(candidate_chains) > 1:
        if not allow_ambiguous:
            raise ValueError(
                f"Entity {entity_id} maps to {len(candidate_chains)} coordinate chains "
                f"({candidate_chains}); chain selection is ambiguous."
            )
        print(
            f"Warning: {record.structure_path.name} entity {entity_id} has "
            f"{len(candidate_chains)} coordinate chains; using {candidate_chains[0]}",
            file=sys.stderr,
        )

    return candidate_chains[0]


def choose_homomer_chains(
    record: InteractionRecord,
    entity_id: str,
    candidate_chains: list[str],
    allow_ambiguous: bool,
) -> tuple[str, str]:
    """Select two distinct coordinate chains for a homomeric entity pair."""

    if len(candidate_chains) < 2:
        raise ValueError(
            f"Entity {entity_id} maps to {len(candidate_chains)} coordinate chains; "
            "two distinct chains are required for a homomer contact map."
        )

    if len(candidate_chains) > 2:
        if not allow_ambiguous:
            raise ValueError(
                f"Entity {entity_id} maps to {len(candidate_chains)} coordinate chains "
                f"({candidate_chains}); chain-pair selection is ambiguous."
            )
        print(
            f"Warning: {record.structure_path.name} entity {entity_id} has "
            f"{len(candidate_chains)} coordinate chains; using {candidate_chains[0]} and {candidate_chains[1]}",
            file=sys.stderr,
        )

    return candidate_chains[0], candidate_chains[1]


def records_within_sequence(
    record: InteractionRecord,
    chain_id: str,
    records: list[ResidueRecord],
    sequence_length: int,
) -> list[ResidueRecord]:
    """Keep only coordinates whose zero-based indices fit the full sequence axis."""

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


def select_chains_for_record(
    record: InteractionRecord,
    mmcif_dict: dict,
    structure,
    representative_atom: str,
    allow_ambiguous: bool,
) -> ChainSelection:
    """Resolve a TSV entity pair to two coordinate chains and sequence axes.

    The TSV stores requested entity IDs. mmCIF metadata may remap those IDs, so
    chain lookup uses the remapped IDs. Sequence lengths are taken from full
    polymer metadata, while residue records contain only resolved positions with
    usable representative coordinates.
    """

    chain_records = residue_records_by_chain(structure, representative_atom)
    chains_by_entity = coordinate_chains_by_entity(mmcif_dict, chain_records)
    sequence_lengths_by_entity = entity_sequence_lengths(mmcif_dict)
    remapped_entities = tuple(remap_entity_id(entity_id, mmcif_dict) for entity_id in record.entity_ids)

    sequence_lengths = []
    for entity_id in remapped_entities:
        if entity_id not in sequence_lengths_by_entity:
            raise ValueError(f"Entity {entity_id} has no full polymer sequence in _entity_poly_seq.")
        sequence_lengths.append(sequence_lengths_by_entity[entity_id])
    sequence_lengths = (sequence_lengths[0], sequence_lengths[1])

    if remapped_entities[0] == remapped_entities[1]:
        candidate_chains = chains_by_entity.get(remapped_entities[0], [])
        selected_chain_ids = choose_homomer_chains(
            record,
            remapped_entities[0],
            candidate_chains,
            allow_ambiguous,
        )
    else:
        selected_chain_ids = (
            choose_chain_for_entity(
                record,
                remapped_entities[0],
                chains_by_entity.get(remapped_entities[0], []),
                allow_ambiguous,
            ),
            choose_chain_for_entity(
                record,
                remapped_entities[1],
                chains_by_entity.get(remapped_entities[1], []),
                allow_ambiguous,
            ),
        )

    first_records = records_within_sequence(
        record,
        selected_chain_ids[0],
        chain_records[selected_chain_ids[0]],
        sequence_lengths[0],
    )
    second_records = records_within_sequence(
        record,
        selected_chain_ids[1],
        chain_records[selected_chain_ids[1]],
        sequence_lengths[1],
    )

    return ChainSelection(
        requested_entity_ids=record.entity_ids,
        remapped_entity_ids=(remapped_entities[0], remapped_entities[1]),
        chain_ids=selected_chain_ids,
        sequence_lengths=sequence_lengths,
        first_chain_records=first_records,
        second_chain_records=second_records,
    )


# ---------------------------------------------------------------------------
# Contact-map construction
# ---------------------------------------------------------------------------

def coordinates_from_records(records: list[ResidueRecord]) -> np.ndarray:
    """Stack residue coordinates into an array shaped [resolved_residues, 3]."""

    return np.array([record.coordinate for record in records], dtype=np.float32)


def build_contact_map(
    row_records: list[ResidueRecord],
    column_records: list[ResidueRecord],
    row_sequence_length: int,
    column_sequence_length: int,
    threshold: float,
    chunk_size: int,
) -> np.ndarray:
    """Build an int8 contact map in full sequence coordinates.

    Output shape is [row_sequence_length, column_sequence_length]. It starts as
    all -1 so every pair involving an unresolved residue remains unknown. For
    resolved row and column residues, Euclidean distances are thresholded:
    <= ``threshold`` Angstrom becomes 1, and > ``threshold`` becomes 0.

    Distances are chunked along the row axis. Each SciPy ``cdist`` result has
    shape [chunk_resolved_rows, resolved_columns], which is then scattered back
    into the full map using the original zero-based sequence indices.
    """

    contact_map = np.full(
        (row_sequence_length, column_sequence_length),
        CONTACT_UNKNOWN,
        dtype=np.int8,
    )
    row_coordinates = coordinates_from_records(row_records)
    column_coordinates = coordinates_from_records(column_records)
    row_indices = np.array([record.seq_index for record in row_records], dtype=np.intp)
    column_indices = np.array([record.seq_index for record in column_records], dtype=np.intp)

    for start in range(0, len(row_coordinates), chunk_size):
        end = min(start + chunk_size, len(row_coordinates))
        distances = cdist(row_coordinates[start:end], column_coordinates, metric="euclidean")
        contact_labels = np.where(distances <= threshold, CONTACT_PRESENT, CONTACT_ABSENT).astype(np.int8)
        contact_map[row_indices[start:end, None], column_indices[None, :]] = contact_labels

    return contact_map


def summarize_contact_map(contact_map: np.ndarray) -> ContactMapStats:
    """Count contact, known, and unknown labels in a finished contact map."""

    return ContactMapStats(
        contacts=int((contact_map == CONTACT_PRESENT).sum()),
        known_pairs=int((contact_map != CONTACT_UNKNOWN).sum()),
        unknown_pairs=int((contact_map == CONTACT_UNKNOWN).sum()),
    )


# ---------------------------------------------------------------------------
# Per-record processing and output
# ---------------------------------------------------------------------------

def process_record(
    record: InteractionRecord,
    args: argparse.Namespace,
) -> tuple[tuple[int, int], ChainSelection, ContactMapStats]:
    """Process one selected interaction row and save its contact map."""

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

    contact_map = build_contact_map(
        selection.first_chain_records,
        selection.second_chain_records,
        selection.sequence_lengths[0],
        selection.sequence_lengths[1],
        args.threshold,
        args.chunk_size,
    )
    stats = summarize_contact_map(contact_map)

    record.output_path.parent.mkdir(parents=True, exist_ok=True)
    np.save(record.output_path, contact_map)
    return contact_map.shape, selection, stats


def summary_row(
    record: InteractionRecord,
    shape: tuple[int, int],
    selection: ChainSelection,
    stats: ContactMapStats,
) -> dict:
    """Build one TSV summary row for chain selection and map-label counts."""

    return {
        "row_index": record.row_index,
        "pdb_id": record.pdb_id,
        "assembly_number": record.assembly_number,
        "structure_file": record.structure_path.name,
        "output_file": f"{CONTACT_MAP_DATA_SUBDIR}/{record.output_path.name}",
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


def write_summary(summary_path: Path, summary_rows: list[dict]) -> None:
    """Write the audit summary with a stable column order."""

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


# ---------------------------------------------------------------------------
# Batch orchestration
# ---------------------------------------------------------------------------

def main() -> None:
    """Run contact-map generation over all selected TSV rows."""

    args = parse_args()

    if args.chunk_size <= 0:
        raise ValueError("--chunk-size must be greater than zero.")

    output_dir = Path(args.output_dir)
    summary_path = Path(args.summary_tsv) if args.summary_tsv else output_dir / "contact_map_summary.tsv"

    print(f"Starting contact map generation at {datetime.now().isoformat(timespec='seconds')}")
    print(f"Selection TSV: {args.selection_tsv}")
    print(f"Input directory: {args.input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Contact-map directory: {output_dir / CONTACT_MAP_DATA_SUBDIR}")

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
            summary_rows.append(summary_row(record, shape, selection, stats))
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
