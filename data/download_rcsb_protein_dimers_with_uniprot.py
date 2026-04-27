#!/usr/bin/env python3
"""
Download non-redundant RCSB PDB biological assemblies that are protein-only dimers.

Criteria implemented:
  - RCSB Search API return_type="assembly"
  - Entry Polymer Types == "Protein (only)"
  - Biological assembly has exactly two generated polymer instances
  - Those two generated polymer instances are proteins
  - Optional: exclude non-polymer ligands/ions/cofactors from the assembly
  - Redundancy reduction: group by the unordered pair of 100% sequence-identity
    cluster IDs of the two generated protein entities; keep one best-ranked
    assembly per pair.

Default output:
  outdir/assemblies/*.cif.gz
  outdir/selected_assemblies.tsv
  outdir/all_candidate_assemblies.tsv

The TSV files include UniProt accessions for each side of the dimer when RCSB
provides a mapping. Missing mappings are written as empty strings. If one PDB
polymer entity maps to multiple UniProt accessions, they are joined by semicolon.

Requires only the Python standard library.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"
FILES_BASE = "https://files.rcsb.org/download"


@dataclass(frozen=True)
class Candidate:
    assembly_id: str
    pdb_id: str
    assembly_number: str
    cluster_pair: Tuple[str, str]
    entity_pair: Tuple[str, str]
    uniprot_pair: Tuple[str, str]
    dimer_type: str
    resolution: float
    modeled_residue_count: int
    method: str
    oligomeric_count: str
    local_filename: str
    download_url: str

    def rank_key(self) -> Tuple[float, int, str, str]:
        # Lower resolution is best. If no resolution exists (for example NMR),
        # use largest modeled-residue count and then stable IDs as tie-breakers.
        return (self.resolution, -self.modeled_residue_count, self.pdb_id, self.assembly_number)


def terminal(attribute: str, operator: str, value: Any) -> Dict[str, Any]:
    return {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": attribute,
            "operator": operator,
            "value": value,
        },
    }


# Requests a json body via POST from a URL, in this case PDB API, and converts to DICT
def http_json_post(url: str, payload: Dict[str, Any], *, timeout: int = 120, retries: int = 4) -> Dict[str, Any]:
    body = json.dumps(payload).encode("utf-8")
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    last_error: Optional[BaseException] = None
    for attempt in range(retries):
        req = urllib.request.Request(url, data=body, headers=headers, method="POST")
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
            if not data:
                return {}
            return json.loads(data.decode("utf-8"))
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last_error = exc
            if isinstance(exc, urllib.error.HTTPError) and 400 <= exc.code < 500 and exc.code not in (408, 429):
                detail = exc.read().decode("utf-8", errors="replace")
                raise RuntimeError(f"HTTP {exc.code} from {url}: {detail}") from exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"Failed POST {url}: {last_error!r}") from last_error


def download_binary(url: str, path: Path, *, timeout: int = 120, retries: int = 4) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    last_error: Optional[BaseException] = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp, tmp.open("wb") as out:
                while True:
                    chunk = resp.read(1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
            tmp.replace(path)
            return
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
            last_error = exc
            if tmp.exists():
                tmp.unlink()
            if isinstance(exc, urllib.error.HTTPError) and 400 <= exc.code < 500 and exc.code not in (408, 429):
                raise RuntimeError(f"HTTP {exc.code} while downloading {url}") from exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"Failed download {url}: {last_error!r}") from last_error


# Get all assemblies with exactly 2 entities AND Protein (only)
def search_candidate_assemblies(no_ligands: bool) -> List[str]:
    nodes: List[Dict[str, Any]] = [
        terminal("rcsb_entry_info.selected_polymer_entity_types", "exact_match", "Protein (only)"),
        terminal("rcsb_assembly_info.polymer_entity_instance_count", "equals", 2),
        # terminal("rcsb_assembly_info.polymer_entity_instance_count_protein", "equals", 2),
        # terminal("rcsb_assembly_info.selected_polymer_entity_types", "exact_match", "Protein (only)"),
        terminal("rcsb_struct_symmetry.kind", "exact_match", "Global Symmetry"),
    ]
    if no_ligands:
        nodes.append(terminal("rcsb_assembly_info.nonpolymer_entity_instance_count", "equals", 0))

    payload = {
        "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
        "return_type": "assembly",
        "request_options": {
            "return_all_hits": True,
            "results_content_type": ["experimental"],
            "results_verbosity": "compact",
        },
    }
    data = http_json_post(SEARCH_URL, payload)
    result_set = data.get("result_set") or []
    ids: List[str] = []
    for hit in result_set:
        if isinstance(hit, str):
            ids.append(hit.upper())
        elif isinstance(hit, dict) and hit.get("identifier"):
            ids.append(str(hit["identifier"]).upper())
    return sorted(set(ids))


ASSEMBLY_QUERY = """
query AssemblyBatch($assembly_ids: [String!]!) {
  assemblies(assembly_ids: $assembly_ids) {
    rcsb_id
    rcsb_assembly_container_identifiers {
      entry_id
      assembly_id
    }
    rcsb_assembly_info {
      modeled_polymer_monomer_count
      nonpolymer_entity_instance_count
      polymer_entity_instance_count
      polymer_entity_instance_count_protein
    }
    pdbx_struct_assembly {
      oligomeric_count
    }
    pdbx_struct_assembly_gen {
      asym_id_list
      oper_expression
    }
    entry {
      rcsb_id
      rcsb_entry_info {
        experimental_method
        resolution_combined
        selected_polymer_entity_types
      }
      polymer_entities {
        rcsb_id
        entity_poly {
          rcsb_entity_polymer_type
        }
        rcsb_polymer_entity_container_identifiers {
          entity_id
          asym_ids
          auth_asym_ids
          reference_sequence_identifiers {
            database_name
            database_accession
          }
        }
        uniprots {
          rcsb_id
          rcsb_uniprot_container_identifiers {
            uniprot_id
          }
        }
        rcsb_polymer_entity_group_membership {
          aggregation_method
          similarity_cutoff
          group_id
        }
      }
    }
  }
}
"""


def batched(seq: Sequence[str], size: int) -> Iterable[Sequence[str]]:
    for i in range(0, len(seq), size):
        yield seq[i:i + size]


# Retrieve metadata for extracted assembly IDS via GraphQL API from PDB
def fetch_assembly_metadata(assembly_ids: Sequence[str], batch_size: int) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for batch_no, ids in enumerate(batched(assembly_ids, batch_size), start=1):
        payload = {"query": ASSEMBLY_QUERY, "variables": {"assembly_ids": list(ids)}}
        data = http_json_post(GRAPHQL_URL, payload)
        if data.get("errors"):
            raise RuntimeError("GraphQL errors:\n" + json.dumps(data["errors"], indent=2))
        out.extend(((data.get("data") or {}).get("assemblies") or []))
        print(f"Fetched metadata batch {batch_no}: {len(ids)} assemblies", file=sys.stderr)
    return out


def parse_resolution(assembly: Dict[str, Any]) -> float:
    vals = (((assembly.get("entry") or {}).get("rcsb_entry_info") or {}).get("resolution_combined") or [])
    if isinstance(vals, (int, float)):
        return float(vals)
    numeric = []
    for x in vals:
        try:
            if x is not None:
                numeric.append(float(x))
        except (TypeError, ValueError):
            pass
    return min(numeric) if numeric else math.inf


def parse_method(assembly: Dict[str, Any]) -> str:
    vals = (((assembly.get("entry") or {}).get("rcsb_entry_info") or {}).get("experimental_method") or [])
    if isinstance(vals, str):
        return vals
    if isinstance(vals, list):
        return ";".join(str(v) for v in vals if v is not None)
    return ""


def parse_assembly_id(assembly_id: str) -> Tuple[str, str]:
    pdb, asm = assembly_id.rsplit("-", 1)
    return pdb.upper(), asm


def split_csv_ids(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, list):
        out: List[str] = []
        for item in value:
            out.extend(split_csv_ids(item))
        return out
    return [x.strip() for x in str(value).split(",") if x.strip()]


def cluster_id_for_entity(polymer_entity: Dict[str, Any]) -> Optional[str]:
    memberships = polymer_entity.get("rcsb_polymer_entity_group_membership") or []
    for m in memberships:
        if not isinstance(m, dict):
            continue
        if m.get("aggregation_method") != "sequence_identity":
            continue
        try:
            cutoff = int(m.get("similarity_cutoff"))
        except (TypeError, ValueError):
            continue
        if cutoff == 100 and m.get("group_id"):
            return str(m["group_id"])
    rcsb_id = polymer_entity.get("rcsb_id")
    return f"singleton:{rcsb_id}" if rcsb_id else None


def uniprot_ids_for_entity(polymer_entity: Dict[str, Any]) -> str:
    """Return UniProt accession(s) for one RCSB polymer entity.

    The preferred source is the integrated RCSB UniProt object. The generic
    reference_sequence_identifiers field is used as a fallback. Multiple
    accessions are returned as a semicolon-separated string to keep the TSV
    unambiguous while still fitting into one cell.
    """
    uniprot_ids = set()

    # Preferred source: integrated UniProt annotations on the polymer entity.
    for uniprot in polymer_entity.get("uniprots") or []:
        if not isinstance(uniprot, dict):
            continue

        container = uniprot.get("rcsb_uniprot_container_identifiers") or {}
        value = container.get("uniprot_id") or uniprot.get("rcsb_id")

        if isinstance(value, list):
            uniprot_ids.update(str(x) for x in value if x)
        elif value:
            uniprot_ids.add(str(value))

    # Fallback source: generic reference-sequence cross references.
    if not uniprot_ids:
        identifiers = polymer_entity.get("rcsb_polymer_entity_container_identifiers") or {}
        for ref in identifiers.get("reference_sequence_identifiers") or []:
            if not isinstance(ref, dict):
                continue
            database_name = str(ref.get("database_name") or "").lower()
            if database_name == "uniprot":
                accession = ref.get("database_accession")
                if isinstance(accession, list):
                    uniprot_ids.update(str(x) for x in accession if x)
                elif accession:
                    uniprot_ids.add(str(accession))

    return ";".join(sorted(uniprot_ids))


def build_asym_to_entity(assembly: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    mapping: Dict[str, Dict[str, Any]] = {}
    for pe in (((assembly.get("entry") or {}).get("polymer_entities")) or []):
        poly_type = ((pe.get("entity_poly") or {}).get("rcsb_entity_polymer_type") or "")
        if "Protein" not in str(poly_type):
            continue
        ids = pe.get("rcsb_polymer_entity_container_identifiers") or {}
        for asym in split_csv_ids(ids.get("asym_ids")):
            mapping[asym] = pe
    return mapping


def participating_unique_protein_entities(assembly: Dict[str, Any]) -> List[Dict[str, Any]]:
    asym_to_entity = build_asym_to_entity(assembly)
    unique_by_rcsb_id: Dict[str, Dict[str, Any]] = {}

    for gen in (assembly.get("pdbx_struct_assembly_gen") or []):
        for asym_id in split_csv_ids(gen.get("asym_id_list")):
            pe = asym_to_entity.get(asym_id)
            if pe is None:
                continue
            rcsb_id = str(pe.get("rcsb_id") or "")
            if not rcsb_id:
                continue
            unique_by_rcsb_id[rcsb_id] = pe

    return [unique_by_rcsb_id[k] for k in sorted(unique_by_rcsb_id)]


def candidate_from_assembly(assembly: Dict[str, Any]) -> Optional[Candidate]:
    assembly_id = str(assembly.get("rcsb_id") or "").upper()
    if not assembly_id or "-" not in assembly_id:
        return None
    pdb_id, asm_no = parse_assembly_id(assembly_id)

    info = assembly.get("rcsb_assembly_info") or {}
    if info.get("polymer_entity_instance_count") != 2 or info.get("polymer_entity_instance_count_protein") != 2:
        return None

    entities = participating_unique_protein_entities(assembly)
    if len(entities) == 1:
        entity_id = str(entities[0].get("rcsb_id") or "")
        cluster_id = cluster_id_for_entity(entities[0])
        if not entity_id or cluster_id is None:
            print(f"Skipping {assembly_id}: missing protein-entity or cluster information", file=sys.stderr)
            return None
        uniprot_id = uniprot_ids_for_entity(entities[0])
        dimer_type = "homo"
        entity_pair = (entity_id, entity_id)
        cluster_pair = (str(cluster_id), str(cluster_id))
        uniprot_pair = (uniprot_id, uniprot_id)
    elif len(entities) == 2:
        entity_records: List[Tuple[str, str, str]] = []
        for pe in entities:
            entity_id = str(pe.get("rcsb_id") or "")
            cluster_id = cluster_id_for_entity(pe)
            if not entity_id or cluster_id is None:
                print(f"Skipping {assembly_id}: missing protein-entity or cluster information", file=sys.stderr)
                return None
            uniprot_id = uniprot_ids_for_entity(pe)
            entity_records.append((entity_id, str(cluster_id), uniprot_id))
        entity_records.sort(key=lambda item: item[0])
        dimer_type = "hetero"
        entity_pair = (entity_records[0][0], entity_records[1][0])
        cluster_pair = tuple(sorted((entity_records[0][1], entity_records[1][1])))
        uniprot_pair = (entity_records[0][2], entity_records[1][2])
    else:
        print(
            f"Skipping {assembly_id}: expected 1 or 2 unique protein entities in a protein dimer assembly, got {len(entities)}",
            file=sys.stderr,
        )
        return None

    resolution = parse_resolution(assembly)
    modeled = int(info.get("modeled_polymer_monomer_count") or 0)
    method = parse_method(assembly)
    oligomeric_count = str(((assembly.get("pdbx_struct_assembly") or {}).get("oligomeric_count")) or "")
    local_filename = f"{pdb_id.lower()}-assembly{asm_no}.cif.gz"
    download_url = f"{FILES_BASE}/{pdb_id.lower()}-assembly{asm_no}.cif.gz"

    return Candidate(
        assembly_id=assembly_id,
        pdb_id=pdb_id,
        assembly_number=asm_no,
        cluster_pair=(str(cluster_pair[0]), str(cluster_pair[1])),
        entity_pair=(str(entity_pair[0]), str(entity_pair[1])),
        uniprot_pair=(str(uniprot_pair[0]), str(uniprot_pair[1])),
        dimer_type=dimer_type,
        resolution=resolution,
        modeled_residue_count=modeled,
        method=method,
        oligomeric_count=oligomeric_count,
        local_filename=local_filename,
        download_url=download_url,
    )


def write_tsv(path: Path, rows: Sequence[Candidate]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "assembly_id",
            "pdb_id",
            "assembly_number",
            "entity_pair",
            "uniprot_pair",
            "uniprot_1",
            "uniprot_2",
            "dimer_type",
            "cluster_pair_100pct",
            "resolution_best_angstrom",
            "modeled_polymer_monomer_count",
            "experimental_method",
            "oligomeric_count",
            "download_url",
            "local_filename",
        ])
        for c in rows:
            writer.writerow([
                c.assembly_id,
                c.pdb_id,
                c.assembly_number,
                ",".join(c.entity_pair),
                "|".join(c.uniprot_pair),
                c.uniprot_pair[0],
                c.uniprot_pair[1],
                c.dimer_type,
                ",".join(c.cluster_pair),
                "" if math.isinf(c.resolution) else c.resolution,
                c.modeled_residue_count,
                c.method,
                c.oligomeric_count,
                c.download_url,
                c.local_filename,
            ])


def select_representatives(candidates: Sequence[Candidate]) -> List[Candidate]:
    best_by_pair: Dict[Tuple[str, str], Candidate] = {}
    for c in candidates:
        old = best_by_pair.get(c.cluster_pair)
        if old is None or c.rank_key() < old.rank_key():
            best_by_pair[c.cluster_pair] = c
    return sorted(best_by_pair.values(), key=lambda c: (c.pdb_id, c.assembly_number))


def main() -> int:
    parser = argparse.ArgumentParser(description="Download non-redundant RCSB protein-only biological-assembly dimers.")
    parser.add_argument("--outdir", default="rcsb_protein_dimers", help="Output directory [default: %(default)s]")
    parser.add_argument("--batch-size", type=int, default=200, help="GraphQL metadata batch size [default: %(default)s]")
    parser.add_argument("--no-ligands", action="store_true", help="Also require zero non-polymer instances in the assembly")
    parser.add_argument("--no-download", action="store_true", help="Write TSV files only; do not download coordinates")
    parser.add_argument("--limit", type=int, default=0, help="Debug limit on candidate assemblies before metadata fetch")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    assembly_dir = outdir / "assemblies"

    # Extract asssemblies
    print("Searching RCSB for candidate protein-only dimer assemblies...", file=sys.stderr)
    assembly_ids = search_candidate_assemblies(no_ligands=args.no_ligands)
    if args.limit:
        assembly_ids = assembly_ids[:args.limit]
    print(f"Search returned {len(assembly_ids)} candidate assemblies", file=sys.stderr)

    # Extract metadata for assemblies
    print("Fetching assembly/entity metadata, UniProt mappings, and 100% sequence-cluster memberships...", file=sys.stderr)
    metadata = fetch_assembly_metadata(assembly_ids, batch_size=args.batch_size)

    # Extract wanted metadata fields (candidates) and select representatives (reps)
    candidates: List[Candidate] = []
    for assembly in metadata:
        c = candidate_from_assembly(assembly)
        if c is not None:
            candidates.append(c)
    candidates = sorted(candidates, key=lambda c: (c.pdb_id, c.assembly_number))
    reps = select_representatives(candidates)
 
    # Write TSV files for all candidates and selected representatives
    write_tsv(outdir / "all_candidate_assemblies.tsv", candidates)
    write_tsv(outdir / "selected_assemblies.tsv", reps)

    print(f"Parsed {len(candidates)} candidate assemblies", file=sys.stderr)
    print(f"Selected {len(reps)} non-redundant representatives at 100% sequence-pair identity", file=sys.stderr)

    # Download coordinate files for representatives, skipping existing files
    if not args.no_download:
        for i, c in enumerate(reps, start=1):
            dest = assembly_dir / c.local_filename
            if dest.exists() and dest.stat().st_size > 0:
                print(f"[{i}/{len(reps)}] exists {dest}", file=sys.stderr)
                continue
            print(f"[{i}/{len(reps)}] downloading {c.assembly_id}", file=sys.stderr)
            download_binary(c.download_url, dest)

    print(f"Done. Representatives: {outdir / 'selected_assemblies.tsv'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
