#!/usr/bin/env python3
"""Create interaction, sequence, GO annotation, and species files."""

from __future__ import annotations

import argparse
import csv
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("--interactions-tsv", default="/nfs/scratch/pdb_dimers/final_filtered_interactions.tsv")
parser.add_argument("--entity-sequences-tsv", default="/nfs/scratch/pdb_dimers/entity_sequences.tsv")
parser.add_argument("--outdir", default="/nfs/scratch/pdb_dimers/ppi_splitting_pipeline")
parser.add_argument("--batch-size", type=int, default=100)
parser.add_argument("--sleep", type=float, default=0.2)
parser.add_argument("--dry-run", action="store_true")
args = parser.parse_args()

outdir = Path(args.outdir)
splitter = re.compile(r"[;|,\s]+")
missing_values = {"", "-", "na", "n/a", "nan", "none", "null"}

print(f"[{datetime.now().isoformat(timespec='seconds')}] Reading interactions")

interactions = []
entity_to_uniprots = defaultdict(Counter)
entity_to_taxa = defaultdict(Counter)

with open(args.interactions_tsv, newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    needed = {"entity_pair", "uniprot_1", "uniprot_2", "taxonomy_1", "taxonomy_2"}
    if not reader.fieldnames or needed - set(reader.fieldnames):
        sys.exit(f"{args.interactions_tsv} is missing columns: {', '.join(sorted(needed - set(reader.fieldnames or [])))}")

    for row in reader:
        pdbid_1, pdbid_2 = [part.strip() for part in row["entity_pair"].split(",")]
        protein1 = next((x.strip().upper() for x in splitter.split(row["uniprot_1"]) if x.strip().lower() not in missing_values), "")
        protein2 = next((x.strip().upper() for x in splitter.split(row["uniprot_2"]) if x.strip().lower() not in missing_values), "")
        interactions.append((pdbid_1, pdbid_2, protein1, protein2))

        if protein1:
            entity_to_uniprots[pdbid_1][protein1.upper()] += 1
        if protein2:
            entity_to_uniprots[pdbid_2][protein2.upper()] += 1
        if row["taxonomy_1"].strip():
            entity_to_taxa[pdbid_1][row["taxonomy_1"].strip()] += 1
        if row["taxonomy_2"].strip():
            entity_to_taxa[pdbid_2][row["taxonomy_2"].strip()] += 1

protein_ids = sorted({item for row in interactions for item in row[:2]})
entity_to_uniprot = {entity: counts.most_common(1)[0][0] for entity, counts in entity_to_uniprots.items()}
entity_to_taxon = {entity: counts.most_common(1)[0][0] for entity, counts in entity_to_taxa.items()}
uniprot_ids = sorted(set(entity_to_uniprot.values()))

print(f"[{datetime.now().isoformat(timespec='seconds')}] Interactions: {len(interactions)}")
print(f"[{datetime.now().isoformat(timespec='seconds')}] Protein IDs: {len(protein_ids)}")
print(f"[{datetime.now().isoformat(timespec='seconds')}] UniProt IDs for GO lookup: {len(uniprot_ids)}")

if args.dry_run:
    sys.exit(0)

outdir.mkdir(parents=True, exist_ok=True)
tmp_interactions = outdir / "pdb_interactions.csv.tmp"
tmp_sequences = outdir / "sequences.fasta.tmp"
tmp_go = outdir / "go_annotations.tsv.tmp"
tmp_species = outdir / "species.tsv.tmp"
tmp_missing = outdir / "missing_uniprot_go_ids.txt.tmp"

with tmp_interactions.open("w", newline="") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(["pdbid_1", "pdbid_2", "protein1", "protein2"])
    writer.writerows(interactions)

print(f"[{datetime.now().isoformat(timespec='seconds')}] Reading sequences")
sequences = {}
with open(args.entity_sequences_tsv, newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        if row["entity_name"] in protein_ids:
            sequences[row["entity_name"]] = row["sequence"]

with tmp_sequences.open("w") as handle:
    for protein_id in protein_ids:
        if protein_id not in sequences:
            print(f"Missing sequence for {protein_id}", file=sys.stderr)
            continue
        handle.write(f">{protein_id}\n{sequences[protein_id]}\n")

go_by_uniprot = {}
missing_uniprot_ids = []
fields = "accession,go_p,go_f,go_c,organism_id"
total_batches = (len(uniprot_ids) + args.batch_size - 1) // args.batch_size

for start in range(0, len(uniprot_ids), args.batch_size):
    batch = uniprot_ids[start : start + args.batch_size]
    print(f"[{datetime.now().isoformat(timespec='seconds')}] UniProt GO batch {start // args.batch_size + 1}/{total_batches}")

    query = urllib.parse.urlencode({"accessions": ",".join(batch), "fields": fields, "format": "tsv"})
    url = f"https://rest.uniprot.org/uniprotkb/accessions?{query}"
    try:
        with urllib.request.urlopen(url, timeout=120) as response:
            text = response.read().decode("utf-8")
    except urllib.error.HTTPError:
        text = "Entry\tGene Ontology (biological process)\tGene Ontology (molecular function)\tGene Ontology (cellular component)\tOrganism (ID)\n"
        for uniprot_id in batch:
            query = urllib.parse.urlencode({"accessions": uniprot_id, "fields": fields, "format": "tsv"})
            url = f"https://rest.uniprot.org/uniprotkb/accessions?{query}"
            try:
                with urllib.request.urlopen(url, timeout=120) as response:
                    lines = response.read().decode("utf-8").splitlines()
                    text += "\n".join(lines[1:]) + "\n"
            except urllib.error.HTTPError:
                missing_uniprot_ids.append(uniprot_id)

    for row in csv.DictReader(text.splitlines(), delimiter="\t"):
        go_by_uniprot[row["Entry"]] = (
            row["Gene Ontology (biological process)"],
            row["Gene Ontology (molecular function)"],
            row["Gene Ontology (cellular component)"],
        )

    if args.sleep and start + args.batch_size < len(uniprot_ids):
        time.sleep(args.sleep)

with tmp_go.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["protein_id", "go_bp", "go_mf", "go_cc"])
    for protein_id in protein_ids:
        writer.writerow([protein_id, *go_by_uniprot.get(entity_to_uniprot.get(protein_id, ""), ("", "", ""))])

with tmp_species.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["protein_id", "taxon_id"])
    for protein_id in protein_ids:
        writer.writerow([protein_id, entity_to_taxon.get(protein_id, "")])

with tmp_missing.open("w") as handle:
    for uniprot_id in sorted(missing_uniprot_ids):
        handle.write(f"{uniprot_id}\n")

tmp_interactions.replace(outdir / "pdb_interactions.csv")
tmp_sequences.replace(outdir / "sequences.fasta")
tmp_go.replace(outdir / "go_annotations.tsv")
tmp_species.replace(outdir / "species.tsv")
tmp_missing.replace(outdir / "missing_uniprot_go_ids.txt")

print(f"[{datetime.now().isoformat(timespec='seconds')}] Wrote {outdir / 'pdb_interactions.csv'}")
print(f"[{datetime.now().isoformat(timespec='seconds')}] Wrote {outdir / 'sequences.fasta'}")
print(f"[{datetime.now().isoformat(timespec='seconds')}] Wrote {outdir / 'go_annotations.tsv'}")
print(f"[{datetime.now().isoformat(timespec='seconds')}] Wrote {outdir / 'species.tsv'}")
print(f"[{datetime.now().isoformat(timespec='seconds')}] Missing UniProt GO IDs: {len(missing_uniprot_ids)}")
