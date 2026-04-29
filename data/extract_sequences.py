import sys
import os

from numpy import block
import pandas as pd
import gemmi

DATA_DIR = "/nfs/scratch/pdb_dimers/"

def main():
    # Load the selected assemblies
    selected_assemblies = pd.read_csv(os.path.join(DATA_DIR, "selected_assemblies.tsv"), sep="\t")

    # Collect all unique entity names from the entity_pair column and the corresponding cluster ID
    entity_names = list()
    cluster_ids = list()
    sequences = list()
    missing_sequences = list()
    binary_masks = list()
    for _, row in selected_assemblies.iterrows():
        # entity
        entity_pair = row["entity_pair"]
        entities = entity_pair.split(",")
        for i, entity in enumerate(entities):
            if entity.strip() in entity_names:
                continue
            
            entity_names.append(entity.strip())

            # cluster
            cluster_pair = row["cluster_pair_100pct"]
            cluster_ids_pair = cluster_pair.split(",")
            if i == 0:
                cluster_ids.append(cluster_ids_pair[0].strip())
            else:
                cluster_ids.append(cluster_ids_pair[1].strip())
            
            # Check if CIF file exists
            CIF_PATH = os.path.join(DATA_DIR, "assemblies", f"{row['pdb_id'].lower()}-assembly{row['assembly_number']}.cif.gz")
            if not os.path.exists(CIF_PATH):
                print(f"Warning: CIF file {CIF_PATH} does not exist for entity {entity} in assembly {row['assembly_number']} of PDB ID {row['pdb_id']}", file=sys.stderr)
                sequences.append("-") # Placeholder for missing sequence
                missing_sequences.append(entity)
                continue

            # read file and set up entities
            st = gemmi.read_structure(CIF_PATH)
            st.setup_entities()
            entity_number = entity.split("_")[1].strip()

            # extract remapping information
            doc = gemmi.cif.read_file(CIF_PATH)
            block = doc.sole_block()

            cat = block.find_mmcif_category("_pdbx_entity_remapping.")

            for remapping_row in cat:
                if remapping_row[1] == entity_number and remapping_row[0] != entity_number:
                    remapping = {
                        "original_entity": remapping_row[1],
                        "remapped_entity": remapping_row[0],
                    }
                    print(f"Debug: Found remapping for entity {entity} in file {CIF_PATH}", file=sys.stderr)
                    print(f"Debug: Remapping details: {remapping}", file=sys.stderr)
                    entity_number = remapping_row[0] # Update entity number based on remapping
                    break
            
            # asym_id mapping
            asym_ids = list()
            table = block.find("_struct_asym.", ["id", "entity_id"])
            for asym_row in table:
                if entity_number == asym_row[1]:
                    asym_ids.append(asym_row[0])

            # extract sequence information
            size_before = len(sequences)
            seq_length = 0
            for st_entity in st.entities:
                if st_entity.name == entity_number:
                    seq_3letter = st_entity.full_sequence

                    if entity == "1BCK_2":
                        print(f"Debug: Processing entity {entity} in file {CIF_PATH}", file=sys.stderr)
                        print(f"Debug: Entity number extracted: {entity_number}", file=sys.stderr)
                        print(f"Debug: Number of entities in structure: {len(st.entities)}", file=sys.stderr)
                        print(f"Debug: Entity names in structure: {[e.name for e in st.entities]}", file=sys.stderr)
                        print(f"Debug: 3-letter sequence: {seq_3letter}", file=sys.stderr)

                    seq_1letter = gemmi.one_letter_code(seq_3letter)
                    seq_length = len(seq_1letter)
                    sequences.append(seq_1letter)

            if len(sequences) == size_before:
                print(f"Warning: No sequence found for entity {entity} in file {CIF_PATH}", file=sys.stderr)
                sequences.append("-") # Placeholder for missing sequence
                missing_sequences.append(entity)
            
            # Structure annotation mask per residue
            asym_id_to_residues_mask = {asym_id: [] for asym_id in asym_ids}
            for atom_row in block.find(
                "_atom_site.",
                [
                    "label_asym_id",
                    "label_seq_id",
                ],
            ):
                asym_id, seq_id = atom_row

                if asym_id in asym_ids:
                    if int(seq_id)-1 in asym_id_to_residues_mask[asym_id]:
                        continue  # Skip duplicate residue entries
                    else:
                        asym_id_to_residues_mask[asym_id].append(int(seq_id)-1)  # Store the index of present residues
        
            # convert number lists to binary masks
            for asym_id in asym_id_to_residues_mask:
                mask = [0] * seq_length  # Initialize mask with all residues marked as missing
                for index in asym_id_to_residues_mask[asym_id]:
                    if index < seq_length:
                        mask[index] = 1  # Mark present residues as True
                    else:
                        print(f"Warning: Residue index {index} exceeds sequence length {seq_length} for entity {entity} in file {CIF_PATH}. Asym_id: {asym_id}", file=sys.stderr)
                asym_id_to_residues_mask[asym_id] = mask
            
            # convert dict to list of masks (one per sequence)
            entry = ""
            for asym_id in asym_ids:
                entry += f"{asym_id}: {"".join(map(str, asym_id_to_residues_mask[asym_id]))}; "
            binary_masks.append(entry.strip())


    
    # print different lengths
    print(f"Number of unique entities: {len(entity_names)}")
    print(f"Number of cluster IDs: {len(cluster_ids)}")
    print(f"Number of sequences: {len(sequences)}")
    print(f"Number of missing sequences: {len(missing_sequences)}")
    print(f"Number of binary masks: {len(binary_masks)}")
    if missing_sequences:
        print("Missing sequences for entities:", missing_sequences, file=sys.stderr)

    # Create a DataFrame to store the entity names, cluster IDs, and sequences
    df = pd.DataFrame({
        "entity_name": list(entity_names),
        "cluster_id": cluster_ids,
        "sequence": sequences,
        "binary_mask": binary_masks
    })

    # Save the DataFrame to a TSV file
    df.to_csv(os.path.join(DATA_DIR, "entity_sequences.tsv"), sep="\t", index=False)
                        


if __name__ == "__main__":
    main()