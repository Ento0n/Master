import sys
import os

import pandas as pd
import gemmi

DATA_DIR = "/nfs/scratch/pdb_dimers/"

def main():
    # Load the selected assemblies
    selected_assemblies = pd.read_csv(os.path.join(DATA_DIR, "selected_assemblies.tsv"), sep="\t")

    # Collect all unique entity names from the entity_pair column and the corresponding cluster ID
    entity_names = set()
    cluster_ids = list()
    sequences = list()
    missing_sequences = list()
    for _, row in selected_assemblies.iterrows():
        # entity
        entity_pair = row["entity_pair"]
        entities = entity_pair.split(",")
        for i, entity in enumerate(entities):
            if entity.strip() in entity_names:
                continue
            
            entity_names.add(entity.strip())

            # cluster
            cluster_pair = row["cluster_pair_100pct"]
            cluster_ids_pair = cluster_pair.split(",")
            if i == 0:
                cluster_ids.append(cluster_ids_pair[0].strip())
            else:
                cluster_ids.append(cluster_ids_pair[1].strip())

            # extract sequence for entity
            cif_file = os.path.join(DATA_DIR, "assemblies", f"{row['pdb_id'].lower()}-assembly{row['assembly_number']}.cif.gz")
            st = gemmi.read_structure(cif_file)
            st.setup_entities()
            entity_number = entity.split("_")[1].strip()

            if entity == "1MKX_3":
                print(f"Debug: Processing entity {entity} in file {cif_file}", file=sys.stderr)
                print(f"Debug: Entity number extracted: {entity_number}", file=sys.stderr)
                print(f"Debug: Number of entities in structure: {len(st.entities)}", file=sys.stderr)
                print(f"Debug: Entity names in structure: {[e.name for e in st.entities]}", file=sys.stderr)

            size_before = len(sequences)
            for st_entity in st.entities:
                if st_entity.name == entity_number:
                    seq_3letter = st_entity.full_sequence
                    seq_1letter = gemmi.one_letter_code(seq_3letter)
                    sequences.append(seq_1letter)
                    print(f"Entity: {entity}, Cluster ID: {cluster_ids[-1]}, Sequence: {seq_1letter}")

            if len(sequences) == size_before:
                print(f"Warning: No sequence found for entity {entity} in file {cif_file}", file=sys.stderr)
                sequences.append("-") # Placeholder for missing sequence
                missing_sequences.append(entity)
    
    # print different lengths
    print(f"Number of unique entities: {len(entity_names)}")
    print(f"Number of cluster IDs: {len(cluster_ids)}")
    print(f"Number of sequences: {len(sequences)}")
    print(f"Number of missing sequences: {len(missing_sequences)}")
    if missing_sequences:
        print("Missing sequences for entities:", missing_sequences, file=sys.stderr)

    # Create a DataFrame to store the entity names, cluster IDs, and sequences
    df = pd.DataFrame({
        "entity_name": list(entity_names),
        "cluster_id": cluster_ids,
        "sequence": sequences
    })

    # Save the DataFrame to a TSV file
    df.to_csv(os.path.join(DATA_DIR, "entity_sequences.tsv"), sep="\t", index=False)
                        


if __name__ == "__main__":
    main()