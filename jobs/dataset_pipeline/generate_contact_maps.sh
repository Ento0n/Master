#!/bin/bash
#
#SBATCH --job-name=generate_contact_maps                     # Job name
#SBATCH --output=logs/generate_contact_maps.%j.out           # Standard output (%j expands to jobId)
#SBATCH --error=logs/generate_contact_maps.%j.err            # Standard error
#SBATCH --partition=shared-cpu                      # not standard otherwise no permissions
#SBATCH --ntasks=1                                  # Run a single task
#SBATCH --cpus-per-task=1                           # Number of CPU cores per task
#SBATCH --mem=64G                                   # Total memory
#SBATCH --time=4-00:00:00                           # Time limit hh:mm:ss
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=antspa@gmx.de                   # Where to send mail

python ../generate_contact_maps.py \
    --selection-tsv /nfs/scratch/pdb_dimers/final_filtered_interactions_with_partitions.tsv \
    --input-dir /nfs/scratch/pdb_dimers/assemblies \
    --output-dir /nfs/scratch/pdb_dimers/contact_maps
