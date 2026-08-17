#!/bin/bash
#
#SBATCH --job-name=13_sample_negatives_keep_homomers_cd_hit                     # Job name
#SBATCH --output=logs/13_sample_negatives_keep_homomers_cd_hit.%j.out           # Standard output (%j expands to jobId)
#SBATCH --error=logs/13_sample_negatives_keep_homomers_cd_hit.%j.err            # Standard error
#SBATCH --partition=shared-cpu                      # not standard otherwise no permissions
#SBATCH --ntasks=1                                  # Run a single task
#SBATCH --cpus-per-task=1                           # Number of CPU cores per task
#SBATCH --mem=16G                                   # Total memory
#SBATCH --time=4-00:00:00                           # Time limit hh:mm:ss
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=antspa@gmx.de                   # Where to send mail

python /nfs/home/students/a.spannagl/master_repository/scripts/data_pipeline/13_sample_negatives.py --input /nfs/scratch/pdb_dimers/final_final_filtered_interactions_with_partitions_cd_hit.tsv --output /nfs/scratch/pdb_dimers/balanced_interactions_keep_homomers_cd_hit.tsv --keep-positive-homomers