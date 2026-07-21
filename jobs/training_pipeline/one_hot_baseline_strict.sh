#!/bin/bash
#
#SBATCH --job-name=one_hot_baseline_strict                     # Job name
#SBATCH --output=logs/one_hot_baseline_strict.%j.out           # Standard output (%j expands to jobId)
#SBATCH --error=logs/one_hot_baseline_strict.%j.err            # Standard error
#SBATCH --partition=shared-gpu                      # not standard otherwise no permissions
#SBATCH --gres=gpu:1
#SBATCH --qos=limitgpus
#SBATCH --exclude=gpu01.exbio.wzw.tum.de
#SBATCH --ntasks=1                                  # Run a single task
#SBATCH --cpus-per-task=4                           # Number of CPU cores per task
#SBATCH --mem=64G                                   # Total memory
#SBATCH --time=4-00:00:00                           # Time limit hh:mm:ss
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=antspa@gmx.de                   # Where to send mail

python ../execute_pipeline.py --interactions /nfs/scratch/pdb_dimers/balanced_interactions_strict.tsv --sequence-file /nfs/scratch/pdb_dimers/entity_sequences.tsv --output-dir /nfs/scratch/pdb_dimers/ML_runs/one_hot_baseline_strict --accelerator gpu