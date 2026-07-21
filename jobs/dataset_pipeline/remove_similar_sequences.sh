#!/bin/bash
#
#SBATCH --job-name=remove_similar_sequences                     # Job name
#SBATCH --output=logs/remove_similar_sequences.%j.out           # Standard output (%j expands to jobId)
#SBATCH --error=logs/remove_similar_sequences.%j.err            # Standard error
#SBATCH --partition=shared-cpu                      # not standard otherwise no permissions
#SBATCH --ntasks=1                                  # Run a single task
#SBATCH --cpus-per-task=1                           # Number of CPU cores per task
#SBATCH --mem=16G                                   # Total memory
#SBATCH --time=4-00:00:00                           # Time limit hh:mm:ss
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=antspa@gmx.de                   # Where to send mail

python ../remove_similar_sequences.py