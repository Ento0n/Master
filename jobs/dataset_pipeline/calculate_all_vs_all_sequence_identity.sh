#!/bin/bash
#
#SBATCH --job-name=calculate_all_vs_all_sequence_identity                     # Job name
#SBATCH --output=logs/calculate_all_vs_all_sequence_identity.%j.out           # Standard output (%j expands to jobId)
#SBATCH --error=logs/calculate_all_vs_all_sequence_identity.%j.err            # Standard error
#SBATCH --partition=shared-cpu                      # not standard otherwise no permissions
#SBATCH --ntasks=1                                  # Run a single task
#SBATCH --cpus-per-task=1                           # Number of CPU cores per task
#SBATCH --mem=512G                                   # Total memory
#SBATCH --time=4-00:00:00                           # Time limit hh:mm:ss
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=antspa@gmx.de                   # Where to send mail

python ../calculate_all_vs_all_sequence_identity.py --local
