#!/bin/bash
#
#SBATCH --job-name=07_create_sequence_identity_graph                     # Job name
#SBATCH --output=logs/07_create_sequence_identity_graph.%j.out           # Standard output (%j expands to jobId)
#SBATCH --error=logs/07_create_sequence_identity_graph.%j.err            # Standard error
#SBATCH --partition=shared-cpu                      # not standard otherwise no permissions
#SBATCH --ntasks=1                                  # Run a single task
#SBATCH --cpus-per-task=1                           # Number of CPU cores per task
#SBATCH --mem=256G                                   # Total memory
#SBATCH --time=4-00:00:00                           # Time limit dd-hh:mm:ss
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=antspa@gmx.de                   # Where to send mail

python /nfs/home/students/a.spannagl/master_repository/scripts/data_pipeline/07_create_graph.py --type sequence_identity