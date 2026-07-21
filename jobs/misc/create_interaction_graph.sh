#!/bin/bash
#
#SBATCH --job-name=create_interaction_graph                     # Job name
#SBATCH --output=logs/create_interaction_graph.%j.out           # Standard output (%j expands to jobId)
#SBATCH --error=logs/create_interaction_graph.%j.err            # Standard error
#SBATCH --partition=shared-cpu                      # not standard otherwise no permissions
#SBATCH --ntasks=1                                  # Run a single task
#SBATCH --cpus-per-task=1                           # Number of CPU cores per task
#SBATCH --mem=128G                                   # Total memory
#SBATCH --time=4-00:00:00                           # Time limit hh:mm:ss
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=antspa@gmx.de                   # Where to send mail

python ../create_graph.py --type interaction