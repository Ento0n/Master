#!/bin/bash
#
#SBATCH --job-name=03_filter_by_seq_att                     # Job name
#SBATCH --output=logs/03_filter_by_seq_att.%j.out           # Standard output (%j expands to jobId)
#SBATCH --error=logs/03_filter_by_seq_att.%j.err            # Standard error
#SBATCH --partition=shared-cpu                      # not standard otherwise no permissions
#SBATCH --ntasks=1                                  # Run a single task
#SBATCH --cpus-per-task=1                           # Number of CPU cores per task
#SBATCH --mem=16G                                   # Total memory
#SBATCH --time=4-00:00:00                           # Time limit hh:mm:ss
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=antspa@gmx.de                   # Where to send mail

python /nfs/home/students/a.spannagl/master_repository/scripts/data_pipeline/03_filter_by_seq_att.py