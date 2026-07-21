#!/bin/bash
#
#SBATCH --job-name=cd_hit_train_val                     # Job name
#SBATCH --output=logs/cd_hit_train_val.%j.out           # Standard output (%j expands to jobId)
#SBATCH --error=logs/cd_hit_train_val.%j.err            # Standard error
#SBATCH --partition=shared-cpu                      # not standard otherwise no permissions
#SBATCH --ntasks=1                                  # Run a single task
#SBATCH --cpus-per-task=1                           # Number of CPU cores per task
#SBATCH --mem=16G                                   # Total memory
#SBATCH --time=4-00:00:00                           # Time limit hh:mm:ss
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=antspa@gmx.de                   # Where to send mail

cd-hit-2d -i /nfs/scratch/pdb_dimers/CD-HIT/train.fasta -i2 /nfs/scratch/pdb_dimers/CD-HIT/val.fasta -o /nfs/scratch/pdb_dimers/CD-HIT/train_val.out -n 2 -c 0.4 -s2 0 -S2 999999