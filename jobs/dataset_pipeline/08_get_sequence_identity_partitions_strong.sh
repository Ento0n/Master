#!/bin/bash
#
#SBATCH --job-name=08_get_seq_ident_partitions_strong                     # Job name
#SBATCH --output=logs/08_get_seq_ident_partitions_strong.%j.out           # Standard output (%j expands to jobId)
#SBATCH --error=logs/08_get_seq_ident_partitions_strong.%j.err            # Standard error
#SBATCH --partition=shared-cpu                      # not standard otherwise no permissions
##SBATCH --nodelist=grover.exbio.wzw.tum.de          # on compms I get illegal instruction error
#SBATCH --ntasks=1                                  # Run a single task
#SBATCH --cpus-per-task=1                           # Number of CPU cores per task
#SBATCH --mem=96G                                   # Total memory
#SBATCH --time=4-00:00:00                           # Time limit hh:mm:ss
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=antspa@gmx.de                   # Where to send mail

/nfs/home/students/a.spannagl/.local/kahip64/bin/kaffpa \
  /nfs/scratch/pdb_dimers/graphs/sequence_identity.graph \
  --k=10 \
  --preconfiguration=strong \
  --imbalance=3 \
  --output_filename=/nfs/scratch/pdb_dimers/KaHIP/seq_ident_partitions_strong_output.txt