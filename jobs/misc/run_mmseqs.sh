#!/bin/bash
#
#SBATCH --job-name=run_mmseqs                       # Job name
#SBATCH --output=logs/run_mmseqs.%j.out             # Standard output (%j expands to jobId)
#SBATCH --error=logs/run_mmseqs.%j.err              # Standard error
#SBATCH --partition=shared-cpu                      # not standard otherwise no permissions
#SBATCH --ntasks=1                                  # Run a single task
#SBATCH --cpus-per-task=1                           # Number of CPU cores per task
#SBATCH --mem=64G                                   # Total memory
#SBATCH --time=4-00:00:00                           # Time limit hh:mm:ss
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=antspa@gmx.de                   # Where to send mail

mmseqs easy-search /nfs/scratch/pdb_dimers/unique_sequences.fasta /nfs/scratch/pdb_dimers/unique_sequences.fasta /nfs/scratch/pdb_dimers/MMseqs2/hits.tsv /nfs/scratch/pdb_dimers/MMseqs2/tmp \
  --format-output query,target,pident,nident,alnlen,qcov,tcov,evalue,bits \
  --min-seq-id 0.0 \
  -c 0.0 \
  --cov-mode 0 \
  -s 7.5 \
  --max-seqs 100000