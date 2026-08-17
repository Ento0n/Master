#!/bin/bash
#
#SBATCH --job-name=create_benchmark_files
#SBATCH --output=logs/create_benchmark_files.%j.out
#SBATCH --error=logs/create_benchmark_files.%j.err
#SBATCH --partition=shared-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=antspa@gmx.de

python /nfs/home/students/a.spannagl/master_repository/scripts/misc/create_benchmark_files.py \
  --interactions-tsv /nfs/scratch/pdb_dimers/dataset_iterations/04_removed_duplicates.tsv \
  --entity-sequences-tsv /nfs/scratch/pdb_dimers/sequences/entity_sequences.tsv \
  --outdir /nfs/scratch/pdb_dimers/ppi_splitting_pipeline \
  --batch-size 100 \
  --sleep 0.2
