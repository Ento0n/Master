#!/bin/bash
#
#SBATCH --job-name=create_benchmark_files
#SBATCH --output=/nfs/home/students/a.spannagl/master_repository/data/slurm/logs/create_benchmark_files.%j.out
#SBATCH --error=/nfs/home/students/a.spannagl/master_repository/data/slurm/logs/create_benchmark_files.%j.err
#SBATCH --partition=shared-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=antspa@gmx.de

set -euo pipefail

REPO_DIR=/nfs/home/students/a.spannagl/master_repository

mkdir -p "$REPO_DIR/data/slurm/logs"

python "$REPO_DIR/data/create_benchmark_files.py" \
  --interactions-tsv /nfs/scratch/pdb_dimers/final_final_filtered_interactions_with_partitions_cd_hit.tsv \
  --entity-sequences-tsv /nfs/scratch/pdb_dimers/entity_sequences.tsv \
  --outdir /nfs/scratch/pdb_dimers/ppi_splitting_pipeline \
  --batch-size 100 \
  --sleep 0.2
