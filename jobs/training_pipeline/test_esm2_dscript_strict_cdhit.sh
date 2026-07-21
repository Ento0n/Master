#!/bin/bash
#
#SBATCH --job-name=esm2_dscript_strict_cdhit                     # Job name
#SBATCH --output=logs/esm2_dscript_strict_cdhit.%j.out           # Standard output (%j expands to jobId)
#SBATCH --error=logs/esm2_dscript_strict_cdhit.%j.err            # Standard error
#SBATCH --partition=shared-gpu                      # not standard otherwise no permissions
#SBATCH --gres=gpu:1
#SBATCH --qos=limitgpus
#SBATCH --exclude=gpu01.exbio.wzw.tum.de,compms-gpu-1.exbio.wzw.tum.de
#SBATCH --ntasks=1                                  # Run a single task
#SBATCH --cpus-per-task=4                           # Number of CPU cores per task
#SBATCH --mem=64G                                   # Total memory
#SBATCH --time=4-00:00:00                           # Time limit hh:mm:ss
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=antspa@gmx.de                   # Where to send mail

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <string>" >&2
    exit 1
fi

argument="$1"

python ../execute_pipeline.py \
  --model dscript \
  --interactions /nfs/scratch/pdb_dimers/balanced_interactions_strict_cd_hit.tsv \
  --interaction-loss-lambda 0.95 \
  --contact-threshold 0.5 \
  --batch-size 8 \
  --max-epochs 1 \
  --output-subdir "$argument"