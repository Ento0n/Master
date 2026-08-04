#!/bin/bash
#
#SBATCH --job-name=esm2_query_patch
#SBATCH --output=/nfs/home/students/a.spannagl/master_repository/jobs/training_pipeline/logs/esm2_query_patch.%j.out
#SBATCH --error=/nfs/home/students/a.spannagl/master_repository/jobs/training_pipeline/logs/esm2_query_patch.%j.err
#SBATCH --partition=shared-gpu
#SBATCH --gres=gpu:1
#SBATCH --qos=limitgpus
#SBATCH --exclude=gpu01.exbio.wzw.tum.de,compms-gpu-1.exbio.wzw.tum.de
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=antspa@gmx.de

set -euo pipefail

# An optional first argument overrides the automatically generated result
# subdirectory, for example: sbatch esm2_query_patch.sh my_debug_run
output_subdir_args=()
if [[ $# -ge 1 ]]; then
    output_subdir_args=(--output-subdir "$1")
fi

cd /nfs/home/students/a.spannagl/master_repository

# These loss/head settings match the current ESM2 D-SCRIPT run so the contact
# generator is the main architectural difference between the experiments.
python -m scripts.training_pipeline.execute_pipeline \
    --model query_patch \
    --num-interaction-queries 32 \
    --query-heads 4 \
    --query-layers 1 \
    --query-dropout 0.1 \
    --query-contact-bias-init -6.0 \
    --interaction-loss-lambda 0.5 \
    --loss-type sparsity_21 \
    --int-mod-type max \
    --max-epochs 89 \
    "${output_subdir_args[@]}"
