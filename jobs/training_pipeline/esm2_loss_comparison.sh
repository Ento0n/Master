#!/bin/bash
#
# Controlled loss selection for the representative QP-Bilin-IC architecture.
# Each adjacent pair uses the same random seed and differs only in the main
# contact-map loss, making the comparison paired across three seeds.
# ESM2 embeddings keep this decision aligned with the historical sweeps; the
# selected loss can then be frozen for the later ProstT5 model comparison.
#
# Primary selection metric: validation interaction AUPRC (higher is better).
# Exact tie-breaker: macro validation contact AUPRC from the same checkpoint.
# Test evaluation is deliberately disabled until the loss has been selected.
#
#SBATCH --job-name=esm2-loss-comparison
#SBATCH --array=0-5%4
#SBATCH --output=/nfs/home/students/a.spannagl/master_repository/jobs/training_pipeline/logs/esm2-loss-comparison-%A_%a.out
#SBATCH --error=/nfs/home/students/a.spannagl/master_repository/jobs/training_pipeline/logs/esm2-loss-comparison-%A_%a.err
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

# Submit the complete paired comparison from the repository root:
#   sbatch jobs/training_pipeline/esm2_loss_comparison.sh
# Run a two-epoch pipeline check without changing the controlled factors:
#   sbatch jobs/training_pipeline/esm2_loss_comparison.sh --max-epochs 2
# Preview one array task without starting training:
#   SLURM_ARRAY_TASK_ID=0 bash jobs/training_pipeline/esm2_loss_comparison.sh --dry-run
dry_run=false
if [[ "${1:-}" == "--dry-run" ]]; then
    dry_run=true
    shift
fi
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" && "$dry_run" == false ]]; then
    echo "Submit this script with sbatch, or use --dry-run to preview a command." >&2
    exit 1
fi

task_id="${SLURM_ARRAY_TASK_ID:-0}"
if [[ ! "$task_id" =~ ^[0-5]$ ]]; then
    echo "SLURM_ARRAY_TASK_ID must be between 0 and 5." >&2
    exit 1
fi

# Tasks 0/1, 2/3, and 4/5 form paired comparisons at seeds 42, 43, and
# 44. Within each pair, the even task uses balanced BCE and the odd task adds
# Dice. All data, architecture, optimization, and stopping settings stay fixed.
losses=(balanced_bce balanced_bce_dice)
seeds=(42 43 44)
loss_type="${losses[task_id % 2]}"
seed="${seeds[task_id / 2]}"

repo_root=/nfs/home/students/a.spannagl/master_repository
interaction_file=/nfs/scratch/pdb_dimers/final_datasets/balanced_interactions_strict_cd_hit.tsv
embedding_manifest="${ESM2_MANIFEST:-/nfs/scratch/pdb_dimers/embeddings/esm2_t33_650M/manifest.tsv}"
python_bin="${PYTHON_BIN:-python}"
run_group="esm2-loss-comparison-${SLURM_ARRAY_JOB_ID:-dry-run}"
run_name="${run_group}-qp-bilin-ic-${loss_type}-seed-${seed}"
run_output="$repo_root/jobs/training_pipeline/results/esm2_loss_comparison/${SLURM_ARRAY_JOB_ID:-dry-run}/${loss_type}_seed_${seed}"

cd "$repo_root"
if [[ ! -r "$interaction_file" ]]; then
    echo "Cannot read interaction dataset: $interaction_file" >&2
    exit 1
fi
if [[ ! -r "$embedding_manifest" ]]; then
    echo "Cannot read ESM2 embedding manifest: $embedding_manifest" >&2
    exit 1
fi

export WANDB_RUN_GROUP="${WANDB_RUN_GROUP:-$run_group}"
export WANDB_NAME="${WANDB_NAME:-$run_name}"

# Extra arguments are inserted before the controlled experiment identifiers.
# This permits short checks such as --max-epochs 2 while preventing an
# accidental override of the loss, seed, model, selection rule, or output path.
command=(
    srun "$python_bin" -m scripts.training_pipeline.execute_pipeline
    --interactions "$interaction_file"
    --embedding-manifest "$embedding_manifest"
    --embedding-dim 1280
    --interaction-loss-lambda 0.5
    --omega 0.5
    --max-epochs 100
    --batch-size 8
    --num-workers 4
    --learning-rate 0.001
    --early-stopping-patience 8
    --early-stopping-min-delta 0.001
    --num-interaction-queries 4
    --query-heads 4
    --query-layers 1
    --query-dropout 0.1
    --query-contact-bias-init -6.0
    --query-gate-bias-init 0.0
    --query-interface-loss-weight 0.2
    --interface-omega 0.2
    --compatibility-rank 64
    --compatibility-scale 10
    --accelerator gpu
    --devices 1
    "$@"
    --model query_patch
    --query-evidence-mode gated
    --compatibility-head bilinear
    --int-mod-type normal
    --max-pooling
    --trainable-final-slope
    --loss-type "$loss_type"
    --seed "$seed"
    --selection-monitor val_interaction_auprc
    --selection-tiebreaker val_contact_auprc_macro
    --skip-test
    --output-dir "$run_output"
    --output-subdir qp_bilin_ic
)

printf 'Task %s: QP-Bilin-IC, loss=%s, seed=%s\n' "$task_id" "$loss_type" "$seed"
printf '%q ' "${command[@]}"
printf '\n'
if [[ "$dry_run" == true ]]; then
    exit 0
fi
exec "${command[@]}"
