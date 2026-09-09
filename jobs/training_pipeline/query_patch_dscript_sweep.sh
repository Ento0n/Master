#!/bin/bash
#
#SBATCH --job-name=query-patch-dscript-sweep
#SBATCH --array=0-26%4
#SBATCH --output=logs/query-patch-dscript-sweep-%A_%a.out
#SBATCH --error=logs/query-patch-dscript-sweep-%A_%a.err
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

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <sweep-id>" >&2
    exit 1
fi

sweep_id="$1"

# One agent gets one configuration, runs one Lightning trial, then exits.
cd /nfs/home/students/a.spannagl/master_repository
srun wandb agent \
    --forward-signals \
    --count 1 \
    "$sweep_id"
