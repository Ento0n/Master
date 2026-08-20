#!/bin/bash
#
#SBATCH --job-name=test_query_patch                     # Job name
#SBATCH --output=logs/test_query_patch.%j.out           # Standard output (%j expands to jobId)
#SBATCH --error=logs/test_query_patch.%j.err            # Standard error
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

# An optional first argument overrides the automatically generated result
# subdirectory, for example: sbatch esm2_query_patch.sh my_debug_run
output_subdir_args=()
if [[ $# -ge 1 ]]; then
    output_subdir_args=(--output-subdir "$1")
fi

cd /nfs/home/students/a.spannagl/master_repository
python -m scripts.training_pipeline.execute_pipeline \
  --model query_patch \
  --int-mod-type max \
  --loss-type balanced_bce \
  --max-epochs 80 \
  "${output_subdir_args[@]}" \
