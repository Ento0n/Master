#!/bin/bash
#
#SBATCH --job-name=dscript-sweep                                    # Job name
#SBATCH --array=0-27%4
#SBATCH --output=logs/sweep-%A_%a.out                               # Standard output (%j expands to jobId)
#SBATCH --error=logs/sweep-%A_%a.err                                # Standard error
#SBATCH --partition=shared-gpu                                      # not standard otherwise no permissions
#SBATCH --gres=gpu:1
#SBATCH --qos=limitgpus
#SBATCH --exclude=gpu01.exbio.wzw.tum.de,compms-gpu-1.exbio.wzw.tum.de
#SBATCH --ntasks=1                                                  # Run a single task
#SBATCH --cpus-per-task=4                                           # Number of CPU cores per task
#SBATCH --mem=64G                                                   # Total memory
#SBATCH --time=4-00:00:00                                           # Time limit hh:mm:ss
#SBATCH --mail-type=END,FAIL                                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=antspa@gmx.de                                   # Where to send mail

set -euo pipefail

if [[ $# -lt 1 ]]; then 
    echo "Usage: $0 <string>" >&2
    exit 1
fi

sweep_id="$1"

# One agent gets one configuration, runs one Lightning trial, then exits.
cd /nfs/home/students/a.spannagl/master_repository
srun wandb agent \
    --forward-signals \
    --count 1 \
    "$sweep_id"