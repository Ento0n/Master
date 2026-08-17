#!/bin/bash
#
#SBATCH --job-name=14_generate_embeddings                     # Job name
#SBATCH --output=logs/14_generate_embeddings.%j.out           # Standard output (%j expands to jobId)
#SBATCH --error=logs/14_generate_embeddings.%j.err            # Standard error
#SBATCH --partition=shared-gpu                      # not standard otherwise no permissions
#SBATCH --gres=gpu:1
#SBATCH --qos=limitgpus
#SBATCH --ntasks=1                                  # Run a single task
#SBATCH --cpus-per-task=1                           # Number of CPU cores per task
#SBATCH --mem=64G                                   # Total memory
#SBATCH --time=4-00:00:00                           # Time limit hh:mm:ss
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=antspa@gmx.de                   # Where to send mail

set -euo pipefail

python -c 'import os, sys, torch; print(f"python={sys.executable}"); print("CUDA_VISIBLE_DEVICES=" + os.environ.get("CUDA_VISIBLE_DEVICES", "<unset>")); print(f"torch={torch.__version__} torch_cuda={torch.version.cuda} cuda_available={torch.cuda.is_available()} cuda_devices={torch.cuda.device_count()}")'
python /nfs/home/students/a.spannagl/master_repository/scripts/data_pipeline/14_generate_embeddings.py --device cuda
