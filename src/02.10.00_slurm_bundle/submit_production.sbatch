#!/bin/bash
#SBATCH --job-name=zh853_prod
#SBATCH --output=prod_%A_%a.out
#SBATCH --array=1-3                 # 3 production replicas (SPECIFICATION QC: replicate spread)
#SBATCH --gres=gpu:1               # TODO(OQ-3): GPU type, e.g. gpu:a100:1
#SBATCH --partition=gpu            # TODO(OQ-3): partition/queue name
#SBATCH --account=CHANGEME         # TODO(OQ-3): allocation/account
#SBATCH --time=48:00:00            # TODO(OQ-3): wall-time limit
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail
module purge
# TODO(OQ-3): module load cuda/<version>   # match environment-cluster.yml build
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate zh853mor-sim

SYS=${SYS:-system}                 # basename of prmtop/rst7 and eq state
REP=${SLURM_ARRAY_TASK_ID:-1}

python 03_production.py \
  --prmtop "${SYS}.prmtop" \
  --state  "${SYS}_eq.xml" \
  --out    "prod_r${REP}" \
  --ns 500 --seed "${REP}"
