#!/bin/bash
#SBATCH --job-name=zh853_eq
#SBATCH --output=eq_%j.out
#SBATCH --gres=gpu:1               # TODO(OQ-3): GPU type
#SBATCH --partition=gpu            # TODO(OQ-3): partition/queue
#SBATCH --account=CHANGEME         # TODO(OQ-3): allocation/account
#SBATCH --time=12:00:00            # TODO(OQ-3): wall-time
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail
module purge
# TODO(OQ-3): module load cuda/<version>
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate zh853mor-sim

SYS=${SYS:-system}
python 02_equilibrate.py --prmtop "${SYS}.prmtop" --inpcrd "${SYS}.rst7" --out "${SYS}_eq"
