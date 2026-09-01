#!/bin/bash
#SBATCH --job-name=zh853_prod
#SBATCH --output=prod_%A_%a.out
#
# Step 4: restraint-free production, one replica per array task. Submit with the wrapper, from the
# BUILD DIRECTORY, after step 3 has finished:
#
#     ./submit.sh prod          # or: ./submit.sh all   (chains eq -> prod with afterok)
#
# Account, partition, GPU, wall-time and the array size (ZH_REPLICAS, SPECIFICATION QC wants >=3
# for a replicate spread) are NOT set here on purpose -- they live in cluster.env and submit.sh
# passes them as sbatch command-line flags, which override any #SBATCH directive.
#
# Reads ${ZH_SYS}_eq.xml, written by submit_equilibrate.sbatch.

set -euo pipefail

# Exported by submit.sh so the job picks up the same settings the submission used.
if [ -n "${ZH_CLUSTER_ENV:-}" ] && [ -f "${ZH_CLUSTER_ENV}" ]; then
  # shellcheck disable=SC1090
  source "${ZH_CLUSTER_ENV}"
fi

module purge
if [ -n "${ZH_MODULES:-}" ]; then
  # shellcheck disable=SC2086 -- ZH_MODULES is a space-separated module list
  module load ${ZH_MODULES}
fi

source "${ZH_CONDA_SH:-$(conda info --base)/etc/profile.d/conda.sh}"
conda activate "${ZH_SIM_ENV:-zh853mor-sim}"

SYS="${ZH_SYS:-system}"          # basename of prmtop/rst7 and the equilibrated state
REP="${SLURM_ARRAY_TASK_ID:-1}"  # replica index, also the RNG seed
NS="${ZH_PROD_NS:-500}"

echo "host: $(hostname)  job: ${SLURM_JOB_ID:-none}  replica: ${REP}  ${NS} ns"
python -c "import openmm; print('openmm', openmm.version.version)"

python 03_production.py \
  --prmtop "${SYS}.prmtop" \
  --state  "${SYS}_eq.xml" \
  --out    "prod_r${REP}" \
  --ns "${NS}" --seed "${REP}"
