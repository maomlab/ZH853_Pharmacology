#!/bin/bash
#SBATCH --job-name=zh853_eq
#SBATCH --output=eq_%j.out
#
# Step 3: restrained equilibration. Submit with the wrapper, from the BUILD DIRECTORY
# (intermediate/02.10.00_build/<D250>_<timestamp>/, staged by 01_build_system.sh):
#
#     ./submit.sh eq
#
# Account, partition, GPU, wall-time, CPUs and memory are NOT set here on purpose -- they live in
# cluster.env and submit.sh passes them as sbatch command-line flags, which override any #SBATCH
# directive. One file to edit per cluster, instead of a CHANGEME in every job script.
# A bare `sbatch submit_equilibrate.sbatch` will therefore inherit your site defaults and probably
# land on a CPU partition with no GPU; use the wrapper.
#
# Writes ${ZH_SYS}_eq.xml, which submit_production.sbatch reads back.

set -euo pipefail

# Exported by submit.sh so the job picks up the same settings the submission used.
if [ -n "${ZH_CLUSTER_ENV:-}" ] && [ -f "${ZH_CLUSTER_ENV}" ]; then
  # shellcheck disable=SC1090
  source "${ZH_CLUSTER_ENV}"
fi

module purge
# OpenMM (conda-forge) bundles its own CUDA runtime and needs only the node's driver, so ZH_MODULES
# is normally empty -- see the CUDA section of README.md before setting it.
if [ -n "${ZH_MODULES:-}" ]; then
  # shellcheck disable=SC2086 -- ZH_MODULES is a space-separated module list
  module load ${ZH_MODULES}
fi

source "${ZH_CONDA_SH:-$(conda info --base)/etc/profile.d/conda.sh}"
conda activate "${ZH_SIM_ENV:-zh853mor-sim}"

SYS="${ZH_SYS:-system}"
echo "host: $(hostname)  job: ${SLURM_JOB_ID:-none}  system: ${SYS}"
python -c "import openmm; print('openmm', openmm.version.version)"

python 02_equilibrate.py --prmtop "${SYS}.prmtop" --inpcrd "${SYS}.rst7" --out "${SYS}_eq"
