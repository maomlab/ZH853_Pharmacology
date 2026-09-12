#!/bin/bash
#SBATCH --job-name=zh853_reduce
#SBATCH --output=reduce_%A_%a.out
#
# CPU stage: reduce ONE production replica per array task. Reading a 500 ns / 5 GB DCD is I/O and
# single-core bound, so the replicas of the whole panel run concurrently instead of serially, and
# no task needs a GPU.
#
#     ./submit_reduce.sh
#
# The array index selects the (build, replica) pair from the SAME listing the submit script sized
# the array with -- 01_reduce_trajectory.py --list -- so the mapping cannot drift between the two.

set -euo pipefail

# The site settings live at the repository root; this stage has no sampling.env of its own.
# shellcheck source=../02.10.00_slurm_bundle/cluster_env.sh
source "${SLURM_SUBMIT_DIR:-$PWD}/../02.10.00_slurm_bundle/cluster_env.sh" || exit 1

STAGE="${SLURM_SUBMIT_DIR:-$PWD}"

module purge
if [ -n "${ZH_MODULES:-}" ]; then
  # shellcheck disable=SC2086 -- ZH_MODULES is a space-separated module list
  module load ${ZH_MODULES}
fi
source "${ZH_CONDA_SH:-$(conda info --base)/etc/profile.d/conda.sh}"
conda activate "${ZH_PREP_ENV:-zh853mor-prep}"   # MDAnalysis lives here; no OpenMM needed

IDX="${SLURM_ARRAY_TASK_ID:-1}"
echo "host: $(hostname)  job: ${SLURM_JOB_ID:-none}  task ${IDX}"
python "$STAGE/01_reduce_trajectory.py" --index "$IDX" "$@"
