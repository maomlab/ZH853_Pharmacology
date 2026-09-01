#!/bin/bash
#SBATCH --job-name=zh853_params
#SBATCH --output=params_%A_%a.out
#
# CPU stage: GAFF2/AM1-BCC parameters, ONE LIGAND PER ARRAY TASK so the four run concurrently
# instead of ~4 x sqm serially on a login node.
#
#     ./submit.sh params
#
# Resources and the array size come from cluster.env at the repository root; submit.sh passes them
# as sbatch command-line flags. Each task parameterizes the ligand at its array index, taken from
# the ligands.py registry so the mapping cannot drift from the panel.

set -euo pipefail

# shellcheck source=cluster_env.sh
source "${SLURM_SUBMIT_DIR:-$PWD}/cluster_env.sh" || exit 1

BUNDLE="${SLURM_SUBMIT_DIR:-$PWD}"
REPO="$(cd "$BUNDLE/../.." && pwd)"

module purge
if [ -n "${ZH_MODULES:-}" ]; then
  # shellcheck disable=SC2086 -- ZH_MODULES is a space-separated module list
  module load ${ZH_MODULES}
fi
source "${ZH_CONDA_SH:-$(conda info --base)/etc/profile.d/conda.sh}"
conda activate "${ZH_PREP_ENV:-zh853mor-prep}"   # antechamber/parmchk2 live here

# Registry order, minus the apo pseudo-system which has no ligand to parameterize. The registry
# does the filtering: `sed 's/\bapo\b//'` looks right but \b is a GNU extension that BSD sed
# ignores, which left apo as array task 1 and pushed the last ligand off the end.
read -r -a LIGANDS <<< "$(python "$BUNDLE/ligands.py" --list-ligands)"
IDX="${SLURM_ARRAY_TASK_ID:-1}"
LIGAND="${LIGANDS[$((IDX - 1))]:-}"
[ -n "$LIGAND" ] || { echo "ERROR: no ligand at array index $IDX (have ${#LIGANDS[@]})."; exit 1; }

echo "host: $(hostname)  job: ${SLURM_JOB_ID:-none}  task ${IDX}/${#LIGANDS[@]}  ligand: ${LIGAND}"
bash "$REPO/src/02.08.00_ligand_parameterize.sh" "$LIGAND"
