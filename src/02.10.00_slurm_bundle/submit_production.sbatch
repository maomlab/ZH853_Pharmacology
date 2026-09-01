#!/bin/bash
#SBATCH --job-name=zh853_prod
#SBATCH --output=prod_%A_%a.out
#
# Step 4: restraint-free production, one replica per array task. Submit with the wrapper, from the
# BUILD DIRECTORY, after step 3 has finished:
#
#     ./submit.sh prod          # or: ./submit.sh all-simulations   (chains eq -> prod)
#
# Account, partition, GPU, wall-time and the array size (ZH_REPLICAS, SPECIFICATION QC wants >=3
# for a replicate spread) are NOT set here on purpose -- they live in cluster.env and submit.sh
# passes them as sbatch command-line flags, which override any #SBATCH directive.
#
# Reads preprod_final.xml if the unrestrained pre-production leg (step 3.5) has run, otherwise
# ${ZH_SYS}_eq.xml straight from the ramp. Set ZH_PROD_STATE to force one.

set -euo pipefail

# cluster.env is REQUIRED: running on built-in defaults would use a different conda env and a
# different system basename and still exit 0. cluster_env.sh fails loudly instead.
# shellcheck source=cluster_env.sh
source "${SLURM_SUBMIT_DIR:-$PWD}/cluster_env.sh" || exit 1

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

# Prefer the state that has actually equilibrated. Starting from the ramp alone means starting
# from a bilayer whose area per lipid has not converged and a backbone still weakly restrained,
# so say loudly which one is being used rather than letting it pass unnoticed.
STATE="${ZH_PROD_STATE:-}"
if [ -z "$STATE" ]; then
  # Newest pre-production leg (preprod, preprod2, ...), since the stage can be re-run to extend.
  NEWEST=""
  L=1
  while :; do
    N=$([ "$L" -eq 1 ] && echo preprod || echo "preprod${L}")
    [ -f "${N}_final.xml" ] || break
    NEWEST="${N}_final.xml"
    L=$((L + 1))
  done
  if [ -n "$NEWEST" ]; then
    STATE="$NEWEST"
  else
    STATE="${SYS}_eq.xml"
    echo "WARNING: preprod_final.xml not found -- starting production straight from the restrained"
    echo "         ramp. The bilayer is probably not equilibrated; see ./submit.sh preprod."
  fi
fi
[ -f "$STATE" ] || { echo "ERROR: state file $STATE not found in $PWD."; exit 1; }

echo "host: $(hostname)  job: ${SLURM_JOB_ID:-none}  replica: ${REP}  ${NS} ns  state: ${STATE}"
python -c "import openmm; print('openmm', openmm.version.version)"

python 03_production.py \
  --prmtop "${SYS}.prmtop" \
  --state  "${STATE}" \
  --out    "prod_r${REP}" \
  --ns "${NS}" --seed "${REP}"
