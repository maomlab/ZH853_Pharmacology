#!/bin/bash
# Single entry point for the GPU stages (0.5, 3, 4) of the ZH853-MOR bundle.
#
#   ./submit.sh check          # step 0.5  pre-flight: can this env run OpenMM on CUDA?
#   ./submit.sh eq             # step 3    restrained equilibration -> ${ZH_SYS}_eq.xml
#   ./submit.sh prod           # step 4    production job array (needs eq to have finished)
#   ./submit.sh all            # eq, then prod chained with --dependency=afterok
#   ./submit.sh <any> -n       # dry run: print the sbatch command, submit nothing
#
# Run it from the BUILD DIRECTORY (intermediate/02.10.00_build/<D250>_<timestamp>/), which
# 01_build_system.sh stages with everything this needs.
#
# All cluster-specific values come from cluster.env (see cluster.env.example) and are passed to
# sbatch on the COMMAND LINE, which takes precedence over the `#SBATCH` directives inside the job
# scripts. That is why the .sbatch files carry no account/partition/time of their own: there is one
# source of truth per cluster, not one per job script per build.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

die() { echo "ERROR: $*" >&2; exit 1; }

usage() {
  sed -n '2,16p' "${BASH_SOURCE[0]}" | sed -e 's/^#$//' -e 's/^# //'
  exit "${1:-0}"
}

# --- parse args -------------------------------------------------------------------------------
STAGE=""
DRY=0
for arg in "$@"; do
  case "$arg" in
    check|eq|prod|all) [ -z "$STAGE" ] || die "give only one stage, got '$STAGE' and '$arg'."; STAGE="$arg" ;;
    -n|--dry-run)      DRY=1 ;;
    -h|--help)         usage 0 ;;
    *)                 echo "ERROR: unknown argument '$arg'." >&2; usage 1 ;;
  esac
done
[ -n "$STAGE" ] || { echo "ERROR: no stage given." >&2; usage 1; }

# --- locate and load cluster.env --------------------------------------------------------------
# Explicit override, then the build dir we were invoked from, then the bundle's own copy.
ENV_FILE=""
for cand in "${ZH_CLUSTER_ENV:-}" "$PWD/cluster.env" "$HERE/cluster.env"; do
  if [ -n "$cand" ] && [ -f "$cand" ]; then ENV_FILE="$(cd "$(dirname "$cand")" && pwd)/$(basename "$cand")"; break; fi
done
if [ -z "$ENV_FILE" ]; then
  echo "ERROR: no cluster.env found (looked in \$ZH_CLUSTER_ENV, $PWD, $HERE)." >&2
  echo "  Create one from the template -- it is gitignored, so it is per-machine and per-user:" >&2
  echo "      cp $HERE/cluster.env.example $HERE/cluster.env && \$EDITOR $HERE/cluster.env" >&2
  exit 1
fi
# shellcheck disable=SC1090
source "$ENV_FILE"
echo "cluster settings: $ENV_FILE"

# --- validate ---------------------------------------------------------------------------------
# Fail here, with a pointer to the file to edit, rather than letting sbatch reject the job with
# "Invalid account or account/partition combination specified" and no indication of where it came from.
[ -n "${ZH_ACCOUNT:-}" ]   || die "ZH_ACCOUNT is empty in $ENV_FILE (see: sacctmgr show assoc user=\$USER)."
[ -n "${ZH_PARTITION:-}" ] || die "ZH_PARTITION is empty in $ENV_FILE (see: sinfo -s)."

: "${ZH_GRES:=gpu:1}"       ; : "${ZH_CPUS:=8}"          ; : "${ZH_MEM:=32G}"
: "${ZH_EQ_TIME:=12:00:00}" ; : "${ZH_PROD_TIME:=48:00:00}" ; : "${ZH_CHECK_TIME:=00:05:00}"
: "${ZH_REPLICAS:=3}"       ; : "${ZH_PROD_NS:=500}"      ; : "${ZH_SYS:=system}"

SBATCH_ARGS=(
  --account="$ZH_ACCOUNT"
  --partition="$ZH_PARTITION"
  --gres="$ZH_GRES"
  --export="ALL,ZH_CLUSTER_ENV=$ENV_FILE"
)
if [ -n "${ZH_QOS:-}" ];        then SBATCH_ARGS+=(--qos="$ZH_QOS"); fi
if [ -n "${ZH_CONSTRAINT:-}" ]; then SBATCH_ARGS+=(--constraint="$ZH_CONSTRAINT"); fi
if [ -n "${ZH_EXTRA_SBATCH:-}" ]; then
  # shellcheck disable=SC2206 -- deliberate word-splitting: ZH_EXTRA_SBATCH holds whole flags
  SBATCH_ARGS+=($ZH_EXTRA_SBATCH)
fi

# --- input checks -----------------------------------------------------------------------------
need() { [ -f "$1" ] || die "$1 not found in $PWD. Steps 3-5 run from the build directory; see README.md."; }

submit() {  # submit <script> <extra sbatch args...>; echoes the job id
  local script="$1"; shift
  need "$script"
  if [ "$DRY" -eq 1 ]; then
    echo "DRY RUN: sbatch ${SBATCH_ARGS[*]} $* $script" >&2
    echo "DRYRUN"
    return 0
  fi
  sbatch --parsable "${SBATCH_ARGS[@]}" "$@" "$script"
}

case "$STAGE" in
  check)
    # The pre-flight is tiny; do not hold a full training-sized allocation for it.
    jid=$(submit check_gpu_env.sh --time="$ZH_CHECK_TIME" --cpus-per-task=2 --mem=8G \
                                  --job-name=zh853_gpucheck --output=gpucheck_%j.out)
    echo "submitted pre-flight: $jid   (watch: tail -f gpucheck_${jid}.out)"
    ;;
  eq)
    need "${ZH_SYS}.prmtop"; need "${ZH_SYS}.rst7"
    jid=$(submit submit_equilibrate.sbatch --time="$ZH_EQ_TIME" --cpus-per-task="$ZH_CPUS" --mem="$ZH_MEM" \
                 --job-name=zh853_eq --output=eq_%j.out)
    echo "submitted equilibration: $jid   -> ${ZH_SYS}_eq.xml"
    echo "then: ./submit.sh prod   (or resubmit with: sbatch --dependency=afterok:$jid ...)"
    ;;
  prod)
    need "${ZH_SYS}.prmtop"
    [ -f "${ZH_SYS}_eq.xml" ] || echo "WARNING: ${ZH_SYS}_eq.xml not present yet -- step 3 must finish first."
    jid=$(submit submit_production.sbatch --time="$ZH_PROD_TIME" --cpus-per-task="$ZH_CPUS" --mem="$ZH_MEM" --array="1-${ZH_REPLICAS}" \
                 --job-name=zh853_prod --output=prod_%A_%a.out)
    echo "submitted production: $jid   ($ZH_REPLICAS replicas x $ZH_PROD_NS ns)"
    ;;
  all)
    need "${ZH_SYS}.prmtop"; need "${ZH_SYS}.rst7"
    eq=$(submit submit_equilibrate.sbatch --time="$ZH_EQ_TIME" --cpus-per-task="$ZH_CPUS" --mem="$ZH_MEM" \
                --job-name=zh853_eq --output=eq_%j.out)
    echo "submitted equilibration: $eq   -> ${ZH_SYS}_eq.xml"
    # afterok, not afterany: production loads the state equilibration writes, so a failed or
    # cancelled eq must not launch replicas that would die on a missing/half-written .xml.
    dep=(--dependency="afterok:$eq")
    if [ "$DRY" -eq 1 ]; then dep=(--dependency="afterok:<eq_jobid>"); fi
    prod=$(submit submit_production.sbatch --time="$ZH_PROD_TIME" --cpus-per-task="$ZH_CPUS" --mem="$ZH_MEM" --array="1-${ZH_REPLICAS}" \
                  --job-name=zh853_prod --output=prod_%A_%a.out "${dep[@]}")
    echo "submitted production:    $prod   (held until $eq succeeds; $ZH_REPLICAS replicas x $ZH_PROD_NS ns)"
    echo "watch: squeue -u \$USER"
    ;;
esac
