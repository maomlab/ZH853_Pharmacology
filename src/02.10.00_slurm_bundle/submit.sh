#!/bin/bash
# Single entry point for every batch stage of the ZH853-MOR workflow.
#
#   ./submit.sh params         # step 5    CPU array: ligand parameters, one task per ligand
#   ./submit.sh build          # step 6    CPU array: system builds, one task per ligand x D2.50
#   ./submit.sh check          # step 0.5  pre-flight: can this env run OpenMM on CUDA?
#   ./submit.sh eq             # step 3    restrained ramp -> ${ZH_SYS}_eq.xml + eq_qc
#   ./submit.sh preprod        # step 3.5  unrestrained equilibration (DISCARDED) -> preprod_final.xml
#   ./submit.sh prod           # step 4    production job array (resumes from preprod_final.xml)
#   ./submit.sh all            # eq -> preprod -> prod, chained with --dependency=afterok
#   ./submit.sh <any> -n       # dry run: print the sbatch command, submit nothing
#
# WHERE TO RUN IT:
#   params, build, check  -- from this bundle directory (src/02.10.00_slurm_bundle)
#   eq, preprod, prod, all -- from a BUILD directory (intermediate/02.10.00_build/<name>/), which
#                             01_build_system.sh stages with everything those stages need.
#
# Site values come from cluster.env at the REPOSITORY ROOT and are passed to sbatch on the COMMAND
# LINE, which takes precedence over the `#SBATCH` directives inside the job scripts -- which is why
# those files carry no account/partition/time of their own. Per-build sampling (run lengths,
# replicas) comes from sampling.env in the build directory. See cluster_env.sh.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

die() { echo "ERROR: $*" >&2; exit 1; }

usage() {
  sed -n '2,21p' "${BASH_SOURCE[0]}" | sed -e 's/^#$//' -e 's/^# //'
  exit "${1:-0}"
}

# --- parse args -------------------------------------------------------------------------------
STAGE=""
DRY=0
for arg in "$@"; do
  case "$arg" in
    params|build|check|eq|preprod|prod|all) [ -z "$STAGE" ] || die "give only one stage, got '$STAGE' and '$arg'."; STAGE="$arg" ;;
    -n|--dry-run)      DRY=1 ;;
    -h|--help)         usage 0 ;;
    *)                 echo "ERROR: unknown argument '$arg'." >&2; usage 1 ;;
  esac
done
[ -n "$STAGE" ] || { echo "ERROR: no stage given." >&2; usage 1; }

# --- locate and load cluster.env --------------------------------------------------------------
# Shared with the job scripts so there is one definition of where cluster.env lives and one error
# message when it is absent. It exits non-zero rather than falling back to defaults.
# shellcheck source=cluster_env.sh
source "$HERE/cluster_env.sh" || exit 1
ENV_FILE="$ZH_CLUSTER_ENV"

# --- validate ---------------------------------------------------------------------------------
# Fail here, with a pointer to the file to edit, rather than letting sbatch reject the job with
# "Invalid account or account/partition combination specified" and no indication of where it came from.
[ -n "${ZH_ACCOUNT:-}" ]   || die "ZH_ACCOUNT is empty in $ENV_FILE (see: sacctmgr show assoc user=\$USER)."
[ -n "${ZH_GPU_PARTITION:-}" ] || die "ZH_GPU_PARTITION is empty in $ENV_FILE (see: sinfo -s)."

: "${ZH_GPU_GRES:=gpu:1}"       ; : "${ZH_GPU_CPUS:=8}"          ; : "${ZH_GPU_MEM:=32G}"
: "${ZH_GPU_EQ_TIME:=12:00:00}" ; : "${ZH_GPU_PROD_TIME:=72:00:00}" ; : "${ZH_GPU_CHECK_TIME:=00:05:00}"
: "${ZH_REPLICAS:=3}"       ; : "${ZH_PROD_NS:=500}"      ; : "${ZH_SYS:=system}"
: "${ZH_PREPROD_NS:=100}"   ; : "${ZH_GPU_PREPROD_TIME:=24:00:00}"
# Fixed by the six-stage schedule in 02_equilibrate.py: 1,125,000 steps x 2 fs = 2.25 ns.
EQ_NS=2.25

: "${ZH_CPU_TIME:=04:00:00}" ; : "${ZH_CPU_CPUS:=8}" ; : "${ZH_CPU_MEM:=32G}"
: "${ZH_CPU_PARTITION:=$ZH_GPU_PARTITION}"   # fall back to the GPU partition if no CPU one is set

# The CPU stages request no --gres: parameterization and building never touch a GPU, and holding
# one for hours of PACKMOL-Memgen would waste the allocation and queue behind GPU demand.
case "$STAGE" in
  params|build) _part="$ZH_CPU_PARTITION"; _gres="" ;;
  *)            _part="$ZH_GPU_PARTITION";     _gres="$ZH_GPU_GRES" ;;
esac
SBATCH_ARGS=(
  --account="$ZH_ACCOUNT"
  --partition="$_part"
  --export="ALL,ZH_CLUSTER_ENV=$ENV_FILE"
)
[ -n "$_gres" ] && SBATCH_ARGS+=(--gres="$_gres")
if [ -n "${ZH_QOS:-}" ];        then SBATCH_ARGS+=(--qos="$ZH_QOS"); fi
if [ -n "${ZH_CONSTRAINT:-}" ]; then SBATCH_ARGS+=(--constraint="$ZH_CONSTRAINT"); fi
if [ -n "${ZH_EXTRA_SBATCH:-}" ]; then
  # shellcheck disable=SC2206 -- deliberate word-splitting: ZH_EXTRA_SBATCH holds whole flags
  SBATCH_ARGS+=($ZH_EXTRA_SBATCH)
fi

# --- wall-clock estimates ----------------------------------------------------------------------
# No separate calibration run is needed: StateDataReporter writes a "Speed (ns/day)" column, so
# every completed stage in this directory measures this cluster's rate for this system size.
# Pure awk deliberately -- submit.sh runs on the login node and activates no conda env.
#
# Two rates matter and they are NOT interchangeable: equilibration runs at 2 fs with positional
# restraints, production and pre-production at 4 fs with HMR, so production covers roughly twice
# the simulated time per wall-clock hour.

_speed() {  # <log> -> mean ns/day over the last 20 samples, or empty
  [ -f "$1" ] || return 0
  awk -F',' '
    NR == 1 { for (i = 1; i <= NF; i++) { h = $i; gsub(/[#"]/, "", h)
                                          if (h ~ /^Speed/) c = i }
              next }
    c && NF >= c && $c + 0 > 0 { v[++n] = $c + 0 }
    END { if (n < 3) exit
          s = 0; k = 0
          for (i = (n > 20 ? n - 19 : 1); i <= n; i++) { s += v[i]; k++ }
          printf "%.1f", s / k }' "$1" 2>/dev/null || true
}

_newest() {  # most recently modified of the given globs (preprod2 beats preprod, etc.)
  # shellcheck disable=SC2012 -- ls -t is fine here; these are our own generated log names
  ls -t $@ 2>/dev/null | head -1 || true
}

_hours() {  # "48:00:00" or "2-12:00:00" -> hours
  awk -v t="$1" 'BEGIN { d = 0
    if (index(t, "-") > 0) { split(t, a, "-"); d = a[1]; t = a[2] }
    n = split(t, b, ":")
    printf "%.3f", d * 24 + (n >= 1 ? b[1] : 0) + (n >= 2 ? b[2] : 0) / 60 + (n >= 3 ? b[3] : 0) / 3600 }'
}

_dur() {  # hours -> "3.4 h" / "1 d 5.2 h"
  awk -v h="$1" 'BEGIN { if (h < 24) printf "%.1f h", h
                         else printf "%d d %.1f h", int(h / 24), h - 24 * int(h / 24) }'
}

# Rates measured in THIS build directory. eq: 2 fs restrained. prod/preprod: 4 fs HMR.
EQ_LOG="${ZH_SYS}_eq.log"
PROD_LOG="$(_newest 'preprod*.log' 'prod_r*.log')"
EQ_RATE="$(_speed "$EQ_LOG")"
PROD_RATE="$(_speed "$PROD_LOG")"
# One stage's rate is NOT convertible into the other's. It is tempting to scale by the timestep
# (2 fs -> 4 fs, so 2x), and this script used to: measured on the H200 nodes the real ratio is
# 6.2x (94.4 vs 587 ns/day). The timestep is only one factor -- equilibration also carries the
# restraint force, and it writes an energy-bearing state report every 2,500 steps against
# production's 25,000, each of which forces a GPU sync and a full energy evaluation. A derived
# number would have been wrong by 3x while looking authoritative, so an unmeasured stage now
# reports "not measured" and the other rate is shown only as context.
if [ -n "$EQ_RATE" ]; then EQ_TXT="${EQ_RATE} ns/day @2fs (${EQ_LOG})"
else EQ_TXT="2 fs not measured yet"; fi
if [ -n "$PROD_RATE" ]; then PROD_TXT="${PROD_RATE} ns/day @4fs (${PROD_LOG})"
else PROD_TXT="4 fs not measured yet"; fi
if [ -n "$EQ_RATE" ] || [ -n "$PROD_RATE" ]; then
  echo "measured rate: ${EQ_TXT}; ${PROD_TXT}"
else
  echo "measured rate: none yet -- the first completed stage writes a Speed column and calibrates this"
fi

estimate() {  # <ns> <rate> <walltime> <label>
  local ns="$1" rate="$2" wall="$3" label="$4"
  [ -n "$rate" ] || { echo "  estimate: unknown -- no completed run at this timestep yet"; return 0; }
  local hrs; hrs="$(awk -v n="$ns" -v r="$rate" 'BEGIN{printf "%.3f", 24*n/r}')"
  local wh;  wh="$(_hours "$wall")"
  local pct; pct="$(awk -v h="$hrs" -v w="$wh" 'BEGIN{printf "%.0f", (w>0? 100*h/w : 0)}')"
  echo "  estimate: ~$(_dur "$hrs") for ${ns} ns ${label} (wall-time ${wall}, ~${pct}% of it)"
  # A job that dies at the limit loses everything since the last checkpoint, so say so up front.
  if [ "$(awk -v h="$hrs" -v w="$wh" 'BEGIN{print (h > 0.85*w) ? 1 : 0}')" = "1" ]; then
    echo "  WARNING: that is over 85% of the requested wall-time -- raise it in cluster.env," \
         "or expect to restart from the checkpoint."
  fi
}

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

N_SYSTEMS=$(python "$HERE/ligands.py" --list | wc -w | tr -d ' ')
N_LIGANDS=$(python "$HERE/ligands.py" --list-ligands | wc -w | tr -d ' ')

case "$STAGE" in
  params)
    need ligands.py
    n=$N_LIGANDS                  # registry entries that have a ligand
    jid=$(submit submit_parameterize.sbatch --time="$ZH_CPU_TIME" \
                 --cpus-per-task="$ZH_CPU_CPUS" --mem="$ZH_CPU_MEM" --array="1-${n}" \
                 --job-name=zh853_params --output=params_%A_%a.out)
    echo "submitted parameterization: $jid   ($n ligands, one per array task)"
    echo "  -> intermediate/02.08.00_ligand_params/<LIGAND>/"
    echo "then: ./submit.sh build"
    ;;
  build)
    need ligands.py
    n=$((N_SYSTEMS * 2))          # every system x {ASP, ASH}
    jid=$(submit submit_build.sbatch --time="$ZH_CPU_TIME" \
                 --cpus-per-task="$ZH_CPU_CPUS" --mem="$ZH_CPU_MEM" --array="1-${n}" \
                 --job-name=zh853_build --output=build_%A_%a.out)
    echo "submitted builds: $jid   ($N_SYSTEMS systems x 2 D2.50 states = $n array tasks)"
    echo "  -> intermediate/02.10.00_build/<LIGAND>_<D250>_<timestamp>/"
    echo "then, from each build directory: ./submit.sh all"
    ;;
  check)
    # The pre-flight is tiny; do not hold a full training-sized allocation for it.
    jid=$(submit check_gpu_env.sh --time="$ZH_GPU_CHECK_TIME" --cpus-per-task=2 --mem=8G \
                                  --job-name=zh853_gpucheck --output=gpucheck_%j.out)
    echo "submitted pre-flight: $jid   (watch: tail -f gpucheck_${jid}.out)"
    ;;
  eq)
    need "${ZH_SYS}.prmtop"; need "${ZH_SYS}.rst7"
    jid=$(submit submit_equilibrate.sbatch --time="$ZH_GPU_EQ_TIME" --cpus-per-task="$ZH_GPU_CPUS" --mem="$ZH_GPU_MEM" \
                 --job-name=zh853_eq --output=eq_%j.out)
    echo "submitted equilibration: $jid   -> ${ZH_SYS}_eq.xml + eq_qc.{json,png}"
    estimate "$EQ_NS" "$EQ_RATE" "$ZH_GPU_EQ_TIME" "over 6 restrained stages at 2 fs"
    echo "then: ./submit.sh preprod"
    ;;
  preprod)
    need "${ZH_SYS}.prmtop"
    [ -f "${ZH_SYS}_eq.xml" ] || echo "WARNING: ${ZH_SYS}_eq.xml not present yet -- step 3 must finish first."
    jid=$(submit submit_preproduction.sbatch --time="$ZH_GPU_PREPROD_TIME" \
                 --cpus-per-task="$ZH_GPU_CPUS" --mem="$ZH_GPU_MEM" \
                 --job-name=zh853_preprod --output=preprod_%j.out)
    echo "submitted pre-production: $jid   ($ZH_PREPROD_NS ns, unrestrained, discarded)"
    estimate "$ZH_PREPROD_NS" "$PROD_RATE" "$ZH_GPU_PREPROD_TIME" "unrestrained at 4 fs, per leg"
    echo "then: ./submit.sh prod"
    ;;
  prod)
    need "${ZH_SYS}.prmtop"
    [ -f "${ZH_SYS}_eq.xml" ] || echo "WARNING: ${ZH_SYS}_eq.xml not present yet -- step 3 must finish first."
    jid=$(submit submit_production.sbatch --time="$ZH_GPU_PROD_TIME" --cpus-per-task="$ZH_GPU_CPUS" --mem="$ZH_GPU_MEM" --array="1-${ZH_REPLICAS}" \
                 --job-name=zh853_prod --output=prod_%A_%a.out)
    echo "submitted production: $jid   ($ZH_REPLICAS replicas x $ZH_PROD_NS ns)"
    estimate "$ZH_PROD_NS" "$PROD_RATE" "$ZH_GPU_PROD_TIME" "at 4 fs, PER REPLICA"
    echo "  (array tasks run concurrently if the queue has $ZH_REPLICAS GPUs free, serially otherwise)"
    ;;
  all)
    need "${ZH_SYS}.prmtop"; need "${ZH_SYS}.rst7"
    # afterok, not afterany, at both links: each stage loads the state the previous one wrote, and
    # each runs check_equilibration.py, which exits non-zero on a FAIL. So a broken system stops
    # the chain here instead of consuming ZH_REPLICAS x ZH_PROD_NS ns of GPU time.
    eq=$(submit submit_equilibrate.sbatch --time="$ZH_GPU_EQ_TIME" --cpus-per-task="$ZH_GPU_CPUS" --mem="$ZH_GPU_MEM" \
                --job-name=zh853_eq --output=eq_%j.out)
    echo "submitted equilibration:  $eq   -> ${ZH_SYS}_eq.xml + eq_qc.{json,png}"
    estimate "$EQ_NS" "$EQ_RATE" "$ZH_GPU_EQ_TIME" "over 6 restrained stages at 2 fs"
    dep_eq=(--dependency="afterok:$eq")
    if [ "$DRY" -eq 1 ]; then dep_eq=(--dependency="afterok:<eq_jobid>"); fi
    pre=$(submit submit_preproduction.sbatch --time="$ZH_GPU_PREPROD_TIME" \
                 --cpus-per-task="$ZH_GPU_CPUS" --mem="$ZH_GPU_MEM" \
                 --job-name=zh853_preprod --output=preprod_%j.out "${dep_eq[@]}")
    echo "submitted pre-production: $pre   ($ZH_PREPROD_NS ns unrestrained, discarded)"
    estimate "$ZH_PREPROD_NS" "$PROD_RATE" "$ZH_GPU_PREPROD_TIME" "unrestrained at 4 fs, per leg"
    dep_pre=(--dependency="afterok:$pre")
    if [ "$DRY" -eq 1 ]; then dep_pre=(--dependency="afterok:<preprod_jobid>"); fi
    prod=$(submit submit_production.sbatch --time="$ZH_GPU_PROD_TIME" --cpus-per-task="$ZH_GPU_CPUS" --mem="$ZH_GPU_MEM" --array="1-${ZH_REPLICAS}" \
                  --job-name=zh853_prod --output=prod_%A_%a.out "${dep_pre[@]}")
    echo "submitted production:     $prod   ($ZH_REPLICAS replicas x $ZH_PROD_NS ns)"
    estimate "$ZH_PROD_NS" "$PROD_RATE" "$ZH_GPU_PROD_TIME" "at 4 fs, PER REPLICA"
    echo
    echo "total to last replica: ~$(awk -v e="$EQ_NS" -v er="${EQ_RATE:-0}" -v p="$ZH_PREPROD_NS" \
        -v q="$ZH_PROD_NS" -v pr="${PROD_RATE:-0}" \
        'BEGIN{ if (er<=0 || pr<=0) { print "unknown"; exit } h=24*e/er + 24*p/pr + 24*q/pr;
                if (h<24) printf "%.1f h", h; else printf "%d d %.1f h", int(h/24), h-24*int(h/24) }') " \
         "(sequential; replicas overlap if GPUs are free)"
    echo "watch: squeue -u \$USER"
    ;;
esac
