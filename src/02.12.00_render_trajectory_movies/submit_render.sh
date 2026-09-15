#!/bin/bash
# Submit the movie rendering as a job array, one exported replica per task.
#
#     ./submit_render.sh                    # both passes, every exported replica
#     ./submit_render.sh -n                 # dry run: print the sbatch command, submit nothing
#     ./submit_render.sh --only dashboard   # matplotlib diagnostics only (no node needed)
#     ./submit_render.sh --only molstar     # cartoon render only
#     ./submit_render.sh -- --fps 30        # anything after -- is passed to 02_render_dashboard.py
#
# Rendering is single-core and serial WITHIN a replica: matplotlib draws every frame in turn and
# pipes it to ffmpeg, and MolStar draws its frames one at a time in software WebGL. Threading one
# replica is not available, so the panel is made fast by rendering the replicas concurrently --
# which is what this array does. Fifteen replicas that took an hour end to end take about as long
# as the slowest one.
#
# Resources come from cluster.env at the repository ROOT (the CPU settings, ZH_CPU_*), the same
# file every other stage uses. The array size is taken from the stage's own listing of exported
# movies, so adding a replica needs no edit here.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
die() { echo "ERROR: $*" >&2; exit 1; }

DRY=0
ONLY=both
while [ $# -gt 0 ]; do
  case "$1" in
    -n) DRY=1; shift ;;
    --only) ONLY="${2:-}"; shift 2 ;;
    --) shift; break ;;
    *) break ;;
  esac
done
case "$ONLY" in
  both|dashboard|molstar) ;;
  *) die "--only takes dashboard, molstar or both (got '${ONLY}')." ;;
esac
# Pass-through flags for 02_render_dashboard.py. Expanded below as ${EXTRA[@]+"${EXTRA[@]}"}:
# "${EXTRA[@]-}" looks equivalent but expands an EMPTY array to one empty string, which sbatch
# forwards to the job script and argparse then rejects with `unrecognized arguments:`.
EXTRA=("$@")

# shellcheck source=../02.10.00_slurm_bundle/cluster_env.sh
source "$HERE/../02.10.00_slurm_bundle/cluster_env.sh" || exit 1
[ -n "${ZH_ACCOUNT:-}" ] || die "ZH_ACCOUNT is empty in ${ZH_CLUSTER_ENV}."
: "${ZH_CPU_TIME:=04:00:00}" ; : "${ZH_CPU_CPUS:=4}" ; : "${ZH_CPU_MEM:=16G}"
: "${ZH_CPU_PARTITION:=${ZH_GPU_PARTITION:-}}"
[ -n "$ZH_CPU_PARTITION" ] || die "set ZH_CPU_PARTITION (or ZH_GPU_PARTITION) in ${ZH_CLUSTER_ENV}."

LISTING="$(python "$HERE/02_render_dashboard.py" --list)" \
  || die "could not list the exported movies. Run ./submit_export.sh first."
N=$(printf '%s\n' "$LISTING" | grep -c .)
[ "$N" -gt 0 ] || die "no exported movies found. Run ./submit_export.sh first."

FMT='  %-3s %-12s %-9s %8s %10s %8s\n'
# shellcheck disable=SC2059 -- FMT is a fixed format string, not user input
printf "$FMT" "#" "system" "replica" "frames" "ns/frame" "atoms"
printf '%s\n' "$LISTING" | awk -F'\t' -v fmt="$FMT" '{ printf fmt, $1, $2, $3, $4, $5, $6 }'
echo "  ${N} exported replicas, passes: ${ONLY}"

# The two things that silently degrade the OUTPUT rather than failing the job, checked here on the
# login node so the answer arrives before the array does rather than in fifteen task logs.
command -v ffmpeg >/dev/null 2>&1 || {
  echo "WARNING: ffmpeg is not on PATH in this environment. Without it the dashboard movies are"
  echo "  written as animated GIFs -- several times larger, and not seekable."
  echo "  conda env update -f environment_zh853mor-prep.yml"; }
if [ "$ONLY" != "dashboard" ]; then
  command -v node >/dev/null 2>&1 || {
    echo "WARNING: node is not on PATH; the MolStar pass will be skipped by every task."
    echo "  conda env update -f environment_zh853mor-prep.yml"; }
  [ -d "$HERE/node_modules/puppeteer" ] || {
    echo "WARNING: node_modules is not installed; the MolStar pass will be skipped by every task."
    echo "  Run this on a LOGIN node (compute nodes usually have no outbound network):"
    echo "      cd ${HERE#"$(cd "$HERE/../.." && pwd)/"} && npm install"; }
fi

# SLURM writes --output relative to the submission directory, which is this stage directory under
# src/ -- where generated files do not belong (D-16). Keep the logs with the data this stage
# produces instead. The directory must exist before sbatch runs.
LOGS="$(cd "$HERE/../.." && pwd)/intermediate/02.12.00_render_trajectory_movies/logs"
[ "$DRY" -eq 1 ] || mkdir -p "$LOGS" || die "cannot create the log directory $LOGS."

SBATCH_ARGS=(
  --account="$ZH_ACCOUNT"
  --output="$LOGS/render_%A_%a.out"
  --partition="$ZH_CPU_PARTITION"
  --time="$ZH_CPU_TIME"
  --cpus-per-task="$ZH_CPU_CPUS"
  --mem="$ZH_CPU_MEM"
  --array="1-${N}"
  --export="ALL,ZH_CLUSTER_ENV=$ZH_CLUSTER_ENV,ZH_RENDER_ONLY=$ONLY"
)
[ -n "${ZH_QOS:-}" ]        && SBATCH_ARGS+=(--qos="$ZH_QOS")
[ -n "${ZH_CONSTRAINT:-}" ] && SBATCH_ARGS+=(--constraint="$ZH_CONSTRAINT")
if [ -n "${ZH_EXTRA_SBATCH:-}" ]; then
  # shellcheck disable=SC2206 -- deliberate word-splitting: ZH_EXTRA_SBATCH holds whole flags
  SBATCH_ARGS+=($ZH_EXTRA_SBATCH)
fi

if [ "$DRY" -eq 1 ]; then
  echo "DRY RUN: sbatch ${SBATCH_ARGS[*]} $HERE/submit_render.sbatch" \
       "${EXTRA[@]+${EXTRA[*]}}"
  exit 0
fi
jid=$(cd "$HERE" && sbatch --parsable "${SBATCH_ARGS[@]}" submit_render.sbatch \
        ${EXTRA[@]+"${EXTRA[@]}"})
echo "submitted movie rendering: $jid   ($N replicas, one per array task)"
echo "  -> product/02.12.00_<system>_<replica>_{dashboard,molstar}_<date>.mp4"
echo "  logs: intermediate/02.12.00_render_trajectory_movies/logs/render_${jid}_*.out"
