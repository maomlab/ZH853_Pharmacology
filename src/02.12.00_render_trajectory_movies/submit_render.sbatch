#!/bin/bash
#SBATCH --job-name=zh853_render
#SBATCH --output=../../intermediate/02.12.00_render_trajectory_movies/logs/render_%A_%a.out
#
# The output path is relative to the SUBMISSION directory, i.e. to this stage directory, and
# points into intermediate/ so no log is written into src/ (D-16). submit_render.sh creates that
# directory and passes an absolute path; submitting this file BY HAND needs
#     mkdir -p ../../intermediate/02.12.00_render_trajectory_movies/logs
# first -- SLURM cannot start a job whose output file it cannot open.
#
# CPU stage: render ONE exported replica per array task. Rendering is single-core and serial per
# replica -- matplotlib draws 300 frames and pipes them to ffmpeg, MolStar draws 300 more in
# software WebGL -- so the way to make the panel fast is to run the replicas concurrently rather
# than to thread any one of them. Nothing here needs a GPU.
#
#     ./submit_render.sh
#
# The array index selects the (system, replica) pair from the SAME listing the submit script sized
# the array with -- 02_render_dashboard.py --list -- so the mapping cannot drift between the two.
#
# ZH_RENDER_ONLY (exported by submit_render.sh) restricts which pass runs: dashboard | molstar |
# both. The MolStar pass additionally needs node and an installed node_modules; when either is
# missing the task says so and still produces the dashboard movie, because a cluster image without
# Chromium's shared libraries is a normal outcome and not a reason to fail the array.

set -euo pipefail

# The site settings live at the repository root; this stage has no sampling.env of its own.
# shellcheck source=../02.10.00_slurm_bundle/cluster_env.sh
source "${SLURM_SUBMIT_DIR:-$PWD}/../02.10.00_slurm_bundle/cluster_env.sh" || exit 1

STAGE="${SLURM_SUBMIT_DIR:-$PWD}"
ONLY="${ZH_RENDER_ONLY:-both}"

module purge
if [ -n "${ZH_MODULES:-}" ]; then
  # shellcheck disable=SC2086 -- ZH_MODULES is a space-separated module list
  module load ${ZH_MODULES}
fi
source "${ZH_CONDA_SH:-$(conda info --base)/etc/profile.d/conda.sh}"
conda activate "${ZH_PREP_ENV:-zh853mor-prep}"   # MDAnalysis + matplotlib + ffmpeg live here

IDX="${SLURM_ARRAY_TASK_ID:-1}"
echo "host: $(hostname)  job: ${SLURM_JOB_ID:-none}  task ${IDX}  passes: ${ONLY}"

# Resolve this task's index to a (system, replica) pair ONCE, so both passes render the same one.
read -r SYSTEM REPLICA <<<"$(python "$STAGE/02_render_dashboard.py" --list \
  | awk -F'\t' -v i="$IDX" '$1 == i { print $2, $3 }')"
[ -n "${SYSTEM:-}" ] || { echo "ERROR: no exported movie at index ${IDX}" >&2; exit 1; }
echo "rendering ${SYSTEM}/${REPLICA}"

# ffmpeg is what separates an mp4 from a silent fallback to a much larger, unseekable GIF. It is
# declared in environment_zh853mor-prep.yml; if it is missing the env predates that and needs
# `conda env update -f environment_zh853mor-prep.yml`.
if ! command -v ffmpeg >/dev/null 2>&1; then
  echo "WARNING: ffmpeg is not on PATH in ${ZH_PREP_ENV:-zh853mor-prep}." >&2
  echo "  conda env update -f environment_zh853mor-prep.yml   # then resubmit" >&2
fi

status=0
if [ "$ONLY" = "both" ] || [ "$ONLY" = "dashboard" ]; then
  echo "--- dashboard"
  python "$STAGE/02_render_dashboard.py" --system "$SYSTEM" --replica "$REPLICA" "$@" || status=1
fi

if [ "$ONLY" = "both" ] || [ "$ONLY" = "molstar" ]; then
  echo "--- molstar"
  if ! command -v node >/dev/null 2>&1; then
    echo "SKIPPED: node is not on PATH. Add nodejs to the env"
    echo "  (conda env update -f environment_zh853mor-prep.yml), or render this pass locally."
  elif [ ! -d "$STAGE/node_modules/puppeteer" ]; then
    echo "SKIPPED: node_modules is not installed. On a LOGIN node (compute nodes usually have no"
    echo "  outbound network): cd src/02.12.00_render_trajectory_movies && npm install"
  else
    # Chromium is launched with --no-sandbox already; a cluster image may still lack the X/NSS
    # libraries it links against, which surfaces as a launch failure rather than a render fault.
    node "$STAGE/03_render_molstar.js" --system "$SYSTEM" --replica "$REPLICA" \
      || { echo "molstar pass failed (often Chromium's shared libraries); the dashboard movie" \
                "for this replica was still written" >&2; status=1; }
  fi
fi
exit "$status"
