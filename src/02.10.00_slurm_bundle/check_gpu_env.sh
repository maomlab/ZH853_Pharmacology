#!/bin/bash
#SBATCH --job-name=zh853_gpucheck
#SBATCH --output=gpucheck_%j.out
#SBATCH --gres=gpu:1               # TODO(OQ-3): GPU type, e.g. gpu:a100:1
#SBATCH --partition=gpu            # TODO(OQ-3): partition/queue
#SBATCH --account=CHANGEME         # TODO(OQ-3): allocation/account
#SBATCH --time=00:05:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#
# Pre-flight diagnostic: confirms a GPU node + the zh853mor-sim env can actually run OpenMM on
# CUDA, BEFORE you spend a real equilibration/production job finding out it can't.
#   sbatch check_gpu_env.sh          # submit as a short GPU job
#   bash   check_gpu_env.sh          # or run directly on an interactive GPU node
# Intentionally does NOT `set -e`: every check runs so you see the full picture.

echo "=== ZH853 GPU environment check ==="
echo "host: $(hostname)    date: $(date)"
echo

echo "--- loaded modules ---"
module list 2>&1 || echo "(no module system)"
echo

echo "--- nvidia-smi (driver + max CUDA + GPU) ---"
if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi
else
  echo "nvidia-smi not found -- are you on a GPU node? (srun --gres=gpu:1 ... --pty bash)"
fi
echo

echo "--- conda env (zh853mor-sim) ---"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate zh853mor-sim || echo "WARNING: could not activate zh853mor-sim (run 00_install.sh?)"
echo "python: $(which python)"
conda list 2>/dev/null | grep -E "^(openmm|openmmforcefields|cuda-version|cudatoolkit) " || true
echo

echo "--- openmm.testInstallation ---"
python -m openmm.testInstallation || true
echo

echo "--- explicit CUDA run (200 steps on the GPU) ---"
python - <<'PY'
import sys
import openmm as mm
from openmm import Vec3, unit

names = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
print("OpenMM", mm.version.version, "| platforms:", names)
if "CUDA" not in names:
    print("RESULT: FAIL -- CUDA platform not available to OpenMM")
    sys.exit(3)

# minimal but real system: 216 Lennard-Jones particles on a grid, 200 dynamics steps
n = 216
system = mm.System()
nb = mm.NonbondedForce()
nb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
for _ in range(n):
    system.addParticle(12.0)
    nb.addParticle(0.0, 0.34, 0.4)
system.addForce(nb)
integ = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 2 * unit.femtoseconds)
plat = mm.Platform.getPlatformByName("CUDA")
try:
    ctx = mm.Context(system, integ, plat, {"Precision": "mixed"})
    pos = [Vec3(0.4 * (i % 6), 0.4 * ((i // 6) % 6), 0.4 * (i // 36)) for i in range(n)] * unit.nanometer
    ctx.setPositions(pos)
    integ.step(200)
    energy = ctx.getState(getEnergy=True).getPotentialEnergy()
except Exception as exc:  # noqa: BLE001 -- diagnostic wants the full failure text
    print(f"RESULT: FAIL -- CUDA context/run error: {exc}")
    sys.exit(4)
try:
    print("CUDA device:", plat.getPropertyValue(ctx, "DeviceName"))
except Exception:
    pass
print(f"200 steps OK, potential energy = {energy}")
print("RESULT: PASS -- OpenMM runs on CUDA")
PY
rc=$?

echo
if [ "$rc" -eq 0 ]; then
  echo "==> PASS: this env can run the equilibration/production sbatch on this GPU."
else
  echo "==> FAIL: CUDA is not usable here. Checklist:"
  echo "    - Are you on a GPU node? nvidia-smi above must show a GPU."
  echo "    - cuDNN is NOT needed for OpenMM -- do not add it."
  echo "    - Load a CUDA module <= the nvidia-smi 'CUDA Version', or pin cuda-version in the env."
  echo "    - JIT compiler error? export OPENMM_CUDA_COMPILER=\$(which nvcc)"
fi
exit "$rc"
