#!/bin/bash
# Build the conda environments on the cluster. Run once on a login node.
# Three task-specific envs (they have mutually incompatible openmm pins, so they must be split):
#   zh853mor-prep    CPU  -- system building + ligand params (AmberTools)
#   zh853mor-sim     GPU  -- equilibration, production, free energy
#   zh853mor-plumed  GPU  -- metadynamics only (openmm-plumed; older openmm)  [optional]
set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
REPO="$(cd "$(dirname "$0")/../.." && pwd)"

# Faster solves if available (classic conda is slow); harmless if already present.
conda install -n base -y -c conda-forge mamba >/dev/null 2>&1 || true

conda env create -f "$REPO/environment_zh853mor-prep.yml"    || conda env update -f "$REPO/environment_zh853mor-prep.yml"
conda env create -f "$REPO/environment_zh853mor-sim.yml"  || conda env update -f "$REPO/environment_zh853mor-sim.yml"
# Metadynamics only -- uncomment when you reach Methods 3.9:
# conda env create -f "$REPO/environment_zh853mor-plumed.yml" || conda env update -f "$REPO/environment_zh853mor-plumed.yml"

# Verify the GPU run env sees CUDA. NOTE: run this on a GPU node (srun --gres=gpu:1 ... --pty bash)
# -- on a login node only CPU/Reference platforms appear.
conda activate zh853mor-sim
python -m openmm.testInstallation || true
echo
echo "If 'CUDA' is not listed above, you are probably on a login node -- re-run"
echo "  python -m openmm.testInstallation   on a GPU node. cuDNN is NOT needed for OpenMM."
