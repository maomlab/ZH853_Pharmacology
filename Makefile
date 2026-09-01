# ==============================================================================
# ZH853-MOR analysis. Environments (which interpreter each target needs):
#
#   LOCAL analysis env  -- EVERY target in this Makefile runs here.
#                          Create with `make env`, then `conda activate zh853mor`.
#                          (`molstar-render` additionally needs Node.js >= 18.)
#
#   CLUSTER envs (zh853mor-prep / -sim / -plumed) -- used ONLY by the SLURM bundle
#                          in src/02.10.00_slurm_bundle/, run ON the cluster, NOT via
#                          make. See that directory's README. `make env-cluster`
#                          just creates those specs if you are on the cluster.
#
# Dependencies are encoded as prerequisites, so `make <target>` first runs whatever
# it needs (e.g. `make interactions` runs `fetch`; `make molstar-render` runs
# `prep-assess`). Prereqs are cheap/idempotent (fetch skips existing files).
# Run `make help` for the grouped target list.
# ==============================================================================

.PHONY: help \
        env env-cluster env-plumed \
        lint format typecheck test check \
        fetch \
        qc interactions mutations analogs design interaction-map depictions analysis \
        prep-assess prep-protonation prep-receptor prep-orient prep-ligand prep-analogs \
        membrane-plot prep \
        molstar-render figures manuscript clean-intermediate

help:  ## Show this grouped target list
	@awk 'BEGIN {FS = ":.*?## "} \
		/^## / {printf "\n\033[1m%s\033[0m\n", substr($$0, 4)} \
		/^[a-zA-Z0-9_-]+:.*?## / {printf "  \033[36m%-18s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

## Environments
env:  ## Create the LOCAL analysis env (this Makefile runs here)
	conda env create -f environment.yml || conda env update -f environment.yml

env-cluster:  ## Create the cluster prep+sim envs (for the SLURM bundle, on the cluster)
	conda env create -f environment-prep.yml    || conda env update -f environment-prep.yml
	conda env create -f environment-cluster.yml || conda env update -f environment-cluster.yml

env-plumed:  ## Create the cluster metadynamics env (optional; openmm-plumed)
	conda env create -f environment-plumed.yml || conda env update -f environment-plumed.yml

## Development & CI  [local env]
lint:  ## Ruff lint (package + tests)
	ruff check src/zh853mor tests

format:  ## Ruff auto-format
	ruff format src tests
	ruff check --fix src/zh853mor tests

typecheck:  ## Mypy type-check
	mypy

test:  ## Run unit tests
	pytest

check: lint typecheck test  ## Lint + typecheck + test (the CI gate)

## Data acquisition  [local env]
fetch:  ## Download comparator PDBs -> data/comparators/ (prereq of: interactions, mutations)
	python src/01.01.00_fetch_comparators.py

## Static analysis - Objectives 1-3  [local env]
qc:  ## Structure QC report -> product/
	python src/01.02.00_qc_structure.py

interactions: fetch  ## Comparative interaction fingerprints + heatmap -> product/ (Obj 1)
	python src/03.01.00_interaction_fingerprints.py

mutations: fetch  ## Ranked ZH853-selective mutation panel -> product/ (Obj 2)
	python src/03.02.00_mutation_panel.py

analogs:  ## Analog physicochemical property panel -> product/ (Obj 3)
	python src/05.01.00_analog_properties.py

design:  ## Structure-guided modification design -> product/ (Obj 3)
	python src/05.02.00_design_modifications.py

interaction-map:  ## PoseView-style 2D interaction map -> manuscript fig5 (Obj 1)
	python src/03.03.00_interaction_map.py

depictions:  ## 2D vector ligand depictions -> manuscript fig4 (Obj 3)
	python src/05.03.00_ligand_depictions.py

analysis: qc interactions mutations analogs design  ## Run the full static-analysis pipeline

## MD system prep - Phase 2  [local env; outputs feed the SLURM bundle]
prep-assess:  ## Assess prep needs + split components -> product/, intermediate/ (prereq of: molstar-render)
	python src/02.01.00_assess_and_split.py

prep-protonation:  ## PROPKA protonation states -> product/
	python src/02.02.00_protonation.py

prep-receptor:  ## Rebuild receptor sidechains (PDBFixer) -> intermediate/
	python src/02.03.00_prepare_receptor.py

prep-orient: prep-receptor  ## Orient receptor to membrane normal (z) for PACKMOL-Memgen --preoriented
	python src/02.05.00_orient_receptor.py

prep-ligand:  ## Protonated ZH853 + parameterization inputs -> intermediate/
	python src/02.04.00_ligand_prep.py

prep-analogs: prep-orient prep-ligand  ## ZH850/ZH831/ZH809 poses by scaffold transfer from ZH853
	python src/02.07.00_analog_poses.py

membrane-plot: prep-orient  ## Membrane-placement determination plot -> manuscript (panel B)
	python src/02.06.00_membrane_placement.py

prep: prep-assess prep-protonation prep-receptor prep-orient prep-ligand prep-analogs  ## Run the full Phase-2 local prep

## Figures & manuscript  [local env; molstar-render also needs Node.js >= 18]
molstar-render: prep-assess prep-orient  ## Headless MolStar 3D renders -> manuscript (complex, pocket, membrane; needs Node.js)
	cd src/03.10.00_molstar_render && npm install --silent \
	  && python build_overview_camera.py && python build_pocket_mvs.py && python build_membrane_mvs.py \
	  && node render.js && python trim_figures.py

figures: interactions analogs design depictions interaction-map membrane-plot molstar-render  ## Regenerate + stage ALL manuscript figures
	cp $$(ls -t product/03.01.00_fingerprint_heatmap_*.png  | head -1) product/manuscript/figures/fig1_interaction_heatmap.png
	cp $$(ls -t product/05.01.00_analog_property_space_*.png | head -1) product/manuscript/figures/fig2_property_space.png
	cp $$(ls -t product/05.02.00_design_property_shifts_*.png | head -1) product/manuscript/figures/fig3_design_shifts.png

manuscript:  ## Compile the LaTeX -> manuscript.pdf (run `make figures` first if figures changed)
	cd product/manuscript && tectonic manuscript.tex

## Housekeeping
clean-intermediate:  ## Remove cached intermediate results
	rm -rf intermediate/*
