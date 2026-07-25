.PHONY: help env env-cluster lint format typecheck test check fetch qc clean-intermediate

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  %-18s %s\n", $$1, $$2}'

env:  ## Create the local analysis conda env
	conda env create -f environment.yml || conda env update -f environment.yml

env-cluster:  ## Create the cluster conda envs (prep + sim; run on the cluster)
	conda env create -f environment-prep.yml    || conda env update -f environment-prep.yml
	conda env create -f environment-cluster.yml || conda env update -f environment-cluster.yml

env-plumed:  ## Create the metadynamics conda env (optional; openmm-plumed)
	conda env create -f environment-plumed.yml || conda env update -f environment-plumed.yml

lint:  ## Ruff lint (package + tests)
	ruff check src/zh853mor tests

format:  ## Ruff auto-format
	ruff format src tests
	ruff check --fix src/zh853mor tests

typecheck:  ## Mypy type-check
	mypy

test:  ## Run unit tests
	pytest

check: lint typecheck test  ## Lint + typecheck + test (Phase-0 gate)

fetch:  ## Download comparator PDBs into data/comparators/
	python src/01.01.00_fetch_comparators.py

qc:  ## Structure QC report -> product/
	python src/01.02.00_qc_structure.py

interactions:  ## Comparative interaction fingerprints + heatmap -> product/ (Objective 1)
	python src/03.01.00_interaction_fingerprints.py

mutations:  ## Ranked ZH853-selective mutation panel -> product/ (Objective 2)
	python src/03.02.00_mutation_panel.py

analogs:  ## Analog physicochemical property panel -> product/ (Objective 3)
	python src/05.01.00_analog_properties.py

design:  ## Structure-guided modification design -> product/ (Objective 3)
	python src/05.02.00_design_modifications.py

analysis: qc interactions mutations analogs design  ## Run the full static-analysis pipeline

prep-assess:  ## Phase 2: assess prep needs + split components -> product/, intermediate/
	python src/02.01.00_assess_and_split.py

prep-protonation:  ## Phase 2: PROPKA protonation states -> product/
	python src/02.02.00_protonation.py

prep-receptor:  ## Phase 2: rebuild receptor sidechains (PDBFixer) -> intermediate/
	python src/02.03.00_prepare_receptor.py

prep-ligand:  ## Phase 2: protonated ZH853 + parameterization inputs -> intermediate/
	python src/02.04.00_ligand_prep.py

prep: prep-assess prep-protonation prep-receptor prep-ligand  ## Run the full Phase-2 local prep

depictions:  ## 2D vector ligand depictions -> product/ + manuscript fig
	python src/05.03.00_ligand_depictions.py

interaction-map:  ## PoseView-style 2D interaction map -> product/ + manuscript fig
	python src/03.03.00_interaction_map.py

molstar-render:  ## Headless MolStar 3D renders (needs npm install first) -> product/ + manuscript figs
	cd src/03.10.00_molstar_render && npm install --silent \
	  && python build_overview_camera.py && python build_pocket_mvs.py && node render.js \
	  && python trim_figures.py

figures: depictions interaction-map  ## Regenerate the vector figures (MolStar via molstar-render)

manuscript: figures  ## Compile the LaTeX manuscript -> product/manuscript/manuscript.pdf
	cd product/manuscript && tectonic manuscript.tex

clean-intermediate:  ## Remove cached intermediate results
	rm -rf intermediate/*
