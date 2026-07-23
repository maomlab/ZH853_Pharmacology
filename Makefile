.PHONY: help env env-cluster lint format typecheck test check fetch qc clean-intermediate

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  %-18s %s\n", $$1, $$2}'

env:  ## Create the local analysis conda env
	conda env create -f environment.yml || conda env update -f environment.yml

env-cluster:  ## Create the SLURM/simulation conda env
	conda env create -f environment-cluster.yml || conda env update -f environment-cluster.yml

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

clean-intermediate:  ## Remove cached intermediate results
	rm -rf intermediate/*
