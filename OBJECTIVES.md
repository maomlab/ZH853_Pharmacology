# ZH853 MOR Docking
The objective of this project is develop a better characterization of
the structure-function relationship for how ZH853 binds to and
modulates the mu opioid receptor. We have used CryoEM to
experimentally characterize the mu opioid recptor in complex with
ZH854 and the Gi g-protein and scfv16, a single-chain variable
fragment derived from mAb16 (Maeda, et al., 2018, DOI:
10.1038/s41467-018-06002-w). The structure is
`data/mor_gi_scfv16_refine-coot-22_real_space_refined_169_edit.pdb`.

I would like to explore the following questions

  1) What are the unique structural interactions between the ZH853 and
  other experimentally characterized mu opioid receptor bound
  complexes.
	 
  2) What are the key molecular interactions that stabilize the ZH853
  binding mode? Specifically, what mutations should be selected to
  test that would abrogate the ZH853 binding/function but not disrupt
  binding/function of other related full agonists?
  
  3) What modifications could be made to ZH853 that would improve the
  drug-like character? Consider recent strategies for improving the
  pharmacokentics and pharmacodynamics of peptide-like drugs such 
  as GLP-1 modulators.
  
  4) Use free energy perturbation methods ABFE, OpenFF, or
  c(t)-based metadynamics (CTMD) to predict the change in affinity
  among ZH853 and analogs:
  
Endomorphin Analog	(Zadina, 2016, 10.1016/j.neuropharm.2015.12.024)		ZH850	Analog 1	Tyr-c[D-Lys-Trp-Phe-Glu]-NH2	NC([C@@H]1CCC(NCCCC[C@H](C(N[C@H](C(N[C@H](C(N1)=O)Cc2ccccc2)=O)Cc(c[nH]3)c4c3cccc4)=O)NC([C@H](Cc5ccc(O)cc5)N)=O)=O)=O
Endomorphin Analog	(Zadina, 2016, 10.1016/j.neuropharm.2015.12.024)		ZH831	Analog 2	Tyr-c[D-Glu-Phe-Phe-Lys]-NH2	NC(=O)[C@@H]4CCCCNC(=O)CC[C@@H](NC(=O)[C@@H](N)Cc1ccc(O)cc1)C(=O)N[C@@H](Cc2ccccc2)C(=O)N[C@@H](Cc3ccccc3)C(=O)N4
Endomorphin Analog	(Zadina, 2016, 10.1016/j.neuropharm.2015.12.024)		ZH809	Analog 3	Tyr-c-[D-Lys-Trp-Phe-Asp]-NH2	NC(=O)[C@@H]5CC(=O)NCCCC[C@@H](NC(=O)[C@@H](N)Cc1ccc(O)cc1)C(=O)N[C@@H](Cc2c[nH]c3ccccc23)C(=O)N[C@@H](Cc4ccccc4)C(=O)N5
Endomorphin Analog	(Zadina, 2016, 10.1016/j.neuropharm.2015.12.024)		ZH853	Analog 4	(Tyr-[D-Lys-Phe-Phe-Asp]2-NH2)2	NC(=O)CNC(=O)[C@@H]5CCC(=O)NCCCC[C@@H](NC(=O)[C@@H](N)Cc1ccc(O)cc1)C(=O)N[C@@H](Cc2c[nH]c3ccccc23)C(=O)N[C@@H](Cc4ccccc4)C(=O)N5	 

To analyze the structure use OpenMM to run molecular dynamics on the
receptor and ligand complex.  To prepare the structure, use best
practices, and compare against relavant simulations of related GPCR
systems.

To analyze molecular interactions use the Biotite, RDKit, and PLIP.

To actually run computationally intensive simulations, I will manually
run them on an academic SLURM based cluster. So please prepare clear,
self-contained simulation tasks that contian the relevant data files,
submission sbatch files and workflow scripts for installing,
preparing, submitting, gathering and post-processing.

To organize the analyses, please use the following project organizing principles

  1) The analysis should take organize raw/external input data in the
  `data/`
  
  2) Workflow and supporting code should go into the `src/` directory
  where for reproducibility, workflow scripts should be indexed using
  the format `src/##.##.##_file_or_folder_name`, so scripts sort
  sequentially and for a hierarchical organization so scripts can be
  added or modified without renumbering all of them.
  
  3) Temporary data and results can be cached into the intermediate
  directory, using naming scheme for the script that generated them
  where approrpate.
  
  4) Output results, figures, tables, etc. should be put into the
  `product/` folder. Use _YYYYMMDD.<extension> date-code format
  to do a simple form of version control. Also organize the results
  using the same naming scheme for the scripts in src and data in
  the intermediate data directory.
  
For the full analysis generate high-quality summary plots and
quantitative analysis summaries. Report on quality-control
tests to quantify simulation quality. Investigate scientificly
relevant hypotheses and test modeling assumptions. Combine
the analyses into a well contructed manuscript with consice
and clear introduction and related works, methods, results,
and conclusions. 

Use version control and coding best practices in writing code
including linting and type checking. Use a Makefile to run
common analysis and housekeeping workflows.

Create a clear README, and maintain in the docs folder a detailed plan
of what has been accomplished, what is planned, and what has been
considered but not persued.

As the analysis progresses, if there is ambiguity or decisions need to
be made, record these decisions and clarifications in the
SPECIFICATION.md file

