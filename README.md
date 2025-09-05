# ppcdc-docking-pipeline
This repository contains the full computational pipeline developed for identifying candidate inhibitors of human PPCDC. 
It includes:   
- Protenix-based protein structure predictions
- Docking pipeline using GNINA with autobox setup
- PLIP-based scoring, ligand positioning and orientation
- Scripts for interaction analysis, ranking
- Enamine analog retrieval

Specific procedure:
1. filter and cluster molecules to get rid of redundancies -> *filter_and_cluster_molecules.py*
- Components-smiles-stereo-oe.csv to filtered_input.csv

2. docking -> *runDockingPar_Autobox.py* + run_docking_autobox.sh + all dependent scripts
- uses filtered_input.csv, 1QZU_prediction_FMN_substrate.pdb

3. combine cleaned protein with docked ligands: *combine_prot_lig.py*

4. post-processing
- filter affinity: *get_and_filter_affinities.py* -> extract affinities, filter molecules with higher affinity than substrate
- PLIP interactions: plp xmls, *plip_plots.py*, *plip_score.py* -> compares interactions to reference (1QZU, 1MVN, 6AIM), calculation of PLIP scores
- COM: *center_of_mass.py* -> calculated the COM distance between docked ligands and reference (needs single ligands without protein)
- MCS-RMSD: *compute_mcs_rmsd_ligands.py* -> computes MCS and then MCS-RMSD between reference (single substrate without protein) and docked ligand
- enamine: Predictions_mostSimENAMINE.csv
- ranking: *compress_and_rank.py* -> total ranking final_ligand_summary_complete.csv
- robustness of ranking: *ranking_robustness.py* -> goes through all combinations of criteria and calculates overlap and Spearson coefficient of the rankings
