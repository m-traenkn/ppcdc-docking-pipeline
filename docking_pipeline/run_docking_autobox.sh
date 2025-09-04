#!/bin/bash
# ----------------SLURM Parameters----------------
#SBATCH --partition all.q
#SBATCH --ntasks 128
#SBATCH --mem 256g
#SBATCH --nodes 1
#SBATCH --mail-type ALL
#SBATCH --job-name GNINA_Docking_autobox_FMN_sorted_0.7
#SBATCH --chdir /group/bioinf/Users/Masine/BioClusterDocking
#SBATCH --output /group/bioinf/Users/Masine/BioClusterDocking/run_logs/GNINA_Docking_autobox_FMN_sorted_0.7.out
#SBATCH --error /group/bioinf/Users/Masine/BioClusterDocking/run_logs/GNINA_Docking_autobox_FMN_sorted_0.7.err
#SBATCH --time 3:00:00

# ----------------Load Modules--------------------
module load apps/singularity
module load apps/anaconda/3.8

# ----------------Commands------------------------
conda run -p /group/bioinf/Users/Masine/conda_env/RDKIT_ENV python3 -u runDockingPar_Autobox.py