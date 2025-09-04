from CreateLig3D import write_smiles_sdf
from preprocess_prot_file import extract_protein

import os
import pandas as pd
import subprocess
from multiprocessing import Pool
from tqdm import tqdm
from functools import partial
import sys
import time


def create_docking_command(protein_file, ligand_file, lig_name, gnina_out_dir, autobox_ligand_file, container_path):
    cmd = (
        f"gnina --cpu=4 --receptor {protein_file}"
        f"--ligand {ligand_file}"
        f"--cnn_scoring=none"
        f"--autobox_ligand {autobox_ligand_file}"
        f"--autobox_add 4" # instead of docking center and size
        f"--out '{gnina_out_dir}/{lig_name}_docked.pdb'"
        f"--log '{gnina_out_dir}/{lig_name}_docking.log'"
    )
    return f"singularity exec {container_path} {cmd}"


def dock_single_ligand(lig, protein_file, out_dir_ligs, gnina_out_dir, autobox_ligand_file, container_path):
    lig_name, smiles = lig
    sdf_path = os.path.join(out_dir_ligs, f"{lig_name}.sdf")

    try:  # 3D generation
        write_smiles_sdf(smiles, sdf_path)
    except Exception as e:
        print(f"[ERROR] Failed 3D gen for {lig_name}: {e}", flush=True)
        return

    cmd = create_docking_command(
        protein_file, sdf_path, lig_name,
        gnina_out_dir, autobox_ligand_file, container_path
    )

    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] GNINA failed for {lig_name}: {e}", flush=True)
    except Exception as e:
        print(f"[ERROR] Unknown error for {lig_name}: {e}", flush=True)


def run_parallel_docking(list_of_ligs, protein_file, out_dir_ligs, gnina_out_dir, autobox_ligand_file, container_path, num_workers=4):
    with Pool(processes=num_workers) as pool:
        work = partial(
            dock_single_ligand,
            protein_file=protein_file,
            out_dir_ligs=out_dir_ligs,
            gnina_out_dir=gnina_out_dir,
            autobox_ligand_file=autobox_ligand_file,
            container_path=container_path
        )

        list(tqdm(pool.imap_unordered(work, list_of_ligs),
                  total=len(list_of_ligs),
                  desc="Docking ligands",
                  file=sys.stderr))


if __name__ == '__main__':
    print('Starting docking script...\n')
    cluster = True

    if cluster:
        input_file = '1QZU_prediction_FMN_substrate.pdb'
        ligand_smiles_file = 'filtered_input.csv'
        container_path = '/group/bioinf/Software/singularity/gnina_d236cb7.sif'
        output_dir = 'output directory' # insert output directory
    else:
        input_file = '/Volumes/bioinf/Users/Masine/GNINA_Docking_Data/Input/1QZU_sub_FMN.pdb'
        ligand_smiles_file = '/Volumes/bioinf/Users/Masine/GNINA_Docking_Data/Input/filtered_output_cleaned, cutoff = 0.7.csv'
        container_path = '/Volumes/bioinf/Software/singularity/gnina_d236cb7.sif'
        output_dir = '/Volumes/bioinf/Users/Masine/GNINA_Docking_Data'

    native_lig = 'l01' # substrate
    lig_chain = 'B'
    structure_pdb = '1QZU'
    lig_lib_name = '1QZU_FMN_autobox'
    restart = False

    # Prepare folders
    lig_coord_files = os.path.join(output_dir, f'lig3DFiles_{lig_lib_name}') # for sdf
    os.makedirs(lig_coord_files, exist_ok=True)

    gnina_out_dir = os.path.join(output_dir, f'GninaOut_autobox_FMN_{structure_pdb}_{native_lig}_{lig_chain}') # main output -> pdb and log
    os.makedirs(gnina_out_dir, exist_ok=True)

    # Clean protein file
    cleaned_protein = os.path.join(gnina_out_dir, f'{structure_pdb}_cleaned.pdb')
    extract_protein(input_file, os.path.join(gnina_out_dir, structure_pdb + '_cleaned.pdb'), chain_id=None)

    # Set autobox reference ligand (must be present in protein PDB)
    autobox_ligand_file = input_file  # assumes ligand l01 is present in 1QZU_sub_FMN.pdb

    # Load ligands
    ligs_df = pd.read_csv(ligand_smiles_file).drop_duplicates(subset=['ID'])
    if restart:
        old_files = [f.replace('.sdf', '') for f in os.listdir(lig_coord_files)]
        ligs_df = ligs_df[~ligs_df["ID"].isin(old_files)]

    list_of_ligs = list(ligs_df[['ID', 'SMILES']].itertuples(index=False, name=None))

    # Start docking
    start = time.time()
    run_parallel_docking(
        list_of_ligs,
        protein_file=cleaned_protein,
        out_dir_ligs=lig_coord_files,
        gnina_out_dir=gnina_out_dir,
        autobox_ligand_file=autobox_ligand_file,
        container_path=container_path,
        num_workers=32
    )
     print('\n'*6, '#'*40, 'end', '#'*40, '\n', f'Docking {len(ligs)} took: {time.time() - start} seconds')
