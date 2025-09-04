import os

# model 2 is the model with best docking affinity
def extract_model_2(pdb_path):
    """Extracts lines between 'MODEL 2' and 'ENDMDL' (inclusive)."""
    lines = []
    in_model = False
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("MODEL") and " 2" in line:
                in_model = True
                continue
            if line.startswith("ENDMDL") and in_model:
                break
            if in_model: # get coordinates from pdb file
                if line.startswith(("ATOM", "HETATM")):
                    lines.append(line)
    return lines

def combine_protein_ligand(protein_pdb, ligand_pdb, output_pdb):
    # Read protein atoms
    with open(protein_pdb, 'r') as f:
        protein_lines = [line for line in f if line.startswith(('ATOM', 'HETATM'))] # protein coordinates

    # Extract MODEL 2 atoms from ligand
    ligand_lines = extract_model_2(ligand_pdb) # ligand coordinates

    # Write combined file
    with open(output_pdb, 'w') as out:
        for line in protein_lines:
            out.write(line)
        for line in ligand_lines:
            out.write(line)
        out.write('END\n')


def combine_all(protein_pdb, ligand_folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)

    for file in os.listdir(ligand_folder): # go through all pdb files of docked ligand
        if not file.endswith(".pdb"):
            continue

        ligand_path = os.path.join(ligand_folder, file)
        output_path = os.path.join(output_folder, f"{os.path.splitext(file)[0]}_complex.pdb") # ID_docked_complex.pdb
        combine_protein_ligand(protein_pdb, ligand_path, output_path)

# Example usage
combine_all(
    protein_pdb="GninaOut_autobox_FMN/1QZU_cleaned.pdb", # cleaned protein
    ligand_folder="GninaOut_autobox_FMN", # docking results
    output_folder="combined_complexes_complete" # folder with combined pdb files
)

