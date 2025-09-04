import os
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd

# --- Configuration ---
ref_file = "reference.pdb"   # single eference ligand withput protein
dock_folder = "GninaOut_autobox_FMN"      # Folder with single docked ligands (docking results)
output_csv = "ligand_com_distances.csv"

# --- Helper function to compute center of mass ---
def center_of_mass(mol):
    conf = mol.GetConformer()
    masses = np.array([atom.GetMass() for atom in mol.GetAtoms()]) # get masses
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]) # get coordinates
    com = np.average(coords, axis=0, weights=masses) # average masses over coordinates = COM
    return com

# --- Load reference ligand ---
ref = Chem.MolFromPDBFile(ref_file, removeHs=False) # load molecule from pdb as rdkit object
if ref is None:
    raise ValueError("Could not load reference ligand.")
ref_com = center_of_mass(ref)

# --- Collect results ---
results = []

for fname in os.listdir(dock_folder):
    if not fname.endswith(".pdb"):
        continue
    path = os.path.join(dock_folder, fname)
    
    dock = Chem.MolFromPDBFile(path, removeHs=False)
    
    dock_com = center_of_mass(dock) # center of mass
    distance = np.linalg.norm(dock_com - ref_com) # differences between reference and docked ligand
    
    # Extract ID from filename (strip "_docked.pdb")
    ligand_id = fname.replace("_docked.pdb", "")
    
    results.append({
        "ID": ligand_id,
        "distance_to_ref": distance
    })

# --- Export to CSV ---
df = pd.DataFrame(results)
df.to_csv(output_csv, index=False)
