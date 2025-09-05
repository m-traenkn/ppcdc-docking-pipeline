import os
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdMolDescriptors
import pandas as pd
from tqdm import tqdm

# --- Configuration ---
reference_file = "reference.pdb"  # single ligand pdb file
docked_folder = "GninaOut_autobox_FMN" # single ligand docking files
output_csv = "rmsd_results.csv"

# --- Load and prepare reference ligand ---
ref_raw = Chem.MolFromPDBFile(reference_file, removeHs=False) # convert pdb molecule to rdkit object
if ref_raw is None:
    raise ValueError("Could not load reference ligand.")

ref = Chem.AddHs(ref_raw) # add explicit hydrogens to reference
if not ref.GetConformer().Is3D():
    AllChem.EmbedMolecule(ref)
    AllChem.UFFOptimizeMolecule(ref)

results = []

# --- Process each docked ligand ---
for fname in sorted(os.listdir(docked_folder)):
    if not fname.endswith(".pdb"):
        continue

    #print(f"Processing: {fname}")
    path = os.path.join(docked_folder, fname)

    dock_raw = Chem.MolFromPDBFile(path, removeHs=False) # convert docked ligands
    if dock_raw is None:
        continue

    dock = Chem.AddHs(dock_raw) # add hydrogens to ligands
    if not dock.GetConformer().Is3D():
        AllChem.EmbedMolecule(dock)
        AllChem.UFFOptimizeMolecule(dock)

    # --- MCS Alignment ---
    mcs = rdFMCS.FindMCS([ref, dock],
                         completeRingsOnly=True, # only consider complete rings
                         ringMatchesRingOnly=True, # match rings to rings only
                         timeout=60)

    if mcs.canceled:
        print(f"[!] MCS search timed out for {fname}")
        continue

    pattern = Chem.MolFromSmarts(mcs.smartsString)  # Convert MCS SMARTS string to RDKit molecule pattern
    
    # Find atom indices in the reference and docked ligand corresponding to the MCS
    ref_match = ref.GetSubstructMatch(pattern)
    dock_match = dock.GetSubstructMatch(pattern)

    if not ref_match or not dock_match:
        print(f"[!] No MCS match found in {fname}")
        continue

    try:
        atom_map = list(zip(dock_match, ref_match)) # Create atom mapping between ref and dock for RMSD calculation
        
        # Align docked ligand to reference using MCS atoms and calculate RMSD
        rmsd = rdMolAlign.AlignMol(dock, ref, atomMap=atom_map)
        #print(f"MCS size: {len(atom_map)} atoms | RMSD: {rmsd:.3f}")
        
        ligand_id = fname.replace("_docked.pdb", "")
        
        results.append({
            "ID": ligand_id,
            "RMSD": rmsd, # mean distance between matched atoms
            "MCS": len(atom_map) # number of matched atoms
        })
    except Exception as e:
        print(f"[!] Alignment failed for {fname}: {e}")
        continue

# --- Save results ---
df = pd.DataFrame(results)
df.to_csv(output_csv, index=False)

