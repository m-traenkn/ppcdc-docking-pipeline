# Imports
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom
from rdkit.Geometry import Point3D

def write_smiles_sdf(Small_Mol_Smiles, out_file):
  # Input SMILES string
  # Convert to molecule
  mol = Chem.MolFromSmiles(Small_Mol_Smiles)
  mol = Chem.AddHs(mol)  # Add hydrogens

  # Generate 3D coordinates
  AllChem.EmbedMolecule(mol, AllChem.ETKDG())  # You can try AllChem.ETKDGv3() too
  AllChem.UFFOptimizeMolecule(mol)

  # Write to SDF file
  w = Chem.SDWriter(out_file)
  w.write(mol)
  w.close()

def move_ligand_to_center(input_sdf, output_sdf, box_center):
    """
    Adjust the ligand's coordinates to move it to the center of the docking box.

    Parameters:
        input_sdf (str): Path to the input .sdf file containing the ligand.
        output_sdf (str): Path to the output .sdf file for the adjusted ligand.
        box_center (tuple): Coordinates of the docking box center (x, y, z).
    """
    # Load the ligand from the SDF file
    ligand = Chem.SDMolSupplier(input_sdf, removeHs=False)[0]  # Load the first molecule
    if not ligand:
      raise ValueError("Could not load ligand from SDF file")

    # Get the current coordinates of the ligand
    conf = ligand.GetConformer()
    coords = [conf.GetAtomPosition(i) for i in range(ligand.GetNumAtoms())]

    # Calculate the geometric center of the ligand
    ligand_center = [
      sum(coord.x for coord in coords) / len(coords),
      sum(coord.y for coord in coords) / len(coords),
      sum(coord.z for coord in coords) / len(coords)
    ]

    #print("Original ligand center:", ligand_center)

    # Calculate the translation vector to move the ligand to the box center
    translation_vector = [
      box_center[0] - ligand_center[0],
      box_center[1] - ligand_center[1],
      box_center[2] - ligand_center[2]
    ]

    #print("Translation vector:", translation_vector)

    # Apply the translation to each atom
    for i in range(ligand.GetNumAtoms()):
      pos = conf.GetAtomPosition(i)
      new_pos = Point3D(float(pos.x + translation_vector[0]),
                        float(pos.y + translation_vector[1]),
                        float(pos.z + translation_vector[2]))
      conf.SetAtomPosition(i, new_pos)

    # Write the adjusted ligand to the output SDF file
    writer = Chem.SDWriter(output_sdf)
    writer.write(ligand)
    writer.close()