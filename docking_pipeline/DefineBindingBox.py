from Bio.PDB import PDBParser
import numpy as np


def compute_bounding_box(coords, padding=2.0):
    min_corner = coords.min(axis=0) - padding
    max_corner = coords.max(axis=0) + padding
    center = (min_corner + max_corner) / 2
    size = max_corner - min_corner
    return center, size

# extract residue position from PDB file to be used as binding box center

def extract_residue_coordinates(pdb_file, chain_id, res_number, padding = 2.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("target", pdb_file)

    coords = []
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    if residue.get_id()[1] == res_number:
                        for atom in residue:
                            coords.append(atom.coord)
                        break

    if not coords:
        raise ValueError(f"Residue {res_number} in chain {chain_id} not found.")

    return compute_bounding_box(np.array(coords), padding)

# extract center of mass of given ligand to be used as binding box center

def extract_ligand_center(pdb_file, chain_id, ligand_code, padding=2):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("target", pdb_file)

    coords = []

    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    hetfield, resseq, _ = residue.get_id()
                    if hetfield.startswith("H") and residue.resname.strip() == ligand_code:
                        for atom in residue:
                            coords.append(atom.coord)
                        break

    if not coords:
        raise ValueError(f"Ligand {ligand_code} in chain {chain_id} not found.")

    coords = np.array(coords)

    return compute_bounding_box(np.array(coords), padding)



if __name__ == '__main__':
    print('main')

    native_lig = 'l01'
    lig_chain = 'B'
    structure_pdb = '1QZU'
    input_file = '/Volumes/bioinf/Users/Philipp/GNINA_Docking_Data/MasineInput/1QZU_sub_FMN.pdb'

    center, size =extract_ligand_center(input_file, lig_chain, native_lig, padding=.1)

    print(center, size)

