import os
from Bio.PDB import PDBParser, PDBIO, Select
from rdkit import Chem
from rdkit.Chem import AllChem

class LigandSelect(Select):
    def __init__(self, resname, chain_id):
        self.resname = resname
        self.chain_id = chain_id

    def accept_residue(self, residue):
        return (
            residue.get_resname().strip() == self.resname
            and residue.get_parent().id == self.chain_id
        )

def extract_ligand_to_sdf(pdb_path, ligand_code, chain_id, output_sdf_path):
    # Temporary path for the extracted ligand
    tmp_pdb_path = output_sdf_path.replace(".sdf", ".ligand.pdb")

    # 1. Extract the ligand using Biopython
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_path)

    io = PDBIO()
    io.set_structure(structure)
    io.save(tmp_pdb_path, LigandSelect(ligand_code.upper(), chain_id))

    # 2. Load into RDKit and write as .sdf
    mol = Chem.MolFromPDBFile(tmp_pdb_path, removeHs=False)
    if mol is None:
        raise ValueError("RDKit could not parse ligand from extracted PDB.")

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())  # Optional: 3D coordinates
    AllChem.UFFOptimizeMolecule(mol)             # Optional: optimize geometry

    w = Chem.SDWriter(output_sdf_path)
    w.write(mol)
    w.close()

    os.remove(tmp_pdb_path)
    return output_sdf_path


if __name__ == '__main__':


    input_file = '/Volumes/bioinf/Users/Philipp/GNINA_Docking_Data/MasineInput/1QZU_sub_FMN.pdb'
    native_lig = 'l01'
    lig_chain = 'B'
    structure_pdb = '1QZU'



    sdf_path = extract_ligand_to_sdf(
        pdb_path=input_file,
        ligand_code=native_lig,
        chain_id=lig_chain,
        output_sdf_path=f"/Volumes/bioinf/Users/Philipp/GNINA_Docking_Data/MasineInput/{native_lig}_LIG_{lig_chain}.sdf"
    )
    print("Saved ligand as:", sdf_path)