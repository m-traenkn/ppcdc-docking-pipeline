from Bio.PDB import PDBParser, PDBIO, Select

class ProteinAndNativeLigandSelect(Select):
    def __init__(self, chain_id=None, native_ligand_resname="l02"): # keep FMN
        self.chain_id = chain_id
        self.native_ligand_resname = native_ligand_resname

    def accept_chain(self, chain):
        return self.chain_id is None or chain.id == self.chain_id

    def accept_residue(self, residue):
        hetfield, _, _ = residue.get_id()
        resname = residue.get_resname()

        # Keep protein residues (hetfield == ' ') and the native ligand
        return hetfield == ' ' or resname == self.native_ligand_resname

    def accept_atom(self, atom):
        return True

def extract_protein(input_pdb, output_pdb, chain_id=None, native_ligand_resname="l02"): # keep FMN
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, select=ProteinAndNativeLigandSelect(chain_id, native_ligand_resname))

if __name__ == '__main__':
    import os
    pdb_id = '3mxf'
    gnina_out_dir = '/Volumes/bioinf/Users/Philipp/BioClusterDocking/test_files'
    extract_protein('/Volumes/bioinf/Users/Philipp/BioClusterDocking/test_files/3mxf.pdb', os.path.join(gnina_out_dir, pdb_id + '_cleaned.pdb'))