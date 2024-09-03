from rdkit import Chem
from rdkit.Chem import AllChem
import os
from Drug_The_Whole_Genome.utils.Dataset import DTWG_DUDEDataset
from tqdm import tqdm
import glob
from Bio.PDB import PDBParser, PDBIO
from rdkit.Chem import SDMolSupplier, MolToPDBBlock



if __name__ == "__main__":

    # relax pdb 
    dataset=DTWG_DUDEDataset()
    exp_name="RealProtein_BFNPocket"
    
    for item in tqdm(dataset):
        protein_file=item['protein_dir']
        # read by Bio
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', protein_file)[0]
        # remove all atoms with element X
        for chain in structure:
            for residue in chain:
                remove_atom_ids = []
                for atom in residue:
                    if atom.element == 'X':
                        print(f"Remove atom element X: {atom}")
                        remove_atom_ids.append(atom.id)
                for atom_id in remove_atom_ids:
                    residue.detach_child(atom_id)
        # save pocket
        io = PDBIO()
        io.set_structure(structure)
        io.save(protein_file.replace('.pdb','_noX.pdb'))