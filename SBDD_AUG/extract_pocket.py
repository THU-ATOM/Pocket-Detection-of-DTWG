import os
import sys
from tqdm import tqdm
from rdkit import Chem
import numpy as np
from multiprocessing import Pool
from Bio.PDB import PDBParser, Selection, PDBIO
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.StructureBuilder import StructureBuilder
from Drug_The_Whole_Genome.utils.Dataset import DTWG_DUDEDataset
import glob
from copy import deepcopy
from Drug_The_Whole_Genome.utils.utils import pdb2dict,write_lmdb,save_pdb
class PocketExtractor():
    
    def __init__(self):
        pass

    def run(self, tasks):
        # task: {protein:"", ligand:"", threshold:"",output:""}
        # Use pool.map to process tasks in parallel
        task_list=[(x['protein'],x['ligand'],x['threshold'],x['output']) for x in tasks]
        with Pool(80) as p:
            results = p.starmap(self._extract_single, task_list)
        # for task in tqdm(task_list):
        #     self._extract_single(*task)
    
    def _extract_single(self,protein_file,ligand_file,threshold,output_file):

        if ligand_file is not None:
            # read ligand
            try:
                if ligand_file.endswith('.mol2'):
                    ligand = Chem.MolFromMol2File(ligand_file,sanitize=False)
                elif ligand_file.endswith('.pdb'):
                    ligand = Chem.MolFromPDBFile(ligand_file,sanitize=False)
                elif ligand_file.endswith('.sdf'):
                    ligand = Chem.MolFromMolFile(ligand_file,sanitize=False)
                else:
                    raise NotImplementedError
                assert ligand is not None
            except:
                print(f"Failed to read ligand {ligand_file}")
                return
                
            conf = ligand.GetConformer()
            ligand_coords = conf.GetPositions()
            
            # read protein
            try:
                protein = PDBParser(QUIET=True).get_structure("protein",protein_file)[0]
                assert protein is not None
            except:
                print(f"Failed to read protein {protein_file}")
                return
        else:
            # read protein
            try:
                protein = PDBParser(QUIET=True).get_structure("protein",protein_file)[0]
                assert protein is not None
            except:
                print(f"Failed to read protein {protein_file}")
                return

            ligand_chain = protein['B']
            ligand_coords = []
            for residue in ligand_chain:
                for atom in residue:
                    ligand_coords.append(atom.coord)

            for chain in protein:
                if chain.id != 'A':
                    protein.detach_child(chain.id)


   
        # extract pocket
        for chain in protein:
            remove_residue_ids=[]
            for residue in chain:
                f=1
                for atom in residue:
                    protein_atom_coords = np.array(atom.coord)
                    for ligand_coord in ligand_coords:
                        if np.linalg.norm(protein_atom_coords - ligand_coord) < threshold:
                            f=0
                            break
                if f:
                    remove_residue_ids.append(residue.id)
            for residue_id in remove_residue_ids:
                chain.detach_child(residue_id)
            
        self.pocket = protein
        # remove all atoms with element X or H
        self._remove_atom_element_X_H()

        # save pocket
        io = PDBIO()
        io.set_structure(protein)
        io.save(output_file)
    
    def  _remove_atom_element_X_H(self):
        for chain in self.pocket:
            for residue in chain:
                remove_atom_ids = []
                for atom in residue:
                    if atom.element == 'X' :
                        print(f"Remove atom element X: {atom}")
                        remove_atom_ids.append(atom.id)
                    if atom.element == 'H' :
                        print(f"Remove atom element H: {atom}")
                        remove_atom_ids.append(atom.id)
                for atom_id in remove_atom_ids:
                    residue.detach_child(atom_id)

def save_docking_data(complex_file,output_dir):
        
    # split protein into receptor.pdb and ligand.sdf from complex_file
    parser = PDBParser()
    structure = parser.get_structure('complex', complex_file)[0]
    receptor=structure['A']
    save_pdb(receptor,os.path.join(output_dir,'receptor.pdb'))
    ligand=structure['B']
    save_pdb(ligand,os.path.join(output_dir,'ligand.pdb'))
    # convert ligand to sdf
    cmd=f'obabel {os.path.join(output_dir,"ligand.pdb")} -O {os.path.join(output_dir,"ligand.sdf")}'
    os.system(cmd)
    #rm ligand.pdb
    os.remove(os.path.join(output_dir,"ligand.pdb"))
    if os.path.exists(os.path.join(output_dir,"ligand.sdf")):
        return True
    else:
        raise Exception


if __name__ == "__main__":
    dataset=DTWG_DUDEDataset()
    exp_name="AF2_fpocket_BFN"
    
    # extract pocket
    tasks=[]
    for item in tqdm(dataset):
        complexes=glob.glob(os.path.join(item['dir'],exp_name,'*complex*.pdb'))
        complexes=[x for x in complexes if 'pocket' not in os.path.basename(x)]
        complexes.sort()
        for complex in complexes:
            output_file=complex.replace('.pdb','_pocket6A.pdb')
            tasks.append({
                'protein':complex,
                'ligand':None,
                'threshold':6,
                'output':output_file
            })

    pocket_extractor=PocketExtractor()
    pocket_extractor.run(tasks)

    
    # generate lmdb
    for item in tqdm(dataset.get_items()):
        if not os.path.exists(item['protein_dir']):
            continue
        pockets=glob.glob(os.path.join(item['dir'],exp_name,'*pocket*.pdb'))
        pockets=[x for x in pockets if 'pocket' in os.path.basename(x)]
        pockets.sort() 
        parser = PDBParser()
        output_path = os.path.join(item['dir'],exp_name,"pockets.lmdb")

        dics=[]
        for pocket in pockets:
            structure = parser.get_structure('pocket', pocket)
            dic=pdb2dict(structure,pocket)
            dics.append(dic)
        write_lmdb(dics,output_path)


    # generate docking_data

    total_cnt=0
    for item in tqdm(dataset.get_items()):
        if not os.path.exists(item['protein_dir']):
            continue
        pockets=glob.glob(os.path.join(item['dir'],exp_name,'*complex*.pdb'))
        pockets=[x for x in pockets if 'pocket' not in os.path.basename(x)]
        pockets.sort()
        parser = PDBParser()
        output_path = os.path.join(item['dir'],exp_name,"docking_data")
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        for pocket in pockets:
            output_dir=os.path.join(output_path,os.path.basename(pocket).replace('.pdb',''))
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            res=save_docking_data(pocket,output_dir)
            if res:
                total_cnt+=1
            print(f"Total: {total_cnt}")
