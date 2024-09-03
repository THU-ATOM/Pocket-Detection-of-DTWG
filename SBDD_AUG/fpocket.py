from Drug_The_Whole_Genome.utils.utils import read_mol2_ligand,merge_chains,get_binding_pockets,getfpockets,IoU,write_lmdb,pdb2dict,save_pdb
import os
import numpy as np
from tqdm import tqdm
from copy import deepcopy
from Bio.PDB import PDBParser
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning,BiopythonWarning
import shutil
import glob
from Drug_The_Whole_Genome.utils.Dataset import DTWG_DUDEDataset
from copy import deepcopy

warnings.filterwarnings(
    action='ignore',
    category=PDBConstructionWarning)
warnings.filterwarnings(
    action='ignore',
    category=BiopythonWarning)

def process_one_dude(home,pdbid,thres,output_file):
    # f = open('/home/jiayinjun/MultiTarget/tmp/maxiou.txt','a')
    workdir = os.path.join(home,pdbid)
    output_dir = None
    os.chdir(workdir)
    lig_mol2 = 'crystal_ligand.mol2'
    # pdbfile = 'receptor_noX.pdb'
    pdbfile = 'AF2_receptor.pdb'
    p = PDBParser()
    pdbstruct = p.get_structure(pdbid+'_mergedprot', os.path.join(workdir,pdbfile))
    processed_prot,_ = merge_chains(pdbstruct,False,workdir)
    if not os.path.isdir(os.path.join(workdir,pdbid+'_mergedprot_out/pockets')):
        cmd=f"fpocket -f {os.path.join(workdir,pdbid+'_mergedprot.pdb')}  > /dev/null 2>&1"
        os.system(cmd)
    mol_coord,_ = read_mol2_ligand(os.path.join(workdir,lig_mol2))

    processed_prot.id = pdbid+f'_{thres}A_ligpocket'
    _,res_ori,anum_ori = get_binding_pockets(processed_prot,mol_coord,thres,output_dir)
    pocpdb_list,res_fpoc,_,_,_ = getfpockets(deepcopy(processed_prot),os.path.join(workdir,pdbid+'_mergedprot_out/pockets'),thres,output_dir,return_score=True)
    max_iou = 0
    best_poc = None
    for r,p in zip(res_fpoc,pocpdb_list):
        if IoU(res_ori,r) > max_iou:
            max_iou = IoU(res_ori,r)
            best_poc = p
    if best_poc is not None:
        # f.write(f'{pdbid},{max_iou}\n')
        save_pdb(best_poc,output_file)
        # pdbdict = pdb2dict(best_poc,pdbid)
        # write_lmdb([pdbdict],os.path.join(workdir,'fpocket.lmdb'),1)

        no_sidechain_poc=deepcopy(best_poc)[0]
        for chain in no_sidechain_poc:
            for residue in chain:
                remove_atom_ids = []
                for atom in residue:
                    #preserve only C, CA, CB, N, O
                    if not atom.name in ['C','CA','CB','N','O']:
                        remove_atom_ids.append(atom.id)
                    if atom.element == 'X':
                        print(f"Remove atom element X: {atom}")
                        remove_atom_ids.append(atom.id)
                for atom_id in remove_atom_ids:
                    residue.detach_child(atom_id)
        save_pdb(no_sidechain_poc,output_file.replace('.pdb','_no_sidechain.pdb'))

    else:
        print(f'No pocket found for {pdbid}')
    # rm -rf {pdbid}_mergedprot_out
    shutil.rmtree(os.path.join(workdir,pdbid+'_mergedprot_out'))
    # rm -rf {pdbid}_mergedprot.pdb
    os.remove(os.path.join(workdir,pdbid+'_mergedprot.pdb'))
    return max_iou

if __name__ == "__main__":

    dataset=DTWG_DUDEDataset()
    exp_name="AF2_fpocket"
    
    for item in tqdm(dataset):
        if not os.path.exists(item['protein_dir']):
            continue
        print(item)
        home="/data/DTWG_DUD-E/raw"
        output_dir=os.path.join(item['dir'],exp_name)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        process_one_dude(home,item['name'],10,output_file=os.path.join(item['dir'],'AF2_fpocket10A.pdb'))
