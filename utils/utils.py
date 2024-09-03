import os
from Bio.PDB import PDBParser,Chain,Model,Structure
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import is_aa
from Bio.PDB.Residue import DisorderedResidue,Residue
from Bio.PDB.Atom import DisorderedAtom
import warnings
from Bio.PDB.StructureBuilder import PDBConstructionWarning
import numpy as np
from rdkit.Chem import Descriptors
import copy
import numpy as np
#from biopandas.mol2 import PandasMol2
from rdkit import Chem
import rdkit.Chem.AllChem as AllChem
from tqdm import tqdm
from io import StringIO
from rdkit.Chem.rdmolfiles import MolToPDBBlock
from Bio.PDB.NeighborSearch import NeighborSearch
import multiprocessing as mp
from freesasa import calcBioPDB
import freesasa
from Bio.PDB.Polypeptide import protein_letters_3to1
from copy import deepcopy
import lmdb
import pickle

warnings.filterwarnings(
    action='ignore',
    category=PDBConstructionWarning)

def res_is_connected(residue1,residue2):
    ca1 = residue1['CA']
    ca2 = residue2['CA']
    distance = ca1-ca2
    if abs(distance-3.8)<=0.2:
        return 1
    else:
        return 0 

def merge_chains(pdbstruct,af2=False,tgt=None):
    break_point = []
    #p = PDBParser()
    model = pdbstruct[0]
    tmp_chain = Chain.Chain('A') 
    rid = 0
    last_res = None
    for chain in model:
        break_point.append(rid)
        for res in chain:
            try:
                res.detach_parent()
                if not(is_aa(res,standard=True)):
                    continue 
                if af2 and res['CA'].bfactor<50:
                    continue
                if res.is_disordered():
                    if isinstance(res,DisorderedResidue):
                        res = res.selected_child
                        res.id = (res.id[0],rid,res.id[2])      
                    else:
                        new_res = Residue(res.id,res.resname,res.segid)
                        for atom in res:
                            if isinstance(atom,DisorderedAtom):  
                                atom.selected_child.disordered_flag = 0
                                new_res.add(atom.selected_child.copy())
                            else:
                                new_res.add(atom)
                        res = new_res
                        res.id = (res.id[0],rid,res.id[2]) 
                else:
                    res.id = (res.id[0],rid,res.id[2])
                
                if last_res is not None and not res_is_connected(res,last_res):
                    break_point.append(rid)
                last_res = copy.deepcopy(res)
                tmp_chain.add(res.copy())
                rid += 1 
            except:
                pass
    tmp_structure = Structure.Structure(pdbstruct.id)
    tmp_model = Model.Model(0)
    tmp_structure.add(tmp_model) 
    tmp_model.add(tmp_chain)
    if tgt is not None:
        io = PDBIO()
        io.set_structure(tmp_structure)
        io.save(os.path.join(tgt,pdbstruct.id+'.pdb'))    
    return tmp_structure,break_point


def read_mol2_ligand(path,return_pdbstring = False):
    '''
    mol = Chem.MolFromMol2File(path, removeHs=True, sanitize=False)
    if mol is None:
        print("cannot read mol2", path)
    coords = mol.GetConformer().GetPositions()
    atom_types = [a.GetSymbol() for a in mol.GetAtoms()]
    return {'coord': np.array(coords), 'atom_type': atom_types, 'mol': mol, 'smi': Chem.MolToSmiles(mol)}
    '''
    # try:
    #     mol2_df = PandasMol2().read_mol2(path ,columns={0:('atom_id', int), 1:('atom_name', str), 2:('x', float), 3:('y', float), 4:('z', float), 5:('atom_type', str), 6:('subst_id', int), 7:('residue_name', str), 8:('useless1', float), 9:('useless2', str)})

    #     coords = mol2_df.df[['x', 'y', 'z']]
    # except:
    mol = Chem.MolFromMol2File(path,sanitize=False)
    try:
        mw = Descriptors.MolWt(mol)
        mol_noH = Chem.RemoveHs(mol)
    except:
        mw = 0
        mol_noH = mol
    coords = mol_noH.GetConformer().GetPositions()
    if return_pdbstring:
        pdb_string = MolToPDBBlock(mol,flavor=2)
        return np.array(coords),mw,pdb_string
    else:
        return np.array(coords),mw

def get_binding_pockets(original_pdb,lig_coord,thres=6,output=None):
    chain = original_pdb[0]['A'] #only deal with A chain
    tmp_chain = Chain.Chain('A') 
    resid = set()
    for res in chain:
        res_coord = np.array([i.get_coord() for i in res.get_atoms()])
        dist = np.linalg.norm(res_coord[:,None,:]-lig_coord[None,:,:],axis=-1).min()
        if dist<=thres:
            tmp_chain.add(res.copy())
            resid.add(res.id[1])
    tmp_structure = Structure.Structure(original_pdb.id)
    tmp_model = Model.Model(0)
    tmp_structure.add(tmp_model) 
    tmp_model.add(tmp_chain)
    if output is not None:
        io = PDBIO()
        io.set_structure(tmp_structure)
        io.save(os.path.join(output,original_pdb.id+'.pdb'))
    atom_num = 0
    for a in tmp_chain.get_atoms():
        if a.element != 'H':
            atom_num += 1
        else:
            #delete hydrogen atoms
            a.detach_parent()
    return tmp_structure,resid,atom_num
    
def pqr_parser(filename,return_score=False):
    with open(filename,'r') as f:
        data = f.readlines()
    coord = []
    for l in data:
        if "Pocket Score" in l:
            score = float(l.split()[-1])
        elif "Real volume" in l:
            volume = float(l.split()[-1])
        elif l[:4] == 'ATOM':
            coord.append(
                [float(l[30:38]),float(l[38:46]),float(l[46:54])]
            )
    coord = np.array(coord)
    if return_score:
        return coord,score,volume
    else:
        return coord

def IoU(setA,setB):
    iou = len(setA & setB)/len(setA | setB)
    return iou


def getfpockets(pdb_struct,pockets_dir,thres,output=None,return_score=False):
    # original_pdb = os.path.join(home,pdbfile)
    # pockets_dir = os.path.join(home,pdbfile).replace('.pdb','_out/pockets')
    # p = PDBParser()
    # pdb_struct = p.get_structure('0', original_pdb)
    pocket_list = os.listdir(pockets_dir)
    pocpdb_list,res_list,anum_list,score_list,volume_list = [],[],[],[],[]
    i = 0
    pdb_struct_id = deepcopy(pdb_struct.id)
    for f in pocket_list:
        if f[-4:] == '.pqr':
            pdb_struct.id = pdb_struct_id.replace('_ligpocket',f'_fpocket{i}')
            i += 1
            lig_coord,score,volume = pqr_parser(os.path.join(pockets_dir,f),return_score=True)
            pocket,resid,atom_num = get_binding_pockets(pdb_struct,lig_coord,thres,output)
            pocpdb_list.append(pocket)
            res_list.append(resid)
            anum_list.append(atom_num)
            score_list.append(score)
            volume_list.append(volume)
    if return_score:
        return pocpdb_list,res_list,np.array(anum_list),np.array(score_list),np.array(volume_list)
    else:
        return pocpdb_list,res_list,anum_list

def save_pdb(pdb_struct,output):
    io = PDBIO()
    io.set_structure(pdb_struct)
    io.save(output)

def pdb2dict(biopy_structure,pocket_name):
    recpt = list(biopy_structure.get_atoms())
    pocket_atom_type = [x.element for x in recpt]
    pocket_coord = [x.coord for x in recpt]
    return {
        'pocket': pocket_name,
        'pocket_atoms': pocket_atom_type,
        'pocket_coordinates': pocket_coord
    }

def write_lmdb(data, lmdb_path,num=0):
    #resume
    env = lmdb.open(lmdb_path, subdir=False, readonly=False, lock=False, readahead=False, meminit=False, map_size=1099511627776)
    with env.begin(write=True) as txn:
        for d in tqdm(data):
            txn.put(str(num).encode('ascii'), pickle.dumps(d))
            num += 1
    return num

if __name__ == '__main__':
    p = PDBParser()
    model = p.get_structure('0', '/home/jiayinjun/Chem_Awared_PockEnc/4b32.pdb')
    model,breakpoint = merge_chains(model,tgt='/home/jiayinjun/Chem_Awared_PockEnc/')

