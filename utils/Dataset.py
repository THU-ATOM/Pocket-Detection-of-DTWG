import os
import abc
import shutil

class Dataset:
    
    def __init__(self, base_dir):
        self.base_dir = base_dir
        self.index=0
        self.items_dir=[]

    def get_base_dir(self):
        return self.base_dir
    
    def get_tmp_file_dir(self):
        return self.tmp_file_dir

    def get_items(self):
        '''
        return a list of dict, each dict contains the following keys:
        - name: the name of the item
        - dir: the directory of the item
        - pocket_dir: the directory of the pocket pdb file
        - ligand_dir: the directory of the ligand sdf file
        - protein_dir: the directory of the protein pdb file
        '''
        return self.items_dir

    def __iter__(self):
        return self
    
    def __next__(self):
        if self.index >= len(self.items_dir):
            raise StopIteration
        item = self.items_dir[self.index]
        self.index += 1
        return item

    def get_name_list(self):
        if not hasattr(self, '_name_list'):
            self._name_list = [item['name'] for item in self.items_dir]
        return self._name_list

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self.items_dir[index]
        elif isinstance(index, int):
            if index < 0:
                index += len(self.items_dir)
            if index >= len(self.items_dir) or index < 0:
                raise IndexError("Index out of range")
            return self.items_dir[index]
        else:
            raise TypeError("Index must be an integer or slice")


class PDBBindDataset(Dataset):

    def __init__(self, 
                base_dir = "/data/pdbbind_2020/combine_set",
    ):
        super().__init__(base_dir)

        self.dataset_name = "PDBBind"
        self.fasta_file = "/data/rag/PDBBind.fasta"
        self._init_item()
    
    def _init_item(self):
        import glob
        dirs = glob.glob(os.path.join(self.base_dir, '*'))
        self.items_dir = []
        for item in dirs:
            name = item.split('/')[-1]
            if name == 'index' or name == 'readme':
                continue
            new_item={
                'name' : name,
                'dir' : item,
                'pocket6A_dir' : os.path.join(item, name+"_pocket6A.pdb"),
                'pocket10A_dir' : os.path.join(item, name+"_pocket10A.pdb"),
                'ligand_dir' : os.path.join(item, name+"_ligand.mol2"),
                'protein_dir' : os.path.join(item, name+"_protein.pdb"), 
            }
            if not os.path.exists(new_item['pocket6A_dir']) or not os.path.exists(new_item['pocket10A_dir']):
                continue
            self.items_dir.append(new_item)

    def generate_targetdiff_index_pkl(self,output_dir):
        import pickle
        from tqdm import tqdm
        output=[]

        for item in tqdm(self.items_dir):
            protein_fn="/".join(item['pocket10A_dir'].split("/")[-2:])
            ligand_fn="/".join(item['ligand_dir'].split("/")[-2:])
            output.append((protein_fn,ligand_fn))
        with open(output_dir, 'wb') as f:
            pickle.dump(output, f)

class PDBBindCBDataset(Dataset):
    def __init__(self, 
                base_dir = "/data/pdbbind_cb/raw",
    ):
        super().__init__(base_dir)
        self.dataset_name = "PDBBind_cb"
        self.fasta_file = "/data/rag/PDBBind_cb.fasta"
        self._init_item()

    def _init_item(self):
        import glob
        dirs = glob.glob(os.path.join(self.base_dir, '*'))
        print(len(dirs))
        self.items_dir = []
        for item in dirs:
            name = item.split('/')[-1]
            if name == 'index' or name == 'readme':
                continue
            new_item={
                'name' : name,
                'dir' : item,
                'pocket6A_dir' : os.path.join(item, name+"_pocket6A.pdb"),
                'pocket10A_dir' : os.path.join(item, name+"_5atompoc10A.pdb"),
                'ligand_dir' : os.path.join(item, name+"_ligand.sdf"),
                'protein_dir' : os.path.join(item, name+"_protein.pdb"), 
            }
            if not os.path.exists(new_item['pocket6A_dir']) or not os.path.exists(new_item['pocket10A_dir']):
                continue
            self.items_dir.append(new_item)

    def generate_targetdiff_index_pkl(self,output_dir):
        import pickle
        from tqdm import tqdm
        output=[]
        for item in tqdm(self.items_dir):
            protein_fn="/".join(item['pocket10A_dir'].split("/")[-2:])
            ligand_fn="/".join(item['ligand_dir'].split("/")[-2:])
            output.append((protein_fn,ligand_fn))
        with open(output_dir, 'wb') as f:
            pickle.dump(output, f)

    
class BioLipDataset(Dataset):
    def __init__(self, 
                base_dir = "/data/BioLip/raw",
    ):
        super().__init__(base_dir)
        self.dataset_name = "BioLip"
        self.fasta_file = "/data/BioLip/BioLip.fasta"
        self.protein_cif_dir = "/data/BioLip/mmcif_files"
        self._init_item()

    def _init_item(self):
        import glob
        dirs = glob.glob(os.path.join(self.base_dir, '*'))
        self.items_dir = []
        for item in dirs:
            name = item.split('/')[-1]
            new_item={
                'name' : name,
                'dir' : item,
                'pocket6A_dir' : os.path.join(item, "pocket6A.pdb"),
                'pocket10A_dir' : os.path.join(item, "pocket10A.pdb"),
                'pocket15A_dir' : os.path.join(item, "pocket15A.pdb"),
                'ligand_dir' : os.path.join(item, "ligand.sdf"),
                'protein_dir' : os.path.join(self.protein_cif_dir, name.split("_")[0]+".cif"),
            }
            if not os.path.exists(new_item['ligand_dir']):
                continue
            self.items_dir.append(new_item)
    
    def generate_targetdiff_index_pkl(self,output_dir):
        import pickle
        output=[]

        for item in self.items_dir:
            protein_fn=item['name']+"/pocket10A.pdb"
            ligand_fn=item['name']+"/ligand.sdf"
            output.append((protein_fn,ligand_fn))
        with open(output_dir, 'wb') as f:
            pickle.dump(output, f)

class DUDEDataset(Dataset):
    def __init__(self, 
                base_dir = "/data/DUD-E/raw",
    ):
        super().__init__(base_dir)
        self._init_item()
        self.dataset_name = "DUD_E"
        self.fasta_file = "/data/rag/DUD_E.fasta"
    
    def _init_item(self):
        import glob
        dirs = glob.glob(os.path.join(self.base_dir, '*'))
        self.items_dir = []
        for item in dirs:
            name = item.split('/')[-1]
            self.items_dir.append({
                'name' : name,
                'dir' : item,
                'pocket6A_dir' : os.path.join(item, "pocket6A.pdb"),
                'pocket10A_dir' : os.path.join(item, "pocket10A.pdb"),
                'ligand_dir' : os.path.join(item, "crystal_ligand.mol2"),
                'protein_dir' : os.path.join(item, "receptor.pdb"), 
            })
    
    def generate_targetdiff_index_pkl(self,output_dir):
        import pickle
        output=[]

        for item in self.items_dir:
            protein_fn=item['name']+"/pocket10A.pdb"
            ligand_fn=item['name']+"/crystal_ligand.mol2"
            output.append((protein_fn,ligand_fn))
        with open(output_dir, 'wb') as f:
            pickle.dump(output, f)

class DTWG_DUDEDataset(Dataset):
    def __init__(self, 
                base_dir = "/data/DTWG_DUD-E/raw",
    ):
        super().__init__(base_dir)
        self._init_item()
        self.dataset_name = "DTWG_DUD-E"
    
    def _init_item(self):
        import glob
        dirs = glob.glob(os.path.join(self.base_dir, '*'))
        self.items_dir = []
        for item in dirs:
            name = item.split('/')[-1]
            self.items_dir.append({
                'name' : name,
                'dir' : item,
                'pocket10A_dir' : os.path.join(item, "AF2_receptor_pocket10A.pdb"),
                # 'fpocket6A_dir' : os.path.join(item, "AF2_fpocket6A.pdb"),
                'ligand_dir' : os.path.join(item, "crystal_ligand.sdf"),
                'protein_dir' : os.path.join(item, "AF2_receptor.pdb"), 
            })
    
    def generate_targetdiff_index_pkl(self,output_dir):
        import pickle
        output=[]

        for item in self.items_dir:
            if not os.path.exists(item['protein_dir']):
                continue
            protein_fn="/".join(item['pocket10A_dir'].split("/")[-2:])
            ligand_fn="/".join(item['ligand_dir'].split("/")[-2:])
            output.append((protein_fn,ligand_fn))
        with open(output_dir, 'wb') as f:
            pickle.dump(output, f)

class CrossDockedDataset(Dataset):
    def __init__(self, 
                base_dir = "/data/CrossDocked/crossdocked_pocket10",
                index_file = "/data/CrossDocked/index.pkl"
    ):
        '''
        index_file: the index file of the dataset, which is a pickle file containing a list of tuples, each tuple contains the following elements:
        - protein_fn: the filename of the protein pdb file, e.g., 'ALBU_HUMAN_25_609_halothaneSite_0/6ezq_A_rec_6ezq_c7k_lig_tt_min_0_pocket10.pdb'
        - ligand_fn: the filename of the ligand mol2 file, e.g., 'ALBU_HUMAN_25_609_halothaneSite_0/6ezq_A_rec_6ezq_c7k_lig_tt_min_0.sdf'
        '''
        super().__init__(base_dir)
        self.index_file = index_file
        self._init_item()
        self.dataset_name = "CrossDocked"
        # self.fasta_file = "/data/rag/CrossDocked.fasta"
    
    def _init_item(self):
        import pickle
        with open(self.index_file, 'rb') as f:
            index = pickle.load(f)
        for protein_fn,ligand_fn in index:
            name = ligand_fn.split('/')[-1].split('.')[0]
            self.items_dir.append({
                'name' : name,
                'pocket10A_dir' : os.path.join(self.base_dir, protein_fn),
                'ligand_dir' : os.path.join(self.base_dir, ligand_fn),
                'pocket6A_dir' : os.path.join(self.base_dir, protein_fn.replace('pocket10', 'pocket6')),
            })

    
class Dekois2Dataset(Dataset):
    def __init__(self, 
                base_dir = "/data/dekois2/",
    ):
        super().__init__(base_dir)
        self._init_item()
        self.dataset_name = "Dekois2"
        self.fasta_file = "/data/rag/Dekois2.fasta"
    
    def _init_item(self):
        import glob
        dirs = glob.glob(os.path.join(self.base_dir, '*'))
        self.items_dir = []
        for item in dirs:
            dekois_name = item.split('/')[-1]
            pdb_dir=glob.glob(os.path.join(item, '*'))[0]
            self.items_dir.append({
                'name' : dekois_name,
                'dir' : pdb_dir,
                'pocket_dir' : "/data/rag/FLAPP/dekois2_pockets/"+dekois_name+".pdb",
                'ligand_dir' : os.path.join(pdb_dir, "ligand.mol2"),
                'protein_dir' : os.path.join(pdb_dir, "protein.pdb"), 
            })


class PCBADataset(Dataset):
    def __init__(self, 
                base_dir = "/data/lit_pcba/raw",
    ):
        super().__init__(base_dir)
        self._init_item()
        self.dataset_name = "PCBA"
        # self.fasta_file = "/data/rag/DUD_E.fasta"
    
    def _init_item(self):
        import glob
        file_list=['PKM2/5x1w',
                'FEN1/5fv7',
                'TP53/2vuk',
                'MAPK1/3sa0',
                'IDH1/6b0z',
                'ESR1_ago/2b1z',
                'ALDH1/5l2o',
                'MTORC1/3fap',
                'KAT2A/5mlj',
                'GBA/2v3d',
                'PPARG/2q5s',
                'VDR/3a2j',
                'OPRK1/6b73',
                'ESR1_ant/5aau',
                'ADRB2/4ldl']
        dirs = [os.path.join(self.base_dir, item) for item in file_list]
        self.items_dir = []
        for item in dirs:
            name = item.split('/')[-2]+"/"+item.split('/')[-1]
            print(name)
            self.items_dir.append({
                'name' : name,
                'dir' : item,
                'pocket6A_dir' : item + "_pocket6A.pdb",
                'pocket10A_dir' : item + "_pocket10A.pdb",
                'pocket20A_dir' : item + "_pocket20A.pdb",
                'ligand_dir' : item + "_ligand.mol2",
                'protein_dir' : item +  "_protein.pdb",
            })
    
    def generate_targetdiff_index_pkl(self,output_dir):
        import pickle
        output=[]

        for item in self.items_dir:
            protein_fn=item['name']+"_pocket10A.pdb"
            ligand_fn=item['name']+"_ligand.mol2"
            output.append((protein_fn,ligand_fn))
        with open(output_dir, 'wb') as f:
            pickle.dump(output, f)
