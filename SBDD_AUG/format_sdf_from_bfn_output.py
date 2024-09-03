import os
import glob
import rdkit 
from rdkit import Chem
from tqdm import tqdm
import random
import shutil
from Drug_The_Whole_Genome.utils.Dataset import DTWG_DUDEDataset

def format_output():
    # exps=glob.glob('/data/pocket2mol_data/sample_output/*')
    input_exps= '/data/bfn_data/tanhaichuan_bfn_sbdd/pdbbind_cb_09/default/test_outputs_v19/good_cases'
    
    dataset=DTWG_DUDEDataset()
    exp_name="AF2_fpocket_BFN"

    sdf_list={}
    samples=glob.glob(input_exps+'/*.sdf')

    for sample in samples:
        case_name="_".join(sample.split('/')[-1].split('_')[:-1])
        if case_name not in sdf_list:
            sdf_list[case_name]=[]
        sdf_list[case_name].append(sample)

    for case_name, cases in tqdm(sdf_list.items()):
        output_dir=os.path.join(dataset.base_dir, case_name, exp_name)
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        for i , case in enumerate(cases):
            shutil.copy(case, os.path.join(output_dir, f'{i}.sdf'))

if __name__ == '__main__':
    format_output()