import os
import glob
import rdkit 
from rdkit import Chem
from tqdm import tqdm
import random
import shutil

def format_pocket2mol_output():
    # exps=glob.glob('/data/pocket2mol_data/sample_output/*')
    input_exps= '/data/bfn_data/tanhaichuan_bfn_sbdd/CrossDocked0.6/default/test_outputs_v4/good_cases'
    output_exps= 'sa'

    if os.path.exists(output_exps):
        shutil.rmtree(output_exps)
    os.makedirs(output_exps,exist_ok=True)
    sdf_list={}
    samples=glob.glob(input_exps+'/*')

    for sample in samples:
        case_name="_".join(sample.split('/')[-1].split('_')[:-1])
        if case_name not in sdf_list:
            sdf_list[case_name]=[]
        sdf_list[case_name].append(sample)

    for case_name, cases in tqdm(sdf_list.items()):
        output_file=os.path.join(output_exps,case_name+'.sdf')
        print(output_file)
        writer=Chem.SDWriter(output_file)
        for case in cases:
            mol=Chem.MolFromMolFile(case)
            if mol is not None:
                writer.write(mol)
            else:
                print('Error: %s'%case)
        writer.close()

if __name__ == '__main__':
    format_pocket2mol_output()