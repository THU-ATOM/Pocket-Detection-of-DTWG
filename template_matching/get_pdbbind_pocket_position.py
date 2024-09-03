import sys
sys.path.append(".")
from utils.Dataset import PDBBindDataset 
import subprocess
from tqdm import tqdm
import os

def remove_gaps(seq1,seq2):
    assert len(seq1)==len(seq2)
    new_seq1=""
    new_seq2=""
    for i in range(len(seq1)):
        if seq1[i]!="-" :
            new_seq1+=seq1[i]
            new_seq2+=seq2[i]
    return new_seq1,new_seq2
    
def process_item(item):

    try:
        out_bytes = subprocess.check_output(['TMalign',item['pocket6A_dir'],item['protein_dir']])
        out_text = out_bytes.decode('utf-8').strip().split("\n")
        TMscore1=float(out_text[12].split(" ")[1])
        TMscore2=float(out_text[13].split(" ")[1])

        pocket_seq = out_text[17]
        protein_seq = out_text[19]

        if TMscore1<0.8:
            print(f"TMscore1 is less than 0.8 for {item['protein_dir']} and {item['pocket6A_dir']}: {TMscore1}")
            print(pocket_seq)
            print(protein_seq)
            return 0
        
        # print(item)
        protein_seq,pocket_seq=remove_gaps(protein_seq,pocket_seq)
        with open (os.path.join(item['dir'],"pocket6A_position.txt"),"w") as f:
            f.write(protein_seq)
            f.write("\n")
            f.write(pocket_seq)

        return 1
    except Exception as e:
        print(e)
        return 0

if __name__ == '__main__':
    pdbbind_dataset = PDBBindDataset() 
    fail_cnt=0
    for item in tqdm(pdbbind_dataset.get_items()):
        res= process_item(item)
        if res==0:
            fail_cnt+=1
    print(f"Fail count: {fail_cnt}")        
    print(f"Succeed count: {len(pdbbind_dataset.get_items())-fail_cnt}")
