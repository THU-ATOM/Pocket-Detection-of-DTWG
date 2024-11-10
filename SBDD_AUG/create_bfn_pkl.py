# create a pkl file for bfn



key_word_list=[
    "A0A5K1K8R3"
    "A0A5K1K8U3",
    "A0A5K1K930",
    "C0H4G1",
    "C0H4N6",
    "O96139",
    "Q7KWI6",
    "Q8I2A8",
    "Q8I2Q0",
    "Q8I3V7",
    "Q8I4R4",
    "Q8I5A2",
    "Q8I5A6",
    "Q8I298",
    "Q8I470",
    "Q8I578",
    "Q8IDR1",
    "Q8IFN7",
    "Q8IIV4",
    "Q8IJI2",
    "Q8IJV0",
    "Q8IJV2",
    "Q8IL23",
    "Q8ILM9",
]



import os
import glob

fpocket_path="/data/Plasmodium_screening/AF2_fpocket"

pkl_output_path="/data/Plasmodium_screening/genpack_result/bfn_output/fpocket.pkl"
ref_ch4_file="/project/ch4.sdf"
pockets=glob.glob(os.path.join(fpocket_path,"*.pdb"))

# filter pockets
pockets=[p for p in pockets if any([key in p for key in key_word_list])]
print(len(pockets))



result=[]

for pocket in pockets:
    pocket_name=os.path.basename(pocket).replace(".pdb","")
    ref_ligand_path=ref_ch4_file
    pocket_path=pocket

    result.append({
        "pocket_name":pocket_name,
        "pocket_path":pocket_path,
        "ref_ligand_path":ref_ligand_path
    })

import json
# print(json.dumps(result,indent=4))
import pickle
with open(pkl_output_path,"wb") as f:
    pickle.dump(result,f)
print("Finish")
