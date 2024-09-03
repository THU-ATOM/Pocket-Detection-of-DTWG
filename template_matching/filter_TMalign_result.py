# this is used to filter out the TMalign result with TMscore2 > 0.5 and save in a new dir
import os
import sys
import glob
import pickle
import numpy as np
from tqdm import tqdm
if __name__ == "__main__":

    output_dir="/project/DTWG_AF2"
    output_filtered_dir="/project/DTWG_AF2_TM05"

    # output_dir="/data/DTWG_AF2"
    # output_filtered_dir="/data/DTWG_AF2_TM05"

    if not os.path.exists(output_filtered_dir):
        os.makedirs(output_filtered_dir)

    for AF2_item in tqdm(glob.glob(output_dir+"/*.pkl")):
        AF2_item_name=os.path.basename(AF2_item)
        output_filtered_item=os.path.join(output_filtered_dir,AF2_item_name)
        # read pkl
        with open(AF2_item,"rb") as f:
            res=pickle.load(f)
        new_res=[]
        for item in res:
            if item["TMscore2"]>0.5:
                new_res.append(item)
        with open(output_filtered_item,"wb") as f:
            pickle.dump(new_res,f)
        print("Number of TMalign result with TMscore2 > 0.5:",len(new_res))