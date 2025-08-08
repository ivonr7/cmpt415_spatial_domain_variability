from pathlib import Path
from sdi_variation.mrdice import multi_region_dice
from itertools import product
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd
import numpy as np
import sys
from tqdm.auto import tqdm
import argparse
import logging



def compare_unique(sample_folder:str):
    mr_files = list(Path(sample_folder).glob("*.tif"))
    processed = set()
    for auto,gt in tqdm(product(mr_files,mr_files)):
        # Cartesian product consider (x,y) and (y,x) 
        # as distinct but dice is symetric
        # no point in calculating when auto and gt are the same
        if (gt,auto) in processed or gt == auto:
            continue
        img1 = plt.imread(auto.resolve())
        img2 = plt.imread(gt.resolve())
        mrd_data = {}
        mrd_data['scores'],mrd_data['map1'],mrd_data['map2'] = multi_region_dice(img1,img2)
        out_file = str(auto.name).split(".")[0] + "_vs_" + str(gt.name).split('.')[0] + ".csv"
        data = pd.DataFrame(
            data=mrd_data
        )
        data.to_csv(Path(sample_folder) / out_file)
        processed.add((auto,gt))

def compare_all(sample_folder:str):
    mr_files = list(Path(sample_folder).glob("*.tif"))
    score_mat = np.zeros(shape=(len(mr_files),len(mr_files)))
    for i, method1 in tqdm(enumerate(mr_files)):
        for j, method2 in enumerate(mr_files):
            img1 = plt.imread(method1.resolve())
            img2 = plt.imread(method2.resolve())
            scores,_,_ = multi_region_dice(img1,img2)
            score_mat[i,j] = np.median(scores)
        
    methods = [
        "_".join(file.stem.split('.')[0].split("_")[2:]) \
        for file in mr_files
        ]
    plt.figure(figsize=(10,10))
    sns.heatmap(
        score_mat,
        xticklabels=methods,
        yticklabels=methods,
        cmap='viridis',
        annot=True,
        fmt=".1f"
    )


    plt.title("MRDice Agreement Between Methods")
    # plt.colorbar()
    # plt.show()
    plt.tight_layout()
    plt.savefig(Path(sample_folder) / "mrdice_matrix.png")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run multi-region dice comparison on a folder of TIFF segmentations.")
    parser.add_argument("sample_folder", type=str, help="Path to folder containing .tif segmentation files.")
    parser.add_argument("-l",type=str,default=None)
    parser.add_argument('--unique',action="store_true",default=False,
                        help="to compare unique pairs or all pairs")
        
    args = parser.parse_args()
    if args.l:
        logging.basicConfig(filename=args.l,filemode='w',level=logging.INFO)
    if args.unique:
        compare_unique(args.sample_folder)
    else:
        compare_all(args.sample_folder)