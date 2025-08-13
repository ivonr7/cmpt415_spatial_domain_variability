from pathlib import Path
from sdi_variation.mrdice import cluster_matching,plot_cluster_match
from itertools import product
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd
import numpy as np
import sys
from tqdm.auto import tqdm
import argparse
import logging
from itertools import repeat


def compare_unique(sample_folder:str):
    mr_files = list(Path(sample_folder).glob("*.tif"))
    processed = set()
    mrdice_folder = Path(sample_folder).parent / "mrdice"
    mrdice_folder.mkdir(exist_ok=True)
    for auto,gt in tqdm(product(mr_files,mr_files)):
        # Cartesian product consider (x,y) and (y,x) 
        # as distinct but dice is symetric
        # no point in calculating when auto and gt are the same
        if (gt,auto) in processed or gt == auto:
            continue
        img1 = plt.imread(auto.resolve())[:,:,0].squeeze()
        img2 = plt.imread(gt.resolve())[:,:,0].squeeze()

        mrd_data,_ = cluster_matching(img1,img2)

        mrd_data['method1'] = list(repeat(auto.stem.strip('_'),mrd_data.shape[0]))
        mrd_data['method2'] = list(repeat(gt.stem.strip('_'),mrd_data.shape[0]))

        out_file = str(auto.name).split(".")[0].strip('_') + "_vs_" + str(gt.name).split('.')[0].strip('_') 
        data = pd.DataFrame(
            data=mrd_data
        )
        plot_cluster_match(img1,img2,mrd_data)
        plt.savefig(mrdice_folder / (out_file + ".png"))
        plt.close()
        data.to_csv(mrdice_folder / (out_file + ".csv"),index=False)
        processed.add((auto,gt))

def compare_all(sample_folder:str,ut:bool=True):
    mr_files = list(Path(sample_folder).glob("*.tif"))
    score_mat = np.zeros(shape=(len(mr_files),len(mr_files)))
    mrdice_folder = Path(sample_folder).parent / "mrdice"
    mrdice_folder.mkdir(exist_ok=True)
    for i, method1 in tqdm(enumerate(mr_files)):
        for j, method2 in enumerate(mr_files):
            img1 = plt.imread(method1.resolve())[:,:,0].squeeze()
            img2 = plt.imread(method2.resolve())[:,:,0].squeeze()
            graph,_ = cluster_matching(img1,img2)
            score_mat[i,j] = graph['dice'].mean()
            if i != j: 
                score_mat[j,i] = graph['dice'].std()

        
    methods = [
        "_".join(file.stem.split('.')[0].split("_")[1:]) \
        for file in mr_files
        ]
    t_mask = np.triu(score_mat)
    plt.figure(figsize=(10,10))
    sns.heatmap(
        score_mat,
        xticklabels=methods,
        yticklabels=methods,
        cmap='viridis',
        annot=True,
        fmt=".1f",
    )
    plt.imshow(score_mat)


    plt.title("MRDice Agreement Between Methods")
    # plt.show()
    plt.tight_layout()
    plt.savefig(mrdice_folder / "mrdice_matrix.png")
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