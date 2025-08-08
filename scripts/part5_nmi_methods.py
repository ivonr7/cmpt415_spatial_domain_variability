from pathlib import Path
from sdi_variation.nmi import mutual_info
from itertools import product
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd
import numpy as np
import sys
from tqdm.auto import tqdm
import argparse
import logging



def compare_all(sample_folder:str):
    mr_files = list(Path(sample_folder).glob("*.tif"))
    score_mat = np.zeros(shape=(len(mr_files),len(mr_files)))
    for i, method1 in tqdm(enumerate(mr_files)):
        for j, method2 in enumerate(mr_files):
            img1 = plt.imread(method1.resolve())
            img2 = plt.imread(method2.resolve())
            score = mutual_info(img1,img2)
            score_mat[i,j] = score
        
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


    plt.title("NMI Score Between Methods")
    # plt.show()
    plt.tight_layout()
    plt.savefig(Path(sample_folder) / "nmi_matrix.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run multi-region dice comparison on a folder of TIFF segmentations.")
    parser.add_argument("sample_folder", type=str, help="Path to folder containing .tif segmentation files.")
    parser.add_argument("-l",type=str,default=None)

        
    args = parser.parse_args()
    if args.l:
        logging.basicConfig(filename=args.l,filemode='w',level=logging.INFO)

    compare_all(args.sample_folder)