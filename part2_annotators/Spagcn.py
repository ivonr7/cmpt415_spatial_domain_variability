import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import SpaGCN as spg
from hest import iter_hest
from PIL import Image
import openslide
from pathlib import Path
import logging
import argparse 
logger = logging.getLogger(__name__)
def main(hest_folder:str,save_folder:str,n_clusters:int = 7):
    ids = [f"MISC{i}" for i in range(1,13)]
    for i,sample in enumerate(iter_hest(hest_dir=hest_folder,id_list= ids)):

        # Pull out important Variables
        slide = sample.meta["subseries"]
        wsi_file = Path(hest_folder) / "wsis" / (ids[i] + ".tif")
        wsi = np.array(Image.open(wsi_file))
        adata = sample.adata
        logging.info(f"{slide} slide id")
        logging.info(f"Hest ID: {ids[i]}")
        logging.info(f"wsi file {wsi_file}")


        '''
            Code to run SpaGCN from their tutorial_ez_mode
        '''
        #Set coordinates
        adata.obs["x_array"]=adata.obs["array_col"]
        adata.obs["y_array"]=adata.obs["array_row"]
        adata.obs["x_pixel"]=adata.obs["pxl_col_in_fullres"]
        adata.obs["y_pixel"]=adata.obs["pxl_row_in_fullres"]
        x_array=adata.obs["x_array"].tolist()
        y_array=adata.obs["y_array"].tolist()
        x_pixel=adata.obs["x_pixel"].astype(np.int64).tolist()
        y_pixel=adata.obs["y_pixel"].astype(np.int64).tolist()
        adata.X = adata.X.toarray()
        print(x_pixel)
        #Run SpaGCN
        #TODOset as parameters to do full sweeps
        adata.obs["pred"]= spg.detect_spatial_domains_ez_mode(
            adata, wsi, x_array, y_array, x_pixel, y_pixel,
            n_clusters=n_clusters, histology=True, s=1, b=49, p=0.5, 
            r_seed=100, t_seed=100, n_seed=100
        )
        adata.obs["pred"]=adata.obs["pred"].astype('category')
        #Refine domains (optional)
        #shape="hexagon" for Visium data, "square" for ST data.
        adata.obs["refined_pred"]=spg.spatial_domains_refinement_ez_mode(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), x_array=x_array, y_array=y_array, shape="hexagon")
        adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')

        adata.obs["refined_pred"].reset_index().to_csv(
            Path(save_folder) / (slide + "_" + ids[i] + f"spagcn_{n_clusters}_clusters.csv"),
            header = ["Barcode","Cluster"],
            index = False
        )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run SpaGCN on HEST data")
    parser.add_argument(
        "hest_folder",
        type=str,
        help="Path to the folder containing HEST data (and 'wsis' subfolder)"
    )
    parser.add_argument(
        "save_folder",
        type=str,
        help="Folder to save SpaGCN results"
    )
    parser.add_argument(
        '-k', default=7,type=int,help="number of clusters"
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    main(args.hest_folder, args.save_folder,n_clusters=args.k)