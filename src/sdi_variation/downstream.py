import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import anndata as ad
import scanpy as sc
import squidpy as sq

print(f"squidpy=={sq.__version__}")

def nhood(adata:ad.AnnData,clust:str = "cluster"):
    sq.gr.spatial_neighbors(adata)
    sq.gr.nhood_enrichment(
        adata,cluster_key=clust
    )
def Morans_I(adata:ad.AnnData):
    sq.gr.spatial_neighbors(adata)
    sq.gr.spatial_autocorr(
        adata,
        mode = 'moran',
        n_perms = 100
    )


if __name__ == "__main__":

    adata = sq.datasets.visium_hne_adata()
    sq.gr.spatial_neighbors(adata)

    Morans_I(adata)
    print(adata.uns["moranI"].head())
    sq.pl.spatial_scatter(adata, color = "cluster")
    plt.show()