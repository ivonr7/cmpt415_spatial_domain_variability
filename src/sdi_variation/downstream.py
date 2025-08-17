import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import anndata as ad
import scanpy as sc
import squidpy as sq
import logging
import numpy as np
from scipy.stats import norm
# print(f"squidpy=={sq.__version__}")
logger= logging.getLogger(__name__)
'''
    Indexes df columns with substring
'''
def method_cols(columns:list,substr:str):
    for col in columns:
        if substr in col:
            yield col


'''
    For multi regions
'''



def get_cols(columns,keys):
    for col in columns:
        for key in keys:
            if key in col:
                yield col
    
def n_square(n:int):
    return np.ceil(np.sqrt(n)).astype('int')



def merge_mask(args:list):
    if len(args) > 2:
        mask = np.logical_or(args[0],args[1])
        for arg in args[2:]:
            mask = np.logical_or(mask,arg)
    elif len(args) > 1:
        mask = np.logical_or(args[0],args[1])
    elif len(args) > 0:
        mask = args[0]
    else: return None
    return mask

def merged_c_masks(method:pd.Series,keys:list):
    masks = [cluster_mask(method,key) for key in keys]
    return merge_mask(masks)



'''
    Plots a barplot of the gene distribution
'''
def plot_gene_dist(
        count_vector:np.ndarray,
        panel_size:int,
        label:str = "",
        colour = (0,0,1),
        ax = None
    ):
    log_counts = np.log(count_vector + 1)
    ax.bar(
        list(range(count_vector.shape[0])),
        log_counts,
        alpha = 0.5,color = colour,label = label
    )
    mu = np.mean(log_counts)
    r = np.std(log_counts)
    x= np.linspace(mu - 3*r, mu + 3*r, count_vector.shape[0])
    ax.plot(list(range(count_vector.shape[0])),norm.pdf(x,mu,r),color= colour)
    ax.set_xlabel("Gene ID")
    ax.set_ylabel(f"log(gene counts + 1)")
    ax.set_xlim((0,panel_size))

def plot_gene_pdf(
    count_vector:np.ndarray,
    panel_size:int,
    label:str = "",
    colour = (0,0,1),
    ax = None,
    scale = 1
):
    log_counts = np.log(count_vector + 1)
    mu = np.mean(log_counts)
    r = np.std(log_counts)
    x= np.linspace(mu - 3*r, mu + 3*r, count_vector.shape[0])
    ax.plot(
        list(range(count_vector.shape[0])),
        norm.pdf(x,mu,r) * scale,color= colour,label=label
    )
    ax.set_xlim((0,panel_size))
    ax.set_xlabel("Gene ID")
    ax.set_ylabel(f"log(gene counts + 1)")


'''
    Index Gene Expression Matrix
'''
def cluster_mask(method:pd.Series, cluster_id:int):
    return method == cluster_id

def index_genes(gene_mat:np.ndarray,row_indexer:np.ndarray,col_indexer = None):
    if type(col_indexer) == np.ndarray:
        return gene_mat[row_indexer]
    else:
        return gene_mat[row_indexer,:][:,col_indexer]

def get_panel(sample:ad.AnnData,t:float):
    if 'moranI' in list(sample.uns.keys()):
        mI = sample.uns['moranI'] 
        goi = mI['I'] >= t 
        panel_size = goi.sum()  
    else:
        mI = None
        goi = None
        panel_size = sample.X.shape[0]
    logger.info(
        f"Gene Filtering={panel_size != sample.X.shape[0]} \
            Gene Panel={panel_size}"
        )
    return mI, goi, panel_size

'''
    Distributions Distance
'''

def bhattacharyya_distance(counts1, counts2):
    # Convert to probability distributions
    d1 = np.sum(counts1) if np.sum(counts1) > 0 else 1
    d2 = np.sum(counts2) if np.sum(counts1) > 0 else 1

    p = counts1 / d1
    q = counts2 / d1

    # Bhattacharyya coefficient
    bc = np.sum(np.sqrt(p * q))
    # Distance
    return -np.log(bc)
'''
    AnnData SVG methods
'''
def nhood(adata:ad.AnnData,clust:str = "cluster"):
    sq.gr.spatial_neighbors(adata)
    sq.gr.nhood_enrichment(
        adata,cluster_key=clust
    )
def Morans_I(adata:ad.AnnData,workers:int = 1):
    sq.gr.spatial_neighbors(adata)
    sq.gr.spatial_autocorr(
        adata,
        mode = 'moran',
        n_perms = 100,
        n_jobs = workers
    )

def sepal(adata:ad.AnnData,workers:int = 1):
    sq.gr.spatial_neighbors(adata)
    sq.gr.sepal(
        adata,
        max_neighs = 6,
        n_jobs = workers
    )
if __name__ == "__main__":

    adata = sq.datasets.visium_hne_adata()
    sq.gr.spatial_neighbors(adata)

    Morans_I(adata)
    print(adata.uns["moranI"].head())
    sq.pl.spatial_scatter(adata, color = "cluster")
    plt.show()