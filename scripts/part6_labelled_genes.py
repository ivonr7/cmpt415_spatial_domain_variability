import matplotlib.pyplot as plt
import anndata as ad
import numpy as np
import logging
import pandas as pd
from pathlib import Path
from sdi_variation.downstream import method_cols,index_genes,\
    cluster_mask,plot_gene_dist,get_panel,plot_gene_pdf
logger = logging.getLogger(__name__)


'''
    plot Gene Distribution 
    across methods per cluster
'''
def per_label(
    sample_file:str, 
    t:float = 0.3,
    col_filter:str = '.csv'):

    # Get Data and ouputs
    sample = ad.read_h5ad(sample_file)
    sample_folder = Path(sample_file).parent
    plot_folder = sample_folder / "gene_dists"
    plot_folder.mkdir(exist_ok=True)
    # If Choosing A subset of genes to visualize
    _,goi,panel_size = get_panel(sample,t) 
    # Methods
    methods = list(method_cols(sample.obs.columns,substr=col_filter))
    logger.info(f"SDI Methods:\n {methods}")
    max_c = max([sample.obs[method].unique().size for method in methods])
    logger.info(f"max cluster detected {max_c}")
    cmap = plt.get_cmap('tab20')
    for i in range(max_c):
        plt.figure(figsize=(6,10))
        for j, method in enumerate(methods):
            logger.info(f"Plotting Label Distribution\n {method} Cluster {i}")
            clust = cluster_mask(sample.obs[method],i)
            clust_genes = index_genes(
                sample.X.toarray(),
                row_indexer=clust,
                col_indexer=goi
            )
            if clust_genes.size != 0:
                plot_gene_dist(
                    clust_genes.sum(axis=0),
                    panel_size,
                    label=f"{method.strip('.csv')} Cluster {i}",
                    colour=cmap(j / len(methods))
                )
        plt.title(f"Cluster {i} Label Distribution (n={panel_size})")
        plt.legend(loc='upper left', bbox_to_anchor=(0, -0.1))
        plt.tight_layout()
        plt.savefig(plot_folder / f"gene_variation_Cluster_{i}.png")
        # plt.show()

'''
    Plot Clustered Gene Distribution
    across clusters per method
'''
def per_method(sample_file:str, 
         t:float = 0.3,
         col_filter:str = '.csv',
         plot_func = plot_gene_dist):
    sample = ad.read_h5ad(sample_file)
    sample_folder = Path(sample_file).parent
    plot_folder = sample_folder / "gene_dists"
    plot_folder.mkdir(exist_ok=True)
    # If Choosing A subset of genes to visualize
    _,goi,panel_size = get_panel(sample,t) 

    # Job List
    methods = list(method_cols(sample.obs.columns,substr=col_filter))
    logger.info(f"SDI Methods:\n {methods}")
    cmap = plt.get_cmap('tab20')

    for method in methods:
        labels = sample.obs[method]
        logger.info(f"{method} found {labels.unique().size} clusters")
        logger.info("Plotting Gene Distribution Variation")
        ax=plt.figure(figsize=(10,10)).add_subplot()
        for i,label in enumerate(labels.unique()):
            logger.info(f"Plotting Cluster {label}")
            clust = cluster_mask(
                sample.obs[method],label
            )
            clust_genes = index_genes(
                sample.X.toarray(),
                row_indexer=clust,
                col_indexer=goi
            ).sum(axis=0)
            plot_func(
                clust_genes,
                panel_size,
                label=f"{method.strip('.csv')} Cluster {label}",
                colour=cmap(i / len(labels.unique())),ax=ax
            )
        plt.title(
            f"Clustered Gene Distribution(n={panel_size})\n{method.strip(".csv")}"
        )
        plt.legend(loc='upper left', bbox_to_anchor=(0, -0.1))

        plt.tight_layout()
        plt.savefig(plot_folder / f"gene_variation_{method.strip(".csv")}.png")
    



logging.basicConfig(level=logging.INFO)
per_method(
    'methods/MISC3_151674/MISC3.h5',
    col_filter="csv",
    t=0.4,
    plot_func=plot_gene_pdf
    )
