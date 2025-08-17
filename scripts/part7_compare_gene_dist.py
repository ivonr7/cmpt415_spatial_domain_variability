import matplotlib.pyplot as plt
import anndata as ad
import numpy as np
import logging
import pandas as pd
from pathlib import Path
from tqdm.auto import tqdm
from sdi_variation.downstream import method_cols,index_genes,\
    cluster_mask,plot_gene_dist,get_panel,bhattacharyya_distance,get_cols,n_square,merged_c_masks
import ast
logger = logging.getLogger(__name__)



def plot_all(
        sample:ad.AnnData,region:pd.DataFrame,
        goi:np.ndarray,panel_size:int,save_folder:Path
):
    # logger.info(region.head())
    # MRDICE Mapped Methods
    methods = list(
        get_cols(sample.obs.columns,region.loc[0,['method1','method2']])
    )
    logger.info(region.loc[0,['method1','method2']])
    logger.info(methods)
    m1 = methods[0]
    m2 = methods[1]
    # Plot paremeters
    under_col = region[['s1','s2']].max(axis = 0).idxmin()
    over_col = region[['s1','s2']].max(axis = 0).idxmax()
    x,y = region[['s1','s2']].max(axis = 0)

def plot_match(
        sample:ad.AnnData,region:pd.DataFrame,
        goi:np.ndarray,panel_size:int,save_folder:Path
    ):
    # logger.info(region.head())
    # MRDICE Mapped Methods
    methods = list(
        get_cols(sample.obs.columns,region.loc[0,['method1','method2']])
    )
    logger.info(region.loc[0,['method1','method2']])
    logger.info(methods)
    m1 = methods[0]
    m2 = methods[1]
    # Plot paremeters
    under_col = region[['s1','s2']].max(axis = 0).idxmin()
    over_col = region[['s1','s2']].max(axis = 0).idxmax()
    region['bd'] = np.float64(0)
    for r in tqdm(region[under_col].unique()):
        merged = region.loc[region[under_col] == r, over_col]
        if np.isnan(r) or np.isnan(merged).any():
            logger.error("Nan values encountered")
            continue
        ax = plt.figure(figsize=(10,10)).add_subplot()


        sr1_genes = merged_c_masks(
            sample.obs[m1],
            region.loc[merged,over_col].values.tolist())
        sr2_genes = merged_c_masks(
            sample.obs[m2],
            [region.loc[r,under_col]]
        )
        
        gene_counts = sample.X.toarray()
        sr1_counts = index_genes(gene_counts,sr1_genes,goi).sum(axis=0)
        sr2_counts = index_genes(gene_counts,sr2_genes,goi).sum(axis = 0)
        d = bhattacharyya_distance(sr1_counts,sr2_counts)
        region.loc[region[under_col] == r,'bd'] = d

        plot_gene_dist(
            sr1_counts,
            panel_size=panel_size,
            colour=(0,0,1),
            label=f"Cluster {merged.values.tolist()}",
            ax = ax
        )
        plot_gene_dist(
            sr2_counts,
            panel_size=panel_size,
            colour=(1,0,0),
            label=f"Cluster {r}",
            ax = ax
        )
        ax.legend()
        ax.set_title(f"Region {merged.values.tolist()} Region {r} Gene Distribution (Optimal BD={d:.4f})")
        plt.tight_layout()
        output_folder = save_folder / f"{m1}_vs_{m2}"
        output_folder.mkdir(exist_ok=True)
        plt.savefig(output_folder / f"Optimal_region{r}.png")
        region.to_csv(output_folder / f"{m1}_vs_{m2}_scored.csv")




'''
    Plot Matches for all methods for sample
'''    
def plot_matches(sample_file:str,mrd_folder:str,t:float = 0.4):
    sample = ad.read_h5ad(sample_file)
    _, goi,panel_size = get_panel(sample,t = t)
    # Getting Files together
    region_maps = Path(mrd_folder).glob('*.csv')
    sample_folder = Path(sample_file).parent
    out_folder = sample_folder / "pair_genes"
    out_folder.mkdir(exist_ok=True)
    for region_file in region_maps:
        region = pd.read_csv(
            region_file
        )
        plot_match(sample,region,
                   goi,panel_size,save_folder=out_folder)

logging.basicConfig(level=logging.INFO)
folder = Path("clgraph/MISC3_151674")
plot_matches(folder / "MISC3.h5",folder / "mrdice")
