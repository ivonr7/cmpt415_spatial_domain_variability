import matplotlib.pyplot as plt
import anndata as ad
import numpy as np
import logging
import pandas as pd
from pathlib import Path
from tqdm.auto import tqdm
from sdi_variation.downstream import method_cols,index_genes,\
    cluster_mask,plot_gene_dist,get_panel
import ast
logger = logging.getLogger(__name__)

def get_cols(columns,keys):
    for col in columns:
        for key in keys:
            if key in col:
                yield col
    
def n_square(n:int):
    return np.ceil(np.sqrt(n)).astype('int')

def str2list(col):
    items = []
    for row in col:
        try:
            logging.info(row) 
            items.append(ast.literal_eval(row))
        except Exception as e:
            logger.error(f'{e}')
    return items

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
def plot_match(
        sample:ad.AnnData,region:pd.DataFrame,
        goi:np.ndarray,panel_size:int):
    # logger.info(region.head())
    # MRDICE Mapped Methods
    region[['map1','map2']] = region[['map1','map2']].apply(str2list,axis=0)
    methods = list(
        get_cols(sample.obs.columns,region.loc[0,['method1','method2']])
    )
    logger.info(region.loc[0,['method1','method2']])
    logger.info(methods)
    m1 = methods[0]
    m2 = methods[1]
    # Plot paremeters
    s_len = n_square(region.shape[0])
    fig, axes = plt.subplots(s_len,s_len, figsize = (5*s_len,5 * s_len))
    for i,row in tqdm(enumerate(region.itertuples(index = False))):
        
        c1 = row.map1
        c2 = row.map2
        score = row.scores
        # Generate Mask of Spots Assigned to Each Cluster
        cl_1 = merged_c_masks(sample.obs[m1],c1)
        cl_2 = merged_c_masks(sample.obs[m2],c2)
        
        # Index Gene Expression Matrix to get cluster gene dist
        gene_mat = sample.X.toarray()
        g_dist_1 = index_genes(gene_mat,cl_1,goi).sum(axis = 0)
        g_dist_2 = index_genes(gene_mat, cl_2,goi).sum(axis = 0)

        # Plot
        x,y = i % s_len, i // s_len
        ax = axes[y,x]
        plot_gene_dist(
            g_dist_1,
            panel_size=panel_size,
            colour=(0,0,1),
            label=f"Cluster {c1}",
            ax=ax
        )
        plot_gene_dist(
            g_dist_2,
            panel_size=panel_size,
            colour=(1,0,0),
            label=f"Cluster {c2}",
            ax=ax
        )
        ax.set_title(f"Cluster {c1} Cluster {c2} Gene Distribution\n (Optimal with MRD={score:.3f})")
        ax.legend()
    plt.suptitle(f"{m1} vs {m2}")
    fig.tight_layout()
    return f"{m1} vs {m2}"



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
        fname = plot_match(sample,region,
                   goi,panel_size)
        plt.savefig(out_folder / (fname + '.png'))
        plt.close()

logging.basicConfig(level=logging.INFO)
folder = Path("methods/MISC1_151676")
plot_matches(folder / "MISC1.h5",folder / "mrdice")