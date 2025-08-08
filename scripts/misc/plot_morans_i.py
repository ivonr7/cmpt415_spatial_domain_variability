from multiprocessing import Pool
from os import cpu_count
from pathlib import Path
import argparse
import logging
from hest import iter_hest
import squidpy as sq
from anndata import AnnData
import matplotlib.pyplot as plt

def moransI(id:str,adata:AnnData,workers:int):
    sq.gr.spatial_neighbors(adata)
    sq.gr.spatial_autocorr(
        adata,
        mode="moran",
        n_perms=100,
        n_jobs=workers,
    )
    return id, adata.uns['moranI'].fillna(-1)

'''
    Plot Morans I statistic for All SpatialLIBD samples
'''
def main(output:str,hest_dir:str = '../hest_data',workers:int = 12):
    ids = [f"MISC{i}" for i in range(1,13)]
    samples = iter_hest(hest_dir=hest_dir,id_list=ids)
    i_stats = []
    for sample in samples:
        i_stats.append(
            moransI(sample.meta['id'],sample.adata,workers)
        )
    plt.violinplot(
        [i[1]['I'] for i in i_stats], vert=False
    )
    plt.yticks(
        list(range(len(i_stats) + 1)), 
        labels=[''] + [i[0] for i in i_stats]
    )
    plt.title("Spatial Autocorrelation across Spatial Lib D")
    plt.tight_layout()
    plt.savefig(Path(output) / "morans_i_spatialliBD.png")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process some inputs.")
    parser.add_argument('output', type=str, help='Output file or directory path')
    parser.add_argument('--hest_dir', type=str, default='../hest_data', help='Directory containing hest data (default: ../hest_data)')
    parser.add_argument('--workers', type=int, default=12, help='Number of worker threads (default: 12)')

    args = parser.parse_args()
    main(args.output, args.hest_dir, args.workers)