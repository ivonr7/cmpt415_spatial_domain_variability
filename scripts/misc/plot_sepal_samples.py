from multiprocessing import Pool
from os import cpu_count
from pathlib import Path
import argparse
import logging
from hest import iter_hest
import squidpy as sq
from anndata import AnnData, read_h5ad
import matplotlib.pyplot as plt

from sdi_variation.downstream import sepal

def common_genes(samples:list):
    index = samples[0].index
    for sample in samples[1:]:
        index = index.intersection(sample.index)
    return index

'''
    Plot Morans I statistic for All SpatialLIBD samples
'''
def main(input_dir:str,output:str,workers:int = 12):
    sample_files = list(Path(input_dir).rglob("*.h5"))
    i_stats = []
    for sample_file in sample_files:
        sample = read_h5ad(sample_file)
        if not 'sepal_score' in list(sample.uns.keys()):
            i_stats.append(
                sepal(sample,workers)
            )
        else:
            i_stats.append(
                sample.uns['sepal_score']
            )
    index = common_genes(i_stats)
    plt.violinplot(
        [i.loc[index,'sepal_score'] for i in i_stats], vert=False
    )
    plt.yticks(
        list(range(len(i_stats) + 1)), 
        labels=[''] + [sf.stem for sf in sample_files]
    )
    plt.title("Sepal Score across Spatial Lib D")
    plt.tight_layout()
    plt.savefig(Path(output) / "sepal_score_spatialliBD.png")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process some inputs.")
    parser.add_argument('input', type=str, help='Method Folder')

    parser.add_argument('output', type=str, help='Output file or directory path')
    parser.add_argument('--workers', type=int, default=12, help='Number of worker threads (default: 12)')

    args = parser.parse_args()
    main(args.input,args.output, args.workers)