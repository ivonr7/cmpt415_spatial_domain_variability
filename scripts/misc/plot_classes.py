import scanpy as sc
import squidpy as sq
import pandas  as pd
from pathlib import Path
import argparse
def plot_annotations(h5_data:str,output:Path):
    adata = sc.read_h5ad(h5_data)
    cols = pd.Series(adata.obs.columns)
    for annotator in cols[cols.str.contains("cluster")]:

        fig = sq.pl.spatial_scatter(
                adata,
                color = annotator,return_fig=True
            )
        fig.savefig(Path(output) / (annotator + ".png"))


def parse_args():
    parser = argparse.ArgumentParser(description="Plot annotation clusters from h5ad data.")
    
    parser.add_argument(
        '--h5-data',
        type=str,
        required=True,
        help='Path to the input .h5ad file containing annotated spatial data.'
    )

    parser.add_argument(
        '--output',
        type=Path,
        required=True,
        help='Directory where the output plots will be saved.'
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    plot_annotations(h5_data=args.h5_data, output=args.output)