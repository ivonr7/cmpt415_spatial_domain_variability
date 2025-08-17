import matplotlib.pyplot as plt
from sdi_variation.downstream import n_square
from pathlib import Path
import logging
import argparse
logger = logging.getLogger(__name__)

def main(segmentation_folder:str,ext:str = 'tif'):
    files = sorted(list(Path(segmentation_folder).glob(f"*.{ext}")))
    fig,axes = plt.subplots(ncols = len(files),figsize = (40,5))
    for ax,file in zip(axes,files):
        im = plt.imread(file)[:,:,0].squeeze()
        ax.imshow(im,cmap = 'tab20')
        ax.set_title(file.stem)
    fig.tight_layout()
    output_folder = Path(segmentation_folder).parent / "all_segmentations"
    output_folder.mkdir(exist_ok=True)
    fig.savefig(output_folder / f"all_segs.png")


if __name__ =="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('mask_folder',type=str,required=True,help="path to folder containing SDI Segmentations")
    args = parser.parse_args()
    main(args.mask_folder)