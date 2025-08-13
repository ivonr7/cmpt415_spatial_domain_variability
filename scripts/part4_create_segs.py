from sdi_variation.convert_to_image import spots_to_segmentation
from PIL import Image
import re
from pathlib import Path
from anndata.io import read_h5ad
import argparse
def csv_columns(columns:list):
    for column in columns:
        match = re.match(r".*\.csv",column)
        if match:
            yield column

def main(sample:str,output:str, ext:str = ".tif"):
    adata = read_h5ad(sample)
    annotators = list(csv_columns(adata.obs.columns))
    id = Path(sample).stem
    mask_folder = Path(output) / "masks"
    if not mask_folder.exists(): 
        mask_folder.mkdir()
    for annotater in annotators:
        img,cols,rows = spots_to_segmentation(adata.obs, annotater)
        print(cols,rows)
        Image.fromarray(img).save(
            mask_folder / "_".join((annotater.strip(".csv"),ext))
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a sample file and output with specified extension.")
    
    parser.add_argument(
        "sample",
        type=str,
        help="Path to the input sample file"
    )

    parser.add_argument(
        "output",
        type=str,
        help="Path to the output file or folder"
    )

    parser.add_argument(
        "--ext",
        type=str,
        default=".tif",
        help="File extension to use (default: .tif)"
    )

    args = parser.parse_args()
    main(args.sample, args.output, args.ext)