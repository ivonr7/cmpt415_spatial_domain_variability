from sdi_variation import add_annotation as aa
from hest import iter_hest
from pathlib import Path
from tqdm.auto import tqdm
from scanpy import anndata as an
from sdi_variation.downstream import Morans_I
import argparse
import logging

logger = logging.getLogger(__name__)


def main(annotation_folder:str, output:str,
         hest_dir:str = '../../hest_data', 
         svg:bool = False, workers:int =  12):
    ids = [f"MISC{i}" for i in range(1,13)]
    samples = iter_hest(id_list=ids, hest_dir=hest_dir)
    annotators = Path(annotation_folder)

    '''
        Iterate over samples and find all annotations 
        labeled with the same id ie (151507)
    '''
    for id,sample in tqdm(zip(ids,samples), desc = "adding annotations"):
        logger.info(f"Annotating {id}")
        slide = sample.meta['subseries']
        methods = [method.resolve() for method in annotators.glob(slide + '*')] # regex search for slide id
        if len(methods) == 0: 
            logger.error(f"No methods found for {id} {slide}")
            continue
        annotations = aa.join_annotators(methods) # join seperate dataframes
        sample.adata.obs = aa.join_annotation(sample.adata,annotations)

        if svg:
            logger.info("Running SVG Detection")
            Morans_I(sample.adata,workers=workers)
        sample_folder = Path(output) / "_".join([id,slide])
        if not sample_folder.exists():
            logger.info(f"Creating Sample Folder {id}")
            sample_folder.mkdir()
        sample.adata.write(sample_folder / (id+".h5"))








'''
    Chatgpt bs
'''
def parse_args():
    parser = argparse.ArgumentParser(description="Add annotations to HEST samples.")
    

    
    parser.add_argument(
        '-a',
        required=True,
        help='Path to the folder containing annotation files.'
    )
    
    parser.add_argument(
        '-o',
        required=True,
        help='Path to the output file.'
    )
    
    parser.add_argument(
        '-s',
        default='../../hest_data',
        help='Path to the HEST dataset directory (default: ../../hest_data).'
    )
    parser.add_argument(
        '-g', action='store_true',
        default=False, help="Calculate Moran's I for SVG detection"
    )
    parser.add_argument(
        '-w', type=int, default=12,
        help="Number of Workers to use for SVG"
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    logging.basicConfig(level=logging.INFO)
    if args.g:
        main(annotation_folder=args.a,
            output=args.o,
            hest_dir=args.s,
            svg=True,workers=args.w)       
    else:
        main(annotation_folder=args.a,
            output=args.o,
            hest_dir=args.s)