from sdi_variation import add_annotation as aa
from hest import iter_hest
from pathlib import Path
from tqdm.auto import tqdm
from scanpy import anndata as an
import argparse
def main(annotation_folder:str, output:str,hest_dir:str = '../../hest_data'):
    ids = [f"MISC{i}" for i in range(1,13)]
    samples = iter_hest(id_list=ids, hest_dir=hest_dir)
    annotators = Path(annotation_folder)

    '''
        Iterate over samples and find all annotations 
        labeled with the same id ie (151507)
    '''
    for id,sample in tqdm(zip(ids,samples), desc = "adding annotations"):
        
        slide = sample.meta['subseries']
        methods = [method.resolve() for method in annotators.glob(slide + '*')] # regex search for slide id
        annotations = aa.join_annotators(methods) # join seperate dataframes
        sample.adata.obs = aa.join_annotation(sample.adata,annotations)
        sample.adata.write(Path(output) / (id+".h5"))








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
    
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(annotation_folder=args.a,
         output=args.o,
         hest_dir=args.s)