import pandas as pd
import numpy as np


'''
    Convert Anndata obs table into mask to use for segmentation variation analysis
'''
def spots_to_segmentation(vars:pd.DataFrame,clusters:pd.DataFrame|str):
    # Get Clusters 
    if type(clusters) != str:
        vars = vars.join(clusters)
        clusters = "Cluster"
     
    #Get 'array_col' 'array_row' to calculate adjacency matrix
    labeled_indicies = vars.loc[:,['array_col', 'array_row',clusters]]
    # calculate extents to get size allocate matrix
    col_range = labeled_indicies['array_col'].max() 
    row_range = labeled_indicies['array_row'].max() 
    
    # # Offset in case array positions can be negative
    # labeled_indicies['array_col'] += labeled_indicies['array_col'].abs().min()
    # labeled_indicies['array_row'] += labeled_indicies['array_row'].abs().min()

    spot_image = np.zeros(shape=(row_range + 1,col_range + 1))

    # index matrix with indexes and set to cluster label
    spot_image[
        labeled_indicies['array_row'],
        labeled_indicies['array_col']
    ] = labeled_indicies[clusters]
    return spot_image



if __name__ == "__main__":
    from hest import iter_hest
    import matplotlib.pyplot as plt
    id = ["MISC12"]
    hest_dir = "../hest_data"
    sample = iter_hest(hest_dir=hest_dir,id_list=id).__next__()
    clusters = pd.read_csv(
        "/home/isaac/dev/sfu/cmpt415/spatialLIBD_annotations/151507_analysis__clustering_kmeans_2_clusters_clusters.csv",
        index_col="Barcode")
    plt.imshow(spots_to_segmentation(sample.adata.obs,clusters))
    plt.title("Spot Mask")
    plt.show()