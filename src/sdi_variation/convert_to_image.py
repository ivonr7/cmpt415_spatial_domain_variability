import pandas as pd
from hest import iter_hest


id = ["MISC12"]
hest_dir = "../hest_data"
sample = iter_hest(hest_dir=hest_dir,id_list=id).__next__()
clusters = pd.read_csv("/home/isaac/dev/sfu/cmpt415/spatialLIBD_annotations/151507_analysis__clustering_kmeans_2_clusters_clusters.csv")
print(clusters.head())
print(sample.adata.obs.columns)
print(sample.adata.var.head())


def splots_to_segmentation(vars:pd.DataFrame,clusters:pd.DataFrame|str):
    # Get Clusters 

    #Get 'array_col' 'array_row' to calculate adjacency matrix

    # calculate extents to get size allocate matrix

    #join clusters + indexes into 1 df

    # index matrix with indexes and set to cluster label

    pass
