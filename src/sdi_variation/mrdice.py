from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import dice
import numpy as np
from typing import Callable
import logging
from skimage.draw import random_shapes
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import pandas as pd

logger = logging.getLogger(__name__)



'''
    Implementation of multi-region dice method
    for matching segmentation regions in images
    this implements the discrete dice metric
    Author: Isaac von Riedemann
    Email: imv@sfu.ca
'''



'''
    Calculate Cost matrix based on distance function
    dist must return values where worse matches are higher
'''
def compute_cost(img:np.ndarray,gt:np.ndarray,dist:Callable = dice):
    if img.shape != gt.shape:
            logging.critical(f"auto segmentaion is {img.shape} gt is {gt.shape}")

    auto_labels = np.unique(img)
    gt_labels = np.unique(gt)
    cost_matrix = np.zeros(shape=(auto_labels.size,gt_labels.size))
    auto_mask = np.zeros(shape=img.shape)
    gt_mask = np.zeros(shape=gt.shape)
    for i in range(cost_matrix.shape[0]):
        auto_mask[:,:] = 0
        auto_mask[img == auto_labels[i]] = 1
        for j in range(cost_matrix.shape[1]):
            gt_mask[:,:] = 0
            gt_mask[gt == gt_labels[j]] = 1
            cost_matrix[i,j] = np.abs(dist(auto_mask.flatten(),gt_mask.flatten()))
    return cost_matrix


'''
    Brute Force Search
    for finding if merging any unmatched 
    regions improves score
'''
def merged_score(auto:np.ndarray,gt:np.ndarray,
                 mapping:dict, U:np.ndarray,
                 dist:Callable = dice):
    auto_merged, gt_mask = np.zeros(shape = auto.shape),np.zeros(shape=gt.shape)
    scores = {}
    n_regions = np.unique(auto.flatten()).size
    # For each Unmatched Region
    # Calculate dice to existing Ground Truth (other) seg
    # Save all the scores
    for u in U:
        # scores[u] = {}
        scores[u] = np.ones(n_regions) 
        for r_auto, r_gt in mapping.items():
            auto_merged[auto == int(u)] = 1 # unmerged region
            auto_merged[auto == int(r_auto)] = 1 # matched region

            gt_mask[gt == int(r_gt)] = 1 # gt region

            # Calculate Dice Distance = 1 - DICE
            # scores[u][r_auto] = dist(auto_merged.flatten(),gt_mask.flatten())
            scores[u][r_auto] = dist(auto_merged.flatten(),gt_mask.flatten())
            
            # Clear Masks
            auto_merged[:,:] = 0
            gt_mask[:,:] = 0
    return scores
        



def mapped_cost(img:np.ndarray,gt:np.ndarray,mapping:dict,dist:Callable = dice):
    scores = {}
    img_mask,gt_mask = np.zeros(img.shape),np.zeros(gt.shape)
    for i,g in mapping.items():
        img_mask[img == i] = 1
        gt_mask[gt == g] = 1
        scores[i] = dice(img_mask.flatten(),gt_mask.flatten())
        img_mask[:,:] = 0
        gt_mask[:,:] = 0
    all_region_scores = np.ones(
        np.unique(img.flatten()).size
    )
    all_region_scores[list(mapping.keys())] = [score for score in scores.values()]
    return all_region_scores

def merge(img,gt,mapping,mapping_score,unmerged):
        scores = merged_score(img,gt,mapping,unmerged)
        u_mapping = {}
        for u in scores.keys():
            improvement = mapping_score - scores[u]
            most_imp = np.argmax(improvement)
            u_mapping[u] = most_imp,scores[u][most_imp]
        return u_mapping      

'''
    Implements Cluster Matching Algorithmn using the 
    method given in the paper (link!!)

    computes a cost matrix using dice disimilarity measure and then 
    finds the optimal mapping between the two graphs
'''
def cluster_matching(img:np.ndarray,gt:np.ndarray):
    # Compute optimal mapping between inputs
    logger.info("calculating cost matrix")
    cost_matrix = compute_cost(img,gt)
    logger.info("finding mapping between regions")
    # Bipartite Graph Matching
    auto_map,gt_map = linear_sum_assignment(cost_matrix)
    # Using Dictionaries to store edge information
    auto_outgoing = {int(auto):int(gt) for auto,gt in zip(auto_map,gt_map)}
    gt_outgoing = {int(gt):int(auto) for auto,gt in zip(auto_map,gt_map)}
    # What is the dice score between matched regions
    dice_auto = mapped_cost(img,gt,auto_outgoing)
    dice_gt = mapped_cost(gt,img,gt_outgoing)

    if cost_matrix.max() == 0:
        logging.warning(f"No correspondance Found")
    
    # Unmatched Regions
    U1 = np.setdiff1d(np.unique(img), auto_map, assume_unique=True)
    U2 = np.setdiff1d(np.unique(gt), gt_map, assume_unique=True)

    # If unmatched aut match to gt
    if U1.size > 0: 
        u_mapping = merge(
            img,gt,auto_outgoing,
            dice_auto,U1
        )
        for u,(m_auto,score) in u_mapping.items():
            auto_outgoing[u] = auto_outgoing[m_auto]
            dice_auto[m_auto] = score

    if U2.size > 0: 
        u_mapping = merge(
            gt,img,gt_outgoing,
            dice_gt,U2
        )
        for u,(m_gt,score) in u_mapping.items():
           gt_outgoing[u] = gt_outgoing[m_gt]
           dice_gt[m_gt] = score

    # Generate Dataframe of the Mapping
    graph = pd.DataFrame()
    n_auto  = len(list(auto_outgoing.keys()))
    n_gt = len(list(gt_outgoing.keys()))

    if n_auto > n_gt:
        region_col = 's1'
        edge_col = 's2'
        regions = list(auto_outgoing.keys())
        edges = list(auto_outgoing.values())
        scores = 1- dice_auto
    elif n_gt > n_auto:
        region_col = 's2'
        edge_col = 's1'
        regions = list(gt_outgoing.keys())
        edges = list(gt_outgoing.values())
        scores = 1- dice_gt
    else:
        region_col = 's1'
        edge_col = 's2'
        regions = list(auto_outgoing.keys())
        edges = list(auto_outgoing.values())
        scores = 1- dice_auto
    graph[region_col] = regions
    graph[edge_col] = edges
    graph['dice'] = scores


    
    return graph,n_auto >= n_gt



def plot_cluster_match(img1,img2,graph:pd.DataFrame,
    titles:list = ["Method 1","Method 1 Matched","Method 2","Method 2 Matched"]):
    matched_1, matched_2= np.zeros(shape=img1.shape), np.zeros(shape=img2.shape)
    base_cmap = plt.cm.get_cmap('tab10',matched_1.shape[0])
    cmap = ListedColormap(base_cmap.colors[:matched_1.shape[0]])

    n_regions = min(
        graph['s1'].unique().size,
        graph['s2'].unique().size
    )

    for r in range(n_regions):
        over_r = graph.loc[graph['s2'] == r,'s1']
        for o_r in over_r:
            matched_1[img1 == o_r] = r
        matched_2[img2 == r] = r

    # Plot Matched masks
    plt.subplot(2,2,1)
    plt.imshow(img1,cmap=cmap)
    plt.title(titles[0])
    plt.subplot(2,2,2)
    plt.imshow(matched_1,cmap=cmap)
    plt.title(titles[1])
    plt.subplot(2,2,3)
    plt.imshow(img2,cmap=cmap)
    plt.title(titles[2])
    plt.subplot(2,2,4)
    plt.title(titles[3])
    plt.imshow(matched_2,cmap=cmap)
if __name__ == "__main__":
    clust1 = plt.imread("methods/MISC1_151676/masks/151676_deepst_3_cluster_.tif")[:,:,0].squeeze()
    clust2 = plt.imread("methods/MISC1_151676/masks/151676_deepst_4_cluster_.tif")[:,:,0].squeeze()


    mrdice_graph,order_by = cluster_matching(clust1,clust2)
    # print(auto,gt)
    print(mrdice_graph['dice'])
    plot_cluster_match(clust1,clust2,mrdice_graph)
    plt.show()