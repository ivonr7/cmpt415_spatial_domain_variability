from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import dice
import numpy as np
from typing import Callable
import logging
from skimage.draw import random_shapes
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


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
                 M:np.ndarray, U:np.ndarray,
                 dist:Callable = dice):
    merged_scores = np.zeros(shape = (U.size,M.shape[0]))
    auto_merged, gt_mask = np.zeros(shape = auto.shape),np.zeros(shape=gt.shape)

    for i,region in enumerate(range(M.shape[0])):
        m1,m2 = M[region,:]
        gt_mask[:,:] = np.False_ 
        gt_mask[gt == m2] = np.True_
        for j,u in enumerate(U):
            auto_merged[:,:] = np.False_
            auto_merged[auto == m1] = np.True_
            auto_merged[auto == u] = np.True_

            merged_scores[j,i] = np.abs(dice(auto_merged.flatten(), gt.flatten()))
    return merged_scores
        


def merge(auto, gt, M, col):
    U = np.setdiff1d(M[col,:],np.unique(auto))
    scores = merged_score(auto,gt,M,U)
    return U,np.argmin(scores, axis=1)



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
    auto_map,gt_map = linear_sum_assignment(cost_matrix)
    if cost_matrix.max() == 0:
        logging.warning(f"No correspondance Found")
    
    # Merging
    M1,M2 = [],[]
    U1 = np.setdiff1d(np.unique(img), auto_map, assume_unique=True)
    U2 = np.setdiff1d(np.unique(gt), gt_map, assume_unique=True)
    M = np.array(list(zip(auto_map,gt_map)))

    auto_map,gt_map = np.expand_dims(auto_map,axis = 1).tolist(),np.expand_dims(gt_map,axis = 1).tolist()
    # If unmatched aut match to gt
    if U1.size > 0: 
        scores = merged_score(img,gt,M,U1)
        for i,u in enumerate(U1):
            for region, score in enumerate(scores[i,:]):
                m1,m2 = M[region,:]
                scores[i,region] = cost_matrix[m1,m2] - score
            # If there was an improvement ie cost_mat > merged_score
            # then add
            if scores[i,:].max() > 0:
                logger.info(f"Merging U1 score improvement! {scores[i,:].max()}")
                auto_map[np.argmax(scores[i,:])] = (*auto_map[region],u)
                

    if U2.size > 0: 
        scores = merged_score(gt,img,M,U2)
        for i,u in enumerate(U2):
            for region, score in enumerate(scores[i,:]):
                m1,m2 = M[region,:]
                scores[i,region] = cost_matrix[m1,m2] - score
            # If there was an improvement ie cost_mat > merged_score
            # then add
            if scores[i,:].max() > 0:
                logger.info(f"Merging U2 score improvement! {scores[i,:].max()}")
                gt_map[np.argmax(scores[i,:])] = [*gt_map[region],int(u)]
    return auto_map, gt_map

'''
    Returns Dice Score between Clusters
'''
def multi_region_dice(img1,img2):
    matched1, matched2 = cluster_matching(img1,img2)
    scores = []
    for region in zip(matched1,matched2):
        region_mask1 = np.zeros(img1.shape)
        region_mask2 = np.zeros(img2.shape)
        for clust in region[0]:
            region_mask1[img1 == int(clust)] = 1
        for clust in region[1]:
            region_mask2[img2 == int(clust)] = 1
        # Divide by zero protection
        if region_mask1.max() == 0 and region_mask2.max() == 0:
            scores.append(0)
        else:
            scores.append(
                1 - dice(region_mask1.flatten(),region_mask2.flatten())
            )
    return np.array(scores), matched1, matched2

def plot_cluster_match(img1,img2,map1,map2):
    matched_1, matched_2= np.zeros(shape=img1.shape), np.zeros(shape=img2.shape)
    base_cmap = plt.cm.get_cmap('tab10',matched_1.shape[0])
    cmap = ListedColormap(base_cmap.colors[:matched_1.shape[0]])
    for i,(region1,region2) in enumerate(zip(map1,map2)):
        for clust in region1:
            matched_1[img1 == clust] = i
        for clust in region2: 
            matched_2[img2 == clust] = i 
    plt.subplot(2,2,1)
    plt.imshow(img1,cmap=cmap)
    plt.title("Method 1")
    plt.subplot(2,2,2)
    plt.imshow(matched_1,cmap=cmap)
    plt.title("Method 1 Matched")
    plt.subplot(2,2,3)
    plt.imshow(img2,cmap=cmap)
    plt.title("Method 2")
    plt.subplot(2,2,4)
    plt.title("Method 2 Matched")
    plt.imshow(matched_2,cmap=cmap)
    plt.show()
if __name__ == "__main__":
    clust1 = plt.imread("methods/MISC1_151676/masks/151676_deepst_3_cluster_.tif")[:,:,0].squeeze()
    clust2 = plt.imread("methods/MISC1_151676/masks/151676_deepst_6_cluster_.tif")[:,:,0].squeeze()


    auto, gt = cluster_matching(clust1,clust2)
    print(auto,gt)
    plot_cluster_match(clust1,clust2,auto,gt)