from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import dice
import numpy as np
from typing import Callable
import logging
from skimage.draw import random_shapes
import matplotlib.pyplot as plt
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
            cost_matrix[i,j] = dist(auto_mask.flatten(),gt_mask.flatten())
    return cost_matrix


'''
    Brute Force Search
    for finding if merging any unmatched 
    regions improves score
'''
def merged_score(auto:np.ndarray,gt:np.ndarray,
                 M:np.ndarray, U:np.ndarray,
                 dist:Callable = dice):
    merged_scores = np.zeros(shape = (M.shape[0], U.size))
    print(merged_scores)
    auto_merged, gt_mask = np.zeros(shape = auto.shape),np.zeros(shape=gt.shape)

    for i,region in enumerate(range(M.shape[0])):
        m1,m2 = M[region,:]
        gt_mask[:,:] = 0 
        gt_mask[gt == m2] = 1
        for j,u in enumerate(U):
            auto_merged[:,:] = 0
            auto_merged[auto == m1] = 1
            auto_merged[auto == u] = 1

            merged_scores[i,j] = dice(auto_merged, gt)
    print(merged_scores)
    return merged_scores
        


def merge(auto, gt, M, col):
    U = np.setdiff1d(M[col,:],np.unique(auto))
    scores = merged_score(auto,gt,M,U)
    return U,np.argmin(scores, axis=1)



'''
    Implements Multi region dice using the 
    method given in the paper (link!!)

    computes a cost matrix using dice disimilarity measure and then 
    finds the optimal mapping between the two graphs
'''
def multi_region_dice(img:np.ndarray,gt:np.ndarray):
    # Compute optimal mapping between inputs
    logger.info("calculating cost matrix")
    cost_matrix = compute_cost(img,gt)
    logger.info("finding mapping between regions")
    auto_map,gt_map = linear_sum_assignment(cost_matrix)
    if cost_matrix.max() == 0:
        logging.warning(f"No correspondance Found")
    
    # Merging
    M = np.array(list(zip(auto_map,gt_map)))


    # If unmatched aut match to gt
    if auto_map.size != np.unique(img).size: 
        U = np.setdiff1d(M[0,:],np.unique(img))
        scores = merged_score(img,gt,M,U)
        best = np.argmin(scores, axis=1)
        auto_map = []
        for i,region in enumerate(best):
            m1,m2 = M[region,:]
            if cost_matrix[m1,m2] < scores[region,i]: auto_map.append((m1,U[i]))
    
    # if unmatched gt match to auto
    if gt_map.size != np.unique(gt).size: 
        U = np.setdiff1d(M[1,:],np.unique(gt))
        scores = merged_score(gt,img,M,U)
        print(scores)
        best = np.argmin(scores, axis=1)
        gt_map = []
        for i,region in enumerate(best):
            m1,m2 = M[region,:]
            if cost_matrix[m1,m2] < scores[region,i]: auto_map.append((m2,U[i]))
    
    return auto_map, gt_map
    




if __name__ == "__main__":
    pass
