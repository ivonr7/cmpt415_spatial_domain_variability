from scipy.optimize import linear_sum_assignment
import numpy as np
'''
    Implementation of multi-region dice method
    for matching segmentation regions in images
    this implements the discrete dice metric
    Author: Isaac von Riedemann
    Email: imv@sfu.ca
'''


def discrete_dice(seg:np.ndarray, gt:np.ndarray):
    intersection = np.logical_and(seg,gt)
    union = np.logical_or(seg,gt)
    return 2 * intersection.sum() / union.sum()

def compute_cost(img:np.ndarrary,gt:np.ndarray):
    auto_labels = np.unique(img)
    gt_labels = np.unique(gt)
    cost_matrix
