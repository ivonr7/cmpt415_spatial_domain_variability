from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import dice
import numpy as np
from typing import Callable
import logging

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
    logging.critical(f"auto segmentaion is {img.shape} gt is {gt.shape}")
    assert img.shape == gt.shape
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
    Implements Multi region dice using the 
    method given in the paper (link!!)

    computes a cost matrix using dice disimilarity measure and then 
    finds the optimal mapping between the two graphs
'''
def multi_region_dice(img:np.ndarray,gt:np.ndarray):
    # Compute optimal mapping between inputs
    logging.info("calculating cost matrix")
    cost_matrix = compute_cost(img,gt)
    logging.info("finding mapping between regions")
    auto_map,gt_map = linear_sum_assignment(cost_matrix)
    if auto_map.size() < np.unique(img).size():
        logging.warning(f"Method Under Segmented with \
                        {np.unique(img).size() - auto_map.size()} regions left")
    elif auto_map.size() > np.unique(img).size():
        logging.warning(f"Method Over Segmented with \
                        {auto_map.size() - np.unique(img).size()} regions extra")

    # Merging
    # For visualization
    optimal_map = np.zeros(shape=img.shape)
    for idx,(i,j) in enumerate(zip(auto_map,gt_map)):
        optimal_map[img == i] = idx
        optimal_map[gt == j] = idx
    plt.figure()
    plt.imshow(optimal_map)
    plt.title("OR between matched regions")
    

if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    from skimage.draw import disk

    # Create an empty 2D array
    shape = (100, 100)
    ground_truth = np.zeros(shape, dtype=int)

    # Define 3 regions in the ground truth using circles
    rr1, cc1 = disk((30, 30), 15)
    rr2, cc2 = disk((70, 30), 15)
    rr3, cc3 = disk((50, 70), 20)

    ground_truth[rr1, cc1] = 1  # Region 1
    ground_truth[rr2, cc2] = 2  # Region 2
    ground_truth[rr3, cc3] = 3  # Region 3

    # Create an automated segmentation that imperfectly matches ground truth
    automated = np.zeros(shape, dtype=int)

    # Shift regions slightly or add noise
    rr1_auto, cc1_auto = disk((32, 32), 15)
    rr2_auto, cc2_auto = disk((68, 32), 15)
    rr3_auto, cc3_auto = disk((52, 68), 18)

    automated[rr1_auto, cc1_auto] = 1
    automated[rr2_auto, cc2_auto] = 2
    automated[rr3_auto, cc3_auto] = 3

    multi_region_dice(automated,ground_truth)
    # plt.imshow(compute_cost(automated,ground_truth))
    # Visualize both
    fig, axs = plt.subplots(1, 2, figsize=(8, 4))
    axs[0].imshow(ground_truth, cmap='tab10')
    axs[0].set_title('Ground Truth')
    axs[1].imshow(automated, cmap='tab10')
    axs[1].set_title('Automated Segmentation')
    for ax in axs: ax.axis('off')
    plt.tight_layout()
    plt.show()