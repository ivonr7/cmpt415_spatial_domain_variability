from scipy.stats import entropy
import numpy as np


'''
    Mutual Information 
    Calculated based on registration lecture from cmpt415
    Calculated sparsely by computing counts of 1D indicies of image pairs

    Works only for grayscale Images 
'''
def mutual_info(img1:np.ndarray,img2:np.ndarray):
    assert img1.shape == img2.shape
    # Compute flattened indexes
    brightest = max(img1.max(),img2.max())
    idx = img1.flatten() * brightest + img2.flatten() 
    j_hist = np.bincount(idx)
    return entropy(j_hist.flatten() / j_hist.sum())


'''
    Normalized Mutual Information by max marginal entropy 
'''
def normalized_mi(img1:np.ndarray,img2:np.ndarray):
    e1 = entropy(img1.flatten())
    e2 = entropy(img2.flatten())
    return mutual_info(img1,img2) / max(e1,e2)

if __name__ == "__main__":
    np.random.seed(7)
    res = 255
    im2 = np.random.uniform(low=0,high=65535,size=(res,res)).astype(np.uint16)
    im1 = np.random.uniform(low=0,high=65535,size=(res,res)).astype(np.uint16)
    print(mutual_info(im1,im2))


    im2 = np.random.uniform(low=0,high=100,size=(res,res)).astype(np.uint16)
    im1 = im2.copy()
    print(mutual_info(im1,im2))
