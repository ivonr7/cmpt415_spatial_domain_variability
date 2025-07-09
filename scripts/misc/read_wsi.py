from PIL import Image
import openslide
import os
import matplotlib.pyplot as plt

data_dir='../hest_data'
wsi_dir = os.path.join(data_dir, 'wsis')

# full res. WSI image path of a given sample, (e.g., 'NCBI641')
id = 'MISC12'
wsi_filename = os.path.join(wsi_dir, id + '.tif')

wsi_img_PIL = Image.open(wsi_filename)
wsi_img_openslide = openslide.OpenSlide(wsi_filename)
plt.imshow(wsi_img_PIL)
plt.show()