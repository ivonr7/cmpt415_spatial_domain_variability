## Analyzing Variability in Domain Identification Methods

This repository contains all the methods for re-creating the results presented in [Variation in Spatial Domain Identification Methods](assets/SFU_CMPT415_Summer2025_ProjectReport_Isaac.pdf)

## Repository Layout
```python

├── assets  # images and pdfs in repo
├── notebooks # jupyter notebooks for waling through analyses
├── part2_annotators # code to run the SDI methods
├── scripts # Where the code to re-create results lives
├── src # methods to create and analyze the SDI masks
├── README.md # where you are :)
├── requirements.txt
├── setup.sh

```

## Installation
1. it is recommended to use venv for the python installation
    ```python -m venv venv && source venv/bin/activate```
2. install the [HEST-1k](https://github.com/mahmoodlab/HEST?tab=readme-ov-file#hest-library-installation) library 
3. Download the HEST-1k dataset (for the spatialLIBD benchmark the ids are [here](scripts/spatialLIBD_ids.csv) in the id column). Use the download [tutorial](https://github.com/mahmoodlab/HEST/blob/main/tutorials/1-Downloading-HEST-1k.ipynb)
4. Once those steps have been completed you can ```pip install -r requirements.txt``` and ```python setup.py```
5. running the [read_hest.py](scripts/misc/read_hest.py) file will verify the download  and install worked

## Usage
```python
from sdi_variation import cluster_matching,plot_cluster_match
import matplotlib.pyplot as plt
s1 = plt.imread("segmentation_mask_1.img")
s2 = plt.imread("segmentation_mask_2.img")

mrd_data, _ = cluster_matching(s1,s2)
print(f"Dice per Region: {"\n".join(mrd_data['dice'].values.tolist())}")
plot_cluster_match(s1,s2,mrd_data)
plt.show()
```


## Reproducing Results

We've added [test data](assets/test_data/MISC3_151674/) for re-producing our results. It's 1 sample from the spatialLIBD dataset with the columns for each method in the columns of the ```sample.adata.obs``` dataframe. This means we start from step 4 of the workflow I did. Running
```bash
python scripts/scripts/part4_create_segs.py assets/test_data/MISC3_151674/MISC3.h5 assets/test_data/MISC3_151674/
``` 
from the top level directory will create the folder ```assets/test_data/MISC3_151674/masks``` which holds the segmentations.  Next to get the MRDice outputs run again from the project directory
```bash
python scripts/part5_mrdice_methods.py assets/test_data/MISC3_151674/masks --unique
python scripts/part5_mrdice_methods.py assets/test_data/MISC3_151674/masks
``` 
This will create the ```assets/test_data/MISC3_151674/mrdice folder``` which contains mapped regions plots, csvs of the regions mapping and scores and the MRDice Agreement matrix.  To get the Normalized Mutual Information Plots run
```bash
python scripts/part5_nmi_methods.py assets/test_data/MISC3_151674/masks 
```

To reproduce the plots of Moran's I statistic and gene distribution use [this](notebooks/clustered_gene_dist.ipynb) to reproduce the investigation into the gene distribution of mapped regions to reproduce use [this](notebooks/pairwise_gene_dist.ipynb)

