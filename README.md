## Analyzing Variability in Domain Identification Methods



## Installation
1. it is recommended to use venv from the python installation
    ```python -m venv venv && source venv/bin/activate```
2. install the [HEST-1k](https://github.com/mahmoodlab/HEST?tab=readme-ov-file#hest-library-installation) library 
3. Download the HEST-1k dataset (for the spatialLIBD benchmark the ids are [here](spatialLIBD_ids.csv) in the id column). Use the download [tutorial](https://github.com/mahmoodlab/HEST/blob/main/tutorials/1-Downloading-HEST-1k.ipynb)
4. Once those steps have been completed you can ```pip install -r requirements.txt``` and ```python setup.py```
5. running the [read_hest.py](scripts/read_hest.py) file will verify the download  and install worked
