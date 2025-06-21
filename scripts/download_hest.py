import datasets
import pandas as pd


ids = pd.read_csv("spatialLIBD_ids.csv")['id'].values.tolist()
local_dir = '../hest_data'

print(ids)
list_patterns = [f"*{id}[_.]**" for id in ids]

# Note that the full dataset is around 1TB of data
dataset = datasets.load_dataset(
    'MahmoodLab/hest', 
    cache_dir=local_dir,
    patterns=list_patterns
)

