"""
One-time export of pickled train/test indices and protein chunk lists to CSV,
so the R elastic net scripts don't need reticulate to read Python pickles.
"""
import pickle
import pandas as pd
import glob
import os

train_index = pickle.load(open("Data/helper_data/train_index.pkl", "rb"))
test_index = pickle.load(open("Data/helper_data/test_index.pkl", "rb"))

pd.DataFrame({"id": list(train_index)}).to_csv("Data/helper_data/train_index.csv", index=False)
pd.DataFrame({"id": list(test_index)}).to_csv("Data/helper_data/test_index.csv", index=False)



# export every protlist_chunk*.pkl found, matching whatever batch scheme is already in use
## soma to olink direction
chunk_files = glob.glob("Data/helper_data/so_protlist_chunk*.pkl")
for f in chunk_files:
    batch_num = os.path.basename(f).replace("so_protlist_chunk", "").replace(".pkl", "")
    protlist = pickle.load(open(f, "rb"))
    pd.DataFrame({"protein": list(protlist)}).to_csv(
        f"Data/helper_data/so_protlist_chunk{batch_num}.csv", index=False
    )

## olink to soma direction
chunk_files = glob.glob("Data/helper_data/os_protlist_chunk*.pkl")
for f in chunk_files:
    batch_num = os.path.basename(f).replace("os_protlist_chunk", "").replace(".pkl", "")
    protlist = pickle.load(open(f, "rb"))
    pd.DataFrame({"protein": list(protlist)}).to_csv(
        f"Data/helper_data/os_protlist_chunk{batch_num}.csv", index=False
    )






