import xarray as xr
from utils import load_CF, load_CAPACITY

from typing import Any
snakemake: Any


############################## Unwrap relevant variables

##### input
file_nc_CAPACITY = snakemake.input["nc_CAPACITY"]
file_nc_CF = snakemake.input["nc_CF"]
##### output
file_df_CF_CAPACITY = snakemake.output["df_CF_CAPACITY"]


############################## Operations

##### Load inputs
CAPACITY = load_CAPACITY(file_nc_CAPACITY)
CF = load_CF(file_nc_CF)


##### Build dataset and dataframe
ds = xr.Dataset({
    "CAPACITY": CAPACITY,
    "CF": CF,
})

df = (
    ds.to_dataframe()
    .reset_index()
    .loc[:, ["CF", "CAPACITY"]]
    .dropna(subset=["CF", "CAPACITY"])
    .loc[lambda d: d["CAPACITY"] > 0]
    .sort_values("CF", ascending=False)
)

df["CAPACITY"] = df["CAPACITY"].round(6)
df["CUM_CAPACITY"] = df["CAPACITY"].cumsum().round(6)


############################## Create outputs

df.to_csv(file_df_CF_CAPACITY, index=False)
