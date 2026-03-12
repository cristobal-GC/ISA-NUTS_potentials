import xarray as xr
from utils import log_xarray_spatial_info

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
CAPACITY = xr.open_dataarray(file_nc_CAPACITY)
CF = xr.open_dataarray(file_nc_CF)
log_xarray_spatial_info(CAPACITY, source_label=file_nc_CAPACITY)
log_xarray_spatial_info(CF, source_label=file_nc_CF)

##### Align by shared coordinates (not necessary if we are sure they already align)
#CAPACITY, CF = xr.align(CAPACITY, CF, join="inner")

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
