import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import logging
import re
from pathlib import Path
from utils import load_gdf_nuts_and_local, plot_dataarray_on_map, log_xarray_spatial_info

from typing import Any
snakemake: Any  # This is to avoid my IDE to complain about snakemake variable not being defined, but it is actually defined when running the script with snakemake


logger = logging.getLogger(__name__)


def _log_and_print(message):
    logger.info(message)
    print(message)



############################## Unwrap relevant variables

##### params
fig_params = snakemake.params["fig_params"]
##### input
file_gdf_NUTS = snakemake.input["gdf_NUTS"]
file_nc_CF = snakemake.input["nc_CF"]
files_df_CF_CAPACITY = snakemake.input["dfs_CF_CAPACITY"]
##### output
file_map_CF = snakemake.output["map_CF"]
##### wildcards
cutout = snakemake.wildcards["cutout"]
year =snakemake.wildcards["year"]
region = snakemake.wildcards["region"]
resource = snakemake.wildcards["resource"]


############################## Operations

##### Load CF
CF = xr.open_dataarray(file_nc_CF)
log_xarray_spatial_info(CF, source_label=file_nc_CF)

##### Load gdf_NUTS
gdf_NUTS, gdf_NUTS_local = load_gdf_nuts_and_local(file_gdf_NUTS, region)


############################## Create outputs

resolution = 'LR'   # Always Low Resolution
size = fig_params["sizes"][resolution]["size"]
linewidth = fig_params["sizes"][resolution]["linewidth"]
fontsize = fig_params["sizes"][resolution]["fontsize"]

cmap = fig_params['CF'][resource]['cmap']

# Compute common CF limits from all df_CF_CAPACITY inputs.
# These files include ISA0..ISA4 and all regions for the same
# {cutout, nuts, resource, year}, so each regional CF map shares
# the same absolute color scale (vmin/vmax).
vmin = None
vmax = None

for file_path in files_df_CF_CAPACITY:
    df = pd.read_csv(file_path, usecols=["CF"])
    if df.empty:
        continue

    file_vmin = float(df["CF"].min())
    file_vmax = float(df["CF"].max())

    vmin = file_vmin if vmin is None else min(vmin, file_vmin)
    vmax = file_vmax if vmax is None else max(vmax, file_vmax)

if vmin is None or vmax is None:
    raise ValueError("Could not compute vmin/vmax from dfs_CF_CAPACITY: no CF values found.")


# Log CF limits and regions included in the plot (extracted from dfs_CF_CAPACITY filenames)
region_pattern = re.compile(
    rf"df_CF_CAPACITY_ISA\d+_{re.escape(resource)}_(.+)_{re.escape(str(year))}\.csv$"
)
regions = []
for file_path in files_df_CF_CAPACITY:
    match = region_pattern.match(Path(file_path).name)
    if match:
        regions.append(match.group(1))

regions = sorted(set(regions))

_log_and_print(
    f"[plot_nc_CF] CF min/max for cutout={cutout}, year={year}, resource={resource}, regions={regions}: vmin={vmin:.3f}, vmax={vmax:.3f}"
)



##### Make plot
plot_dataarray_on_map(
    data=CF,
    gdf_NUTS=gdf_NUTS,
    gdf_NUTS_local=gdf_NUTS_local,
    region=region,
    file_output=file_map_CF,
    x_coord="lon",
    y_coord="lat",
    cmap=cmap,
    cbar_label="Capacity Factor",
    vmin=vmin,
    vmax=vmax,
    size=size,
    linewidth=linewidth,
    fontsize=fontsize,
    #title=
    bounds_type="gdf"
)