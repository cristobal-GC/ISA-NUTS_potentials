import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from utils import load_gdf_nuts, plot_dataarray_on_map, log_xarray_spatial_info

from typing import Any
snakemake: Any  # This is to avoid my IDE to complain about snakemake variable not being defined, but it is actually defined when running the script with snakemake



############################## Unwrap relevant variables

##### params
fig_params = snakemake.params["fig_params"]
##### input
file_gdf_NUTS = snakemake.input["gdf_NUTS"]
file_nc_CAPACITY = snakemake.input["nc_CAPACITY"]
##### output
file_map_CAPACITY = snakemake.output["map_CAPACITY"]
##### wildcards
cutout = snakemake.wildcards["cutout"]
year = snakemake.wildcards["year"]
region = snakemake.wildcards["region"]
resource = snakemake.wildcards["resource"]
try:
    ISA_list = [int(snakemake.wildcards["isa"])]
except (KeyError, AttributeError):
    ISA_list = [0, 1, 2, 3, 4]

############################## Operations

##### Load CAPACITY
CAPACITY = xr.open_dataarray(file_nc_CAPACITY)
log_xarray_spatial_info(CAPACITY, source_label=file_nc_CAPACITY)

##### Load gdf_NUTS
gdf_NUTS, gdf_NUTS_local = load_gdf_nuts(file_gdf_NUTS, region)


############################## Create outputs

resolution = 'LR'   # Always Low Resolution
size = fig_params["sizes"][resolution]["size"]
linewidth = fig_params["sizes"][resolution]["linewidth"]
fontsize = fig_params["sizes"][resolution]["fontsize"]

cmap = fig_params['CAPACITY'][resource]['cmap']
units = fig_params['CAPACITY'][resource]['units']


##### Make plot
plot_dataarray_on_map(
    data=CAPACITY,
    gdf_NUTS=gdf_NUTS,
    gdf_NUTS_local=gdf_NUTS_local,
    region=region,
    file_output=file_map_CAPACITY,
    x_coord="lon",
    y_coord="lat",
    cmap=cmap,
    cbar_label=units,
    #vmin=CAPACITY.min().values.item(),
    #vmax=CAPACITY.max().values.item(),
    size=size,
    linewidth=linewidth,
    fontsize=fontsize,
    title=f"Total capacity: {CAPACITY.sum().values.item():.2f} MW",
    bounds_type="data"
)
