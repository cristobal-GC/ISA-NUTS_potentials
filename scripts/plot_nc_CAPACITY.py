import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from utils import load_gdf_nuts_and_local, plot_dataarray_on_map

from typing import Any
snakemake: Any  # This is to avoid my IDE to complain about snakemake variable not being defined, but it is actually defined when running the script with snakemake



############################## Unwrap relevant variables

##### params
map_params = snakemake.params["map_params"]
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
isa = snakemake.wildcards["isa"]


############################## Operations

##### Load CAPACITY
CAPACITY = xr.open_dataarray(file_nc_CAPACITY)

##### Load gdf_NUTS
gdf_NUTS, gdf_NUTS_local = load_gdf_nuts_and_local(file_gdf_NUTS, region)


############################## Create outputs

resolution = 'LR'   # Always Low Resolution
size = map_params[resolution]["size"]
linewidth = map_params[resolution]["linewidth"]
fontsize = map_params[resolution]["fontsize"]

cmap = map_params['CAPACITY'][resource]['cmap']
units = map_params['CAPACITY'][resource]['units']


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
    #title=f"CAPACITY (ISA code: {isa})",
    bounds_type="data"
)
