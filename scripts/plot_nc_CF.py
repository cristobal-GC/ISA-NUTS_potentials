import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from utils import load_gdf_nuts_and_local, plot_dataarray_on_map

from typing import Any
snakemake: Any  # This is to avoid my IDE to complain about snakemake variable not being defined, but it is actually defined when running the script with snakemake



############################## Unwrap relevant variables

##### params
cutout_params = snakemake.params["cutout_params"]
fig_params = snakemake.params["fig_params"]
##### input
file_gdf_NUTS = snakemake.input["gdf_NUTS"]
file_nc_CF = snakemake.input["nc_CF"]
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

##### Load gdf_NUTS
gdf_NUTS, gdf_NUTS_local = load_gdf_nuts_and_local(file_gdf_NUTS, region)


############################## Create outputs

resolution = 'LR'   # Always Low Resolution
size = fig_params["sizes"][resolution]["size"]
linewidth = fig_params["sizes"][resolution]["linewidth"]
fontsize = fig_params["sizes"][resolution]["fontsize"]

cmap = fig_params['CF'][resource]['cmap']


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
    vmin=0,
    vmax=1,
    size=size,
    linewidth=linewidth,
    fontsize=fontsize,
    #title=
    bounds_type="gdf"
)