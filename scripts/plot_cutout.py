import cartopy.crs as ccrs
from utils import load_and_limit_cutout, load_gdf_nuts, load_context_boundaries, plot_dataarray_on_map, resolve_user_home_path

import matplotlib
matplotlib.use('Agg')  # This enables backend without GUI (there seems to be problems with projection, PlateCarree)
import matplotlib.pyplot as plt

from typing import Any
snakemake: Any  # This is to avoid my IDE to complain about snakemake variable not being defined, but it is actually defined when running the script with snakemake



##############################
# This script reads the cutout for a given region and resource, limits it to the region, computes the mean field over time, and plots it on a map with the NUTS boundaries.



############################## Unwrap relevant variables

##### params
cutout_params = snakemake.params["cutout_params"]
fig_params = snakemake.params["fig_params"]
##### input
file_gdf_NUTS = snakemake.input["gdf_NUTS"]
file_gdf_NUTS_ref = snakemake.input["gdf_NUTS_ref"]
##### output
file_map_cutout = snakemake.output["map_cutout"]
##### wildcards
cutout = snakemake.wildcards["cutout"]
year =snakemake.wildcards["year"]
region = snakemake.wildcards["region"]
resource = snakemake.wildcards["resource"]
nuts = snakemake.wildcards["nuts"]



############################## Operations

##### Load local geometry (black outline) and context boundaries (grey)
_, gdf_NUTS_local = load_gdf_nuts(file_gdf_NUTS, region)
gdf_context = load_context_boundaries(file_gdf_NUTS_ref, nuts, clip_to=gdf_NUTS_local)


##### Load and limit cutout
file_cutout = resolve_user_home_path(cutout_params[f"{cutout}_{year}"]["path"])
c = load_and_limit_cutout(file_cutout, gdf_NUTS_local)


##### Compute fields
if resource == 'onwind':
    field = c.data.wnd100m.mean(dim="time")
elif resource == 'solar':
    field = c.data.influx_direct.mean(dim='time') + c.data.influx_diffuse.mean(dim='time')
else: 
    raise ValueError(f"[plot_cutout] resource must be 'onwind' or 'solar', but received: {resource}")



############################## Create outputs

resolution = 'LR'   # Always Low Resolution
size = fig_params["sizes"][resolution]["size"]
linewidth = fig_params["sizes"][resolution]["linewidth"]
fontsize = fig_params["sizes"][resolution]["fontsize"]
dpi = fig_params["sizes"][resolution]["dpi"]

cmap = fig_params['cutout'][resource]['cmap']
units = fig_params['cutout'][resource]['units']


##### Make plot
plot_dataarray_on_map(
    data=field,
    gdf_context=gdf_context,
    gdf_NUTS_local=gdf_NUTS_local,
    region=region,
    file_output=file_map_cutout,
    x_coord="lon",
    y_coord="lat",
    cmap=cmap,
    cbar_label=units,
    vmin=field.min().values.item(),
    vmax=field.max().values.item(),
    size=size,
    linewidth=linewidth,
    fontsize=fontsize,
    dpi=dpi,
    #title=
    bounds_type="gdf"
)


