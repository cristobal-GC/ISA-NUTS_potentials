import cartopy.crs as ccrs
from utils import load_and_limit_cutout, load_gdf_nuts, plot_dataarray_on_map, resolve_user_home_path

import matplotlib
matplotlib.use('Agg')  # This enables backend without GUI (there seems to be problems with projection, PlateCarree)
import matplotlib.pyplot as plt

from typing import Any
snakemake: Any  # This is to avoid my IDE to complain about snakemake variable not being defined, but it is actually defined when running the script with snakemake



############################## Unwrap relevant variables

##### params
cutout_params = snakemake.params["cutout_params"]
fig_params = snakemake.params["fig_params"]
##### input
file_gdf_NUTS = snakemake.input["gdf_NUTS"]
##### output
file_map_cutout = snakemake.output["map_cutout"]
##### wildcards
cutout = snakemake.wildcards["cutout"]
year =snakemake.wildcards["year"]
region = snakemake.wildcards["region"]
resource = snakemake.wildcards["resource"]



############################## Operations

##### Load gdf_NUTS and gdf_NUTS_local
gdf_NUTS, gdf_NUTS_local = load_gdf_nuts(file_gdf_NUTS, region)


##### Load and limit cutout
file_cutout = resolve_user_home_path(cutout_params[f"{cutout}_{year}"]["path"])
c = load_and_limit_cutout(file_cutout, gdf_NUTS_local)


##### Compute fields
if resource == 'onwind':
    field = c.data.wnd100m.mean(dim="time")
elif resource == 'solar':
    field = c.data.influx_direct.mean(dim='time') + c.data.influx_diffuse.mean(dim='time')
else: 
    raise ValueError(f"resource must be 'onwind' or 'solar', but received: {resource}")



############################## Create outputs

resolution = 'LR'   # Always Low Resolution
size = fig_params["sizes"][resolution]["size"]
linewidth = fig_params["sizes"][resolution]["linewidth"]
fontsize = fig_params["sizes"][resolution]["fontsize"]

cmap = fig_params['cutout'][resource]['cmap']
units = fig_params['cutout'][resource]['units']


##### Make plot
plot_dataarray_on_map(
    data=field,
    gdf_NUTS=gdf_NUTS,
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
    #title=
    bounds_type="gdf"
)


