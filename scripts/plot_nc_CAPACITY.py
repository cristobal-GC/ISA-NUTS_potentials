import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from utils import load_gdf_nuts, plot_dataarray_on_map, load_CAPACITY

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
CAPACITY = load_CAPACITY(file_nc_CAPACITY)

##### Load gdf_NUTS
gdf_NUTS, gdf_NUTS_local = load_gdf_nuts(file_gdf_NUTS, region)


############################## Create outputs

resolution = 'LR'   # Always Low Resolution
size = fig_params["sizes"][resolution]["size"]
linewidth = fig_params["sizes"][resolution]["linewidth"]
fontsize = fig_params["sizes"][resolution]["fontsize"]
dpi = fig_params["sizes"][resolution]["dpi"]

cmap = fig_params['CAPACITY'][resource]['cmap']
units = fig_params['CAPACITY'][resource]['units']

total_capacity_mw = CAPACITY.sum().values.item()
if total_capacity_mw > 1000:
    total_capacity_value = total_capacity_mw / 1000
    total_capacity_unit = "GW"
else:
    total_capacity_value = total_capacity_mw
    total_capacity_unit = "MW"


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
    dpi=dpi,
    title=f"Total capacity: {total_capacity_value:.2f} {total_capacity_unit}",
    bounds_type="gdf"
)
