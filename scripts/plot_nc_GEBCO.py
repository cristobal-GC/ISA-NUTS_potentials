import numpy as np
import xarray as xr
from shapely import contains_xy
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap
from utils import load_gdf_nuts, plot_dataarray_on_map, load_GEBCO

from typing import Any
snakemake: Any  # This is to avoid my IDE to complain about snakemake variable not being defined, but it is actually defined when running the script with snakemake



############################## Unwrap relevant variables

##### params
fig_params = snakemake.params["fig_params"]
##### input
file_gdf_NUTS = snakemake.input["gdf_NUTS"]
file_nc_GEBCO = snakemake.input["nc_GEBCO"]
##### output
file_map_GEBCO = snakemake.output["map_GEBCO"]
##### wildcards
region = snakemake.wildcards["region"]



############################## Operations

##### Load gdf_NUTS
gdf_NUTS, gdf_NUTS_local = load_gdf_nuts(file_gdf_NUTS, region)

##### Load GEBCO
GEBCO = load_GEBCO(file_nc_GEBCO)



############################## Prepare regional subset

x_coord = "lon"
y_coord = "lat"

if x_coord not in GEBCO.coords or y_coord not in GEBCO.coords:
    raise ValueError(
        f"[plot_nc_GEBCO] Expected GEBCO coordinates '{x_coord}' and '{y_coord}', found {list(GEBCO.coords)}"
    )

xmin, ymin, xmax, ymax = gdf_NUTS_local.total_bounds

margin_x = max((xmax - xmin) * 0.05, 0.25)
margin_y = max((ymax - ymin) * 0.05, 0.25)

lon_slice = slice(xmin - margin_x, xmax + margin_x)

lat_values = GEBCO[y_coord].values
if lat_values[0] <= lat_values[-1]:
    lat_slice = slice(ymin - margin_y, ymax + margin_y)
else:
    lat_slice = slice(ymax + margin_y, ymin - margin_y)

GEBCO_local = GEBCO.sel({x_coord: lon_slice, y_coord: lat_slice})

# Keep only cells inside the requested region geometry.
region_geom = gdf_NUTS_local.geometry.union_all()
lon_vals = GEBCO_local[x_coord].values
lat_vals = GEBCO_local[y_coord].values
lon_2d, lat_2d = np.meshgrid(lon_vals, lat_vals)

inside_region = contains_xy(region_geom, lon_2d, lat_2d)
inside_region_da = xr.DataArray(
    inside_region,
    coords={y_coord: GEBCO_local[y_coord], x_coord: GEBCO_local[x_coord]},
    dims=(y_coord, x_coord),
)
GEBCO_local = GEBCO_local.where(inside_region_da)



############################## Create outputs

resolution = "LR"
size = fig_params["sizes"][resolution]["size"]
linewidth = fig_params["sizes"][resolution]["linewidth"]
fontsize = fig_params["sizes"][resolution]["fontsize"]

# Reuse the classic terrain palette while removing the initial blue segment,
# so the scale starts at green and keeps yellow-brown-white at higher values.
terrain_no_blue = colormaps["terrain"](np.linspace(0.35, 1.0, 256))
gebco_cmap = LinearSegmentedColormap.from_list(
    "terrain_no_blue",
    terrain_no_blue,
)



##### Make plot
plot_dataarray_on_map(
    data=GEBCO_local,
    gdf_NUTS=gdf_NUTS,
    gdf_NUTS_local=gdf_NUTS_local,
    region=region,
    file_output=file_map_GEBCO,
    x_coord=x_coord,
    y_coord=y_coord,
    cmap=gebco_cmap,
    cbar_label="Elevation [m]",
    vmin=0,
    vmax=3500,
    size=size,
    linewidth=linewidth,
    fontsize=fontsize,
    title=None,
    bounds_type="gdf"
)