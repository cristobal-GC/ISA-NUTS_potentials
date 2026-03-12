import numpy as np
import pandas as pd

import rasterio
from rasterio.plot import show
from utils import load_gdf_nuts_and_local, log_raster_spatial_info

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

from typing import Any
snakemake: Any  # This is to avoid my IDE to complain about snakemake variable not being defined, but it is actually defined when running the script with snakemake



############################## Unwrap relevant variables

##### params
fig_params = snakemake.params["fig_params"]
##### input
file_gdf_NUTS = snakemake.input["gdf_NUTS"]
file_raster_ISA = snakemake.input["raster_ISA"]
file_df_ISA = snakemake.input["df_ISA"]
##### output
file_map_ISA = snakemake.output["map_ISA"]
##### wildcards
region = snakemake.wildcards["region"]
resource = snakemake.wildcards["resource"]
resolution = snakemake.wildcards["resolution"]
format = snakemake.wildcards["format"]



############################## Operations

##### Load raster_ISA and read required attributes
with rasterio.open(file_raster_ISA) as raster_ISA:
    log_raster_spatial_info(raster_ISA, source_label=file_raster_ISA)

    raster_crs = raster_ISA.crs
    transform = raster_ISA.transform
    band = raster_ISA.read(1)
    nodata = raster_ISA.nodata if raster_ISA.nodata is not None else 65535

##### Load df_ISA
df = pd.read_csv(file_df_ISA, index_col="value")

##### Load gdf_NUTS and change crs to that of the ISA raster 
gdf_NUTS, gdf_NUTS_local = load_gdf_nuts_and_local(file_gdf_NUTS, region)
gdf_NUTS = gdf_NUTS.to_crs(raster_crs)
gdf_NUTS_local = gdf_NUTS_local.to_crs(raster_crs)



##### Prepare for plotting

# Get colors and labels from config
colors_contraste = [fig_params["ISA"]["colors"][cls] for cls in range(5)]
cmap = ListedColormap(colors_contraste)

# Define legend:
labels = {
    cls: f"{cls}: {fig_params['ISA']['labels'][cls]} ({df.loc[str(cls), 'porc']:.2f}%)"
    for cls in range(5)
}

legend_elements = [
    Patch(facecolor=colors_contraste[cls], edgecolor='black', label=labels[cls])
    for cls in labels.keys()
]


# Mask nodata (use nodata read from raster if available)
band_masked = np.ma.masked_equal(band, nodata)





############################## Create outputs

size = fig_params["sizes"][resolution]["size"]
linewidth = fig_params["sizes"][resolution]["linewidth"]
fontsize = fig_params["sizes"][resolution]["fontsize"]

##### Make plot
fig, ax = plt.subplots(figsize=(size, size))

# It seems that current version of the show function normalises the band values between 0 and 1? I need to set vmin=0 and vmax=1 to get proper colour assesment
show(band_masked,
    transform=transform,
    cmap=cmap,
    vmin=0,
    vmax=1,
    ax=ax
)

# Add gdf (only regions for the same NUTS level than local)
gdf_NUTS[gdf_NUTS.index.astype(str).str.len() == len(region)].plot(ax=ax, color="none", edgecolor='grey', linewidth=linewidth)

# Add gdf_local with double linewidth
gdf_NUTS_local.plot(ax=ax, color="none", edgecolor='black', linewidth=linewidth*2)

# Add legend
ax.legend(
    handles=legend_elements, 
    title='ISA code',
    loc='upper right', 
    fontsize=fontsize,
    title_fontsize=fontsize,
    frameon=True, 
    bbox_to_anchor=(1.5, 1),  # legend out of the plot    
)

# Set limits using same logic as plot_dataarray_on_map
xmin, ymin, xmax, ymax = gdf_NUTS_local.total_bounds
km_per_lon = 85
km_per_lat = 111
center_x = (xmax + xmin) / 2
center_y = (ymax + ymin) / 2
delta_x = xmax - xmin
delta_y = ymax - ymin
delta_km = max([km_per_lon*delta_x, km_per_lat*delta_y]) * 1.02
ax.set_xlim(
    center_x - 0.5*delta_km/km_per_lon,
    center_x + 0.5*delta_km/km_per_lon
)
ax.set_ylim(
    center_y - 0.5*delta_km/km_per_lat,
    center_y + 0.5*delta_km/km_per_lat
)

ax.tick_params(axis="both", labelsize=fontsize*0.8)
ax.set_xlabel("Lon", fontsize=fontsize*0.8)
ax.set_ylabel("Lat", fontsize=fontsize*0.8)


##### Save figure
fig.savefig(file_map_ISA,
            bbox_inches="tight",
            pad_inches=0.2)  

plt.close(fig)