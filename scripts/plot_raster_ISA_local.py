import numpy as np
import geopandas as gpd
import pandas as pd

import rasterio
from rasterio.plot import show

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch



############################## Unwrap relevant variables

##### params
map_params = snakemake.params["map_params"]
##### input
file_gdf_NUTS = snakemake.input["gdf_NUTS"]
file_raster_ISA_local = snakemake.input["raster_ISA_local"]
file_df_ISA_local = snakemake.input["df_ISA_local"]
##### output
file_map_ISA_local = snakemake.output["map_ISA_local"]
##### wildcards
region = snakemake.wildcards["region"]
resource = snakemake.wildcards["resource"]
resolution = snakemake.wildcards["resolution"]
format = snakemake.wildcards["format"]



############################## Operations

##### Load raster_ISA_local
raster_ISA_local = rasterio.open(file_raster_ISA_local)

##### Load df_ISA_local
df = pd.read_csv(file_df_ISA_local, index_col="value")

##### Load gdf_NUTS and apply operations
gdf_NUTS = (gpd.read_file(file_gdf_NUTS)
           .set_index("NUTS_ID")            # set index 
           .to_crs(raster_ISA_local.crs)    # change crs to that of the ISA raster
)



##### Prepare raster

# Get ISA band
band = raster_ISA_local.read(1)

# Define colors
colors_contraste = [    
    "#d73027",  #  0: Maximum risk (Red)
    "#fc8d59",  #  1: Very high (Orange)
    "#fee08b",  #  2: High (Yellow)
    "#91bfdb",  #  3: Moderate (Light blue)
    "#1a9850",  #  4: Low (Green)
]
cmap = ListedColormap(colors_contraste)

# Define legend:
labels = {    
        0: f"0: Maximum ({df.loc['0', 'porc']}%)",
        1: f"1: Very high ({df.loc['1', 'porc']}%)",
        2: f"2: High ({df.loc['2', 'porc']}%)",
        3: f"3: Moderate ({df.loc['3', 'porc']}%)",
        4: f"4: Low ({df.loc['4', 'porc']}%)",
}

legend_elements = [
    Patch(facecolor=colors_contraste[cls], edgecolor='black', label=labels[cls])
    for cls in labels.keys()
]


# Mask nodata
nodata = 65535
band_masked = np.ma.masked_equal(band, nodata)



##### Prepare gdf: filter gdf with only one nuts region    
gdf_NUTS_local = gdf_NUTS.loc[[region]]



############################## Create outputs

size = map_params[resolution]["size"]
linewidth = map_params[resolution]["linewidth"]
fontsize = map_params[resolution]["fontsize"]

##### Make plot
fig, ax = plt.subplots(figsize=(size, size))

# It seems that current version of the show function normalises the band values between 0 and 1? I need to set vmin=0 and vmax=1 to get proper colour assesment
show(band_masked,
     transform=raster_ISA_local.transform,
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

# Tight layout
plt.tight_layout()


# Set limits
xmin, ymin, xmax, ymax = gdf_NUTS_local.total_bounds
center_x = (xmax+xmin)/2
center_y = (ymax+ymin)/2
delta_x = xmax-xmin
delta_y = ymax-ymin
delta = max([delta_x, delta_y])
ax.set_xlim(center_x-0.51*delta, center_x+0.51*delta)
ax.set_ylim(center_y-0.51*delta, center_y+0.51*delta)

ax.set_xticks([])
ax.set_yticks([])


##### Save figure
fig.savefig(file_map_ISA_local,
            bbox_inches="tight",
            pad_inches=0)  

plt.close(fig)