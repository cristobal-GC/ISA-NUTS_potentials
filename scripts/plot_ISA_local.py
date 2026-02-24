import numpy as np
import geopandas as gpd
import yaml
import os
import pandas as pd

import rasterio
from rasterio.mask import mask
from rasterio.plot import show

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch



############################## Unwrap relevant variables

##### params
map_params = snakemake.params["map_params"]
##### input
file_gdf_NUTS = snakemake.input["gdf_NUTS"]
file_raster_ISA_local = snakemake.input["raster_ISA_local"]
##### output
file_map_ISA_local = snakemake.output["map_ISA_local"]
##### wildcards
region = snakemake.wildcards["region"]
resource = snakemake.wildcards["resource"]
resolution = snakemake.wildcards["resolution"]
format = snakemake.wildcards["format"]



############################## Operations

##### Define colors
colors_contraste = [    
    "#d73027",  #  0: Maximum risk (Red)
    "#fc8d59",  #  1: Very high (Orange)
    "#fee08b",  #  2: High (Yellow)
    "#91bfdb",  #  3: Moderate (Light blue)
    "#1a9850",  #  4: Low (Green)
    "#ffffff",  #  65536: Ocean
]
cmap = ListedColormap(colors_contraste)

##### Define classes
classes = [0, 1, 2, 3, 4, 65536]
boundaries = [-1.5, 0.5, 1.5, 2.5, 3.5, 4.5, 65599]  # boundaries between classes
norm = BoundaryNorm(boundaries, len(colors_contraste))



##### Load raster_ISA_local
raster_ISA_local = rasterio.open(file_raster_ISA_local)

##### Get ISA band
band = raster_ISA_local.read(1)



##### Get areas and percentages of each ISA code
# Unique values and frequencies
unique, counts = np.unique(band, return_counts=True)
df = pd.DataFrame({'value': unique, 'frequency': counts})
# Remove 65535
df = df[df['value'] != 65535]
# Make it sure all levels 0-4 exist
all_values = pd.DataFrame({'value': np.arange(5)})
df = all_values.merge(df, on='value', how='left').fillna(0)
# Assign area
df['area'] = df['frequency']*0.025*0.025
# Assign percentage
df['porc'] = 100*df['frequency'].div(df['frequency'].sum())

##### Add TOTAL row
# Compute sum of all numeric columns
totals = df.drop(columns=['value']).sum(numeric_only=True)
# Add row with totals
row_total = pd.DataFrame({**{'value': 'TOTAL'}, **totals.to_dict()}, index=[0])
# Add at the end
df = pd.concat([df, row_total], ignore_index=True)
# Round
df = df.round({'area': 2, 'porc': 2})

##### Guarda datos.....



##### Load gdf_NUTS and apply operations
gdf_NUTS = (gpd.read_file(file_gdf_NUTS)
           .set_index("NUTS_ID")            # set index 
           .to_crs(raster_ISA_local.crs)    # change crs to that of the ISA raster
)

# filter gdf with only one nuts region    
gdf_NUTS_local = gdf_NUTS.loc[[region]]



##### Define legend for both LR and HR maps:
labels = {    
        0: f"0: Maximum ({df.loc[0, 'porc']}%)",
        1: f"1: Very high ({df.loc[1, 'porc']}%)",
        2: f"2: High ({df.loc[2, 'porc']}%)",
        3: f"3: Moderate ({df.loc[3, 'porc']}%)",
        4: f"4: Low ({df.loc[4, 'porc']}%)",
        65536: "Ocean"
}

legend_elements = [
    Patch(facecolor=colors_contraste[i], edgecolor='black', label=labels[cls])
    for i, cls in enumerate(classes) if cls != 65536
]





############################## Create outputs

size = map_params[resolution]["size"]
linewidth = map_params[resolution]["linewidth"]
fontsize = map_params[resolution]["fontsize"]


##### Make plot
fig, ax = plt.subplots(figsize=(size, size))

show(band, transform=raster_ISA_local.transform, cmap=cmap, norm=norm, ax=ax)

##### Add gdf
gdf_NUTS.plot(ax=ax, color="none", edgecolor='grey', linewidth=linewidth)

##### Add gdf_local
gdf_NUTS_local.plot(ax=ax, color="none", edgecolor='black', linewidth=linewidth*2)

##### Add legend
loc = 'upper right'

ax.legend(
    handles=legend_elements, 
    title='ISA code',
    loc=loc, 
    fontsize=fontsize,
    title_fontsize=fontsize,
    frameon=True, 
    bbox_to_anchor=(1.5, 1),  # legend out of the plot    
)

plt.tight_layout()


##### Place the region
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