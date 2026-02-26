import xarray as xr
import matplotlib.pyplot as plt
import geopandas as gpd
import cartopy.crs as ccrs

from typing import Any
snakemake: Any  # This is to avoid my IDE to complain about snakemake variable not being defined, but it is actually defined when running the script with snakemake



############################## Unwrap relevant variables

##### params
cutout_params = snakemake.params["cutout_params"]
map_params = snakemake.params["map_params"]
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

##### Load gdf_NUTS and apply operations
gdf_NUTS = (
    gpd.read_file(file_gdf_NUTS)
    .set_index("NUTS_ID")            # set index 
    .to_crs(4326)                    # set 4326
)

##### Filter gdf with only one nuts region    
gdf_NUTS_local = gdf_NUTS.loc[[region]]



############################## Create outputs

resolution = 'LR'   # Always Low Resolution
size = map_params[resolution]["size"]
linewidth = map_params[resolution]["linewidth"]
fontsize = map_params[resolution]["fontsize"]

cmap = map_params['CF'][resource]['cmap']


##### Make plot
fig, ax = plt.subplots(figsize=(size, size))

# Plot field without automatic cbar
mappable = CF.plot(
    ax=ax,
    x="lon",
    y="lat",
    cmap=cmap,
    vmin=0,
    vmax=1,
    add_colorbar=False
)

ax.tick_params(axis="both", labelsize=fontsize*0.8)
ax.set_xlabel("Longitude", fontsize=fontsize*0.8)
ax.set_ylabel("Latitude", fontsize=fontsize*0.8)

 # Add colorbar out of the figure
cbar = fig.colorbar(
    mappable,
    ax=ax,
    orientation='vertical',
    fraction=0.046,
    pad=0.04   # control size and gap
)

cbar.set_label("Capacity Factor", fontsize=fontsize)
cbar.ax.tick_params(labelsize=fontsize)


# Add gdf (only regions for the same NUTS level than local)
gdf_NUTS[gdf_NUTS.index.astype(str).str.len() == len(region)].plot(ax=ax, color="none", edgecolor='grey', linewidth=linewidth)

# Add gdf_local with double linewidth
gdf_NUTS_local.plot(ax=ax, color="none", edgecolor='black', linewidth=linewidth*2)

# Set limits
xmin, ymin, xmax, ymax = gdf_NUTS_local.total_bounds
km_per_lon = 85
km_per_lat = 111
center_x = (xmax+xmin)/2
center_y = (ymax+ymin)/2
delta_x = xmax-xmin
delta_y = ymax-ymin        
delta_km = max([km_per_lon*delta_x, km_per_lat*delta_y])*1.02
ax.set_xlim(center_x-0.5*delta_km/km_per_lon, center_x+0.5*delta_km/km_per_lon)
ax.set_ylim(center_y-0.5*delta_km/km_per_lat, center_y+0.5*delta_km/km_per_lat)


##### Save figure
fig.savefig(
    file_map_CF,
    bbox_inches="tight",
    pad_inches=0.2
)  

plt.close(fig)