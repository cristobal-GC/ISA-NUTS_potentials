import geopandas as gpd
import cartopy.crs as ccrs
import atlite

import matplotlib
matplotlib.use('Agg')  # This enables backend without GUI (there seems to be problems with projection, PlateCarree)
import matplotlib.pyplot as plt




############################## Unwrap relevant variables

##### params
cutout_params = snakemake.params["cutout_params"]
map_params = snakemake.params["map_params"]
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

##### Load cutout 
file_cutout = cutout_params[f"{cutout}_{year}"]["path"]
c = atlite.Cutout(file_cutout)

##### Load gdf_NUTS and apply operations
gdf_NUTS = (gpd.read_file(file_gdf_NUTS)
           .set_index("NUTS_ID")            # set index 
           .to_crs(c.crs)                   # change crs to that of the cutout
)

##### Prepare gdf: filter gdf with only one nuts region    
gdf_NUTS_local = gdf_NUTS.loc[[region]]


##### limit cutout
# Get region bounding box
xmin, ymin, xmax, ymax = gdf_NUTS_local.total_bounds

# Add margin of one cell
margin_x = c.dx
margin_y = c.dy

xmin -= margin_x
xmax += margin_x
ymin -= margin_y
ymax += margin_y

# Limit cutout
c = c.sel(bounds=(xmin, ymin, xmax, ymax))


##### Compute fields
if resource == 'onwind':
    field = c.data.wnd100m.mean(dim="time")
elif resource == 'solar':
    field = c.data.influx_direct.mean(dim='time') + c.data.influx_diffuse.mean(dim='time')
else: 
    raise ValueError(f"resource must be 'onwind' or 'solar', but received: {resource}")



############################## Create outputs

resolution = 'LR'   # Always Low Resolution
size = map_params[resolution]["size"]
linewidth = map_params[resolution]["linewidth"]
fontsize = map_params[resolution]["fontsize"]

cmap = map_params['cutout'][resource]['cmap']
units = map_params['cutout'][resource]['units']


##### Make plot
fig, ax = plt.subplots(figsize=(size, size), subplot_kw={"projection": ccrs.PlateCarree()})

# Plot field without automatic cbar
mappable = field.plot(ax=ax,
                      cmap=cmap,
                      vmin=field.min().values.item(),
                      vmax=field.max().values.item(),
                      add_colorbar=False                     
)

 # Add colorbar out of the figure
cbar = fig.colorbar(mappable,
                    ax=ax,
                    orientation='vertical',
                    fraction=0.046,
                    pad=0.04   # control size and gap
)

cbar.set_label(units, fontsize=fontsize)
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
fig.savefig(file_map_cutout,
            bbox_inches="tight",
            pad_inches=0)  

plt.close(fig)



