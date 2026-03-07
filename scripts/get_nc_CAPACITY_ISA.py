import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import atlite
from atlite.gis import shape_availability
from atlite.gis import ExclusionContainer
from utils import load_and_limit_cutout, load_gdf_nuts_and_local

import rasterio
from rasterio.plot import show

from typing import Any
snakemake: Any  # This is to avoid my IDE to complain about snakemake variable not being defined, but it is actually defined when running the script with snakemake



############################## Unwrap relevant variables

##### params
cutout_params = snakemake.params["cutout_params"]
CF_params = snakemake.params["CF_params"]
##### input
file_gdf_NUTS = snakemake.input["gdf_NUTS"]
file_nc_CF = snakemake.input["nc_CF"]
file_raster_ISA = snakemake.input["raster_ISA"]
##### output
file_CAPACITY_ISA = snakemake.output["file_CAPACITY_ISA"]
file_CAPACITY_CF = snakemake.output["file_CAPACITY_CF"]
file_CAPACITY_CF_ISA = snakemake.output["file_CAPACITY_CF_ISA"]
##### wildcards
cutout = snakemake.wildcards["cutout"]
year =snakemake.wildcards["year"]
region = snakemake.wildcards["region"]
resource = snakemake.wildcards["resource"]
ISA = snakemake.wildcards["ISA"]



############################## Operations

crs_area = 3035


##### Load gdf_NUTS and set crs_area
gdf_NUTS, gdf_NUTS_local = load_gdf_nuts_and_local(file_gdf_NUTS, region)
gdf_NUTS = gdf_NUTS.to_crs(crs_area)
gdf_NUTS_local = gdf_NUTS_local.to_crs(crs_area)


##### Load and limit cutout, and set crs_area
file_cutout = cutout_params[f"{cutout}_{year}"]["path"]
c = load_and_limit_cutout(file_cutout, gdf_NUTS_local)
c.to_crs(crs_area)



##### Create excluder containers for ISA code, with 25x25 m resolution
excluder_ISA = ExclusionContainer(crs=crs_area, res=25)

### Define shape from the region geometry and excluder crs
shape = gdf_NUTS_local.geometry


##### Get A_ISA
# Add ISA code to the excluder_ISA, 
# codes must go in a list, otherwise, the zero index works wrongly
# option 1: crs=ISA_raster.crs (the object)
# option 2: crs=ISA_raster.crs.to_epsg() (just the number, 25830)
# Ambas opciones dan la misma suma de A, incluso también poner crs=3035, como tenía inicialmente
with rasterio.open(file_raster_ISA) as raster_ISA:
    excluder_ISA.add_raster(raster_ISA, codes=[ISA], crs=raster_ISA.crs, invert=True)
    input(f'ISA raster crs is {raster_ISA.crs}, excluder ISAcrs is {excluder_ISA.crs}. Press Enter to continue...')
# Compute availability matrix for ISA criterion
A_ISA = 100*c.availabilitymatrix(shape, excluder_ISA)


##### Generate CAPACITY matrix
# Generate area matrix, AREA, in km2
df = c.grid.to_crs(crs_area)
df['area'] = df.area*1e-6
AREA = xr.Dataset.from_dataframe(df.set_index(["y","x"]))["area"]
# Compute the CAPACITY matrix
CAPACITY_ISA    = AREA * CF_params[resource]["cap_per_sqkm"] * A_ISA * 0.01



############################## Create outputs
CAPACITY_ISA.to_netcdf(file_CAPACITY_ISA)





    







