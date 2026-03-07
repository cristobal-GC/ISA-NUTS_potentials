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


##### Load CF
CF = xr.open_dataarray(file_nc_CF)

##### Validate CF coordinates match cutout grid
# This is important for filtering the CF with the ISA excluder, which is based on the cutout grid. If they don't match, the filtering will not work and the resulting A_CF and A_CF_ISA will be wrong.
print(f"[get_nc_As] Validating coordinate alignment:")
print(f"  Cutout x: [{c.data.x.min().values:.2f}, {c.data.x.max().values:.2f}], shape: {len(c.data.x)}")
print(f"  Cutout y: [{c.data.y.min().values:.2f}, {c.data.y.max().values:.2f}], shape: {len(c.data.y)}")
print(f"  CF x: [{CF.x.min().values:.2f}, {CF.x.max().values:.2f}], shape: {len(CF.x)}")
print(f"  CF y: [{CF.y.min().values:.2f}, {CF.y.max().values:.2f}], shape: {len(CF.y)}")

coords_match = CF.x.equals(c.data.x) and CF.y.equals(c.data.y)
if not coords_match:
    input(f"[get_nc_As] Coordinates mismatch detected. Press ENTER for reindexing CF to cutout grid.")
    CF = CF.sel(x=c.data.x, y=c.data.y, method="nearest")
else:
    input(f"[get_nc_As] Coordinates match perfectly.")



##### Create excluder containers for ISA code and for all ISA codes, with 25x25 m resolution
excluder_ISA = ExclusionContainer(crs=crs_area, res=25)
excluder_allISA = ExclusionContainer(crs=crs_area, res=25)

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


##### Get A_CF
# Add all ISA codes to the excluder_ISA, 
# codes must go in a list, otherwise, the zero index works wrongly
# option 1: crs=ISA_raster.crs (the object)
# option 2: crs=ISA_raster.crs.to_epsg() (just the number, 25830)
# Ambas opciones dan la misma suma de A, incluso también poner crs=3035, como tenía inicialmente
with rasterio.open(file_raster_ISA) as raster_ISA:
    excluder_allISA.add_raster(raster_ISA, codes=[0, 1, 2, 3, 4], crs=raster_ISA.crs, invert=True)
    input(f'ISA raster crs is {raster_ISA.crs}, excluder allISA crs is {excluder_allISA.crs}. Press Enter to continue...')
# Compute availability matrix for CF criterion
A_CF = 100*(
    c.availabilitymatrix(shape, excluder_allISA).where(CF >= CF_params[resource]["CF_threshold"], 0)
)


##### Get A_ISA_CF
# Compute availability matrix for ISA and CFcriterion
A_CF_ISA = 100*(
    c.availabilitymatrix(shape, excluder_ISA).where(CF >= CF_params[resource]["CF_threshold"], 0)
)


##### Generate CAPACITY matrices 
# Generatearea matrix, AREA, in km2
df = c.grid.to_crs(crs_area)
df['area'] = df.area*1e-6
AREA = xr.Dataset.from_dataframe(df.set_index(["y","x"]))["area"]


### Compute the CAPACITY matrices
CAPACITY_ISA    = AREA * CF_params[resource]["cap_per_sqkm"] * A_ISA * 0.01
CAPACITY_CF     = AREA * CF_params[resource]["cap_per_sqkm"] * A_CF * 0.01
CAPACITY_CF_ISA = AREA * CF_params[resource]["cap_per_sqkm"] * A_CF_ISA * 0.01


############################## Create outputs
CAPACITY_ISA.to_netcdf(file_CAPACITY_ISA)
CAPACITY_CF.to_netcdf(file_CAPACITY_CF)
CAPACITY_CF_ISA.to_netcdf(file_CAPACITY_CF_ISA)







    







