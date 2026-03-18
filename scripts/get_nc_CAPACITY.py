import xarray as xr
import logging
import rasterio

from atlite.gis import ExclusionContainer
from utils import (
    load_CF,
    load_and_limit_cutout,
    load_gdf_nuts,
    resolve_user_home_path,
    log_raster_spatial_info,
)

from typing import Any
snakemake: Any  # This is to avoid my IDE to complain about snakemake variable not being defined, but it is actually defined when running the script with snakemake


logger = logging.getLogger(__name__)


def _log_and_print(message):
    logger.info(message)
    print(message)



def get_CAPACITY_matrix(
        file_raster_ISA,
        ISA_list,
        gdf_NUTS_local,
        c,
        CF,
        CF_threshold,
        cap_per_sqkm,
):
    """Compute the CAPACITY matrix for a given region, combining two criteria: ISA_codes and CF threshold."""

    if gdf_NUTS_local.crs is None:
        raise ValueError("gdf_NUTS_local must have a defined CRS.")
    
    ##### Open ISA raster, jsut to retrieve spatial info and log it, and to create the excluder with the correct resolution and CRS info. The actual masking will be done internally by atlite when building the excluder mask, by passing the raster path and its CRS to the add_raster method, so we don't need to keep the raster open after creating the excluder.
    with rasterio.open(file_raster_ISA) as raster_ISA:
        
        # log spatial info
        log_raster_spatial_info(raster_ISA, source_label=file_raster_ISA)

        # Check if the raster has a CRS defined
        if raster_ISA.crs is None:
            raise ValueError(f"Raster {file_raster_ISA} has no CRS. Cannot build ExclusionContainer mask reliably.")
        else:
            raster_crs = raster_ISA.crs

        ##### Create excluder container, with 25x25 m resolution
        excluder = ExclusionContainer(res=raster_ISA.res[0]) # no need to specify CRS here, by default it is 3035, and atlite will handle reprojection internally when adding the raster criterion, as long as we provide the correct input raster CRS.
            


    _log_and_print(
        f"[get_CAPACITY_matrix] CRS info -> gdf: {gdf_NUTS_local.crs}, raster_ISA: {raster_crs}, excluder: {excluder.crs}, cutout: {c.crs}"
    )


    ##### Add ISA criterion to excluder
    # codes must go in a list, otherwise, the zero index works wrongly
    # Pass the raster path (not an already-open DatasetReader) so atlite can manage opening/closing internally without ending up with closed handles.
    excluder.add_raster(file_raster_ISA, codes=ISA_list, invert=True, crs=raster_crs) # crs is the input raster crs, not the output excluder crs, because atlite will handle the reprojection internally when building the mask, and it needs to know the input raster CRS to do it correctly. 
    
    ### Define shape from the region geometry in cutout CRS   
    _log_and_print(
        f"[get_CAPACITY_matrix] Getting shape to use in availability matrix. Reprojecting gdf_NUTS_local from {gdf_NUTS_local.crs} to {c.crs}."
    ) 
    shape = gdf_NUTS_local.to_crs(c.crs).geometry

    ##### Add CF threshold criterion when computing A matrix
    A = c.availabilitymatrix(shape, excluder).where(CF >= CF_threshold, 0)
    

    ##### Generate CAPACITY matrices 
    # Generatearea matrix, AREA, in km2
    _log_and_print(
        f"[get_CAPACITY_matrix] Getting grid from cutout to use in AREA calculation. Reprojecting cutout from {c.crs} to fixed:3035."
    ) 
    df = c.grid.to_crs(3035)
    df['area'] = df.area*1e-6
    AREA = xr.Dataset.from_dataframe(df.set_index(["y","x"]))["area"]


    ##### Compute the CAPACITY matrix
    CAPACITY = AREA * cap_per_sqkm * A


    return CAPACITY.round(6)



############################## Unwrap relevant variables

##### params
cutout_params = snakemake.params["cutout_params"]
cap_per_sqkm = snakemake.params["cap_per_sqkm"]
CF_threshold = snakemake.params["CF_threshold"]
isa = snakemake.params["isa"]
if isa is None:
    ISA_list = [0, 1, 2, 3, 4]
else:
    ISA_list = [int(isa)]
##### input
file_gdf_NUTS = snakemake.input["gdf_NUTS"]
file_nc_CF = snakemake.input["nc_CF"]
file_raster_ISA = snakemake.input["raster_ISA"]
##### output
file_CAPACITY = snakemake.output["nc_CAPACITY"]
##### wildcards
cutout = snakemake.wildcards["cutout"]
year =snakemake.wildcards["year"]
region = snakemake.wildcards["region"]
resource = snakemake.wildcards["resource"]



############################## Operations

##### Load gdf_NUTS_local
_, gdf_NUTS_local = load_gdf_nuts(file_gdf_NUTS, region)


##### Load and limit cutout
file_cutout = resolve_user_home_path(cutout_params[f"{cutout}_{year}"]["path"])
c = load_and_limit_cutout(file_cutout, gdf_NUTS_local)


##### Load CF
CF = load_CF(file_nc_CF)



############################## Create outputs
CAPACITY = get_CAPACITY_matrix(
    file_raster_ISA,
    ISA_list,
    gdf_NUTS_local,
    c,
    CF,
    CF_threshold,
    cap_per_sqkm,
)

CAPACITY.to_netcdf(file_CAPACITY)







    







