import xarray as xr
import logging
import rasterio

from atlite.gis import ExclusionContainer
from utils import (
    load_and_limit_cutout,
    load_gdf_nuts_local,
    resolve_user_home_path,
    log_raster_spatial_info,
    log_xarray_spatial_info,
)

from typing import Any
snakemake: Any  # This is to avoid my IDE to complain about snakemake variable not being defined, but it is actually defined when running the script with snakemake


logger = logging.getLogger(__name__)


def _log_and_print(message):
    logger.info(message)
    print(message)


def validate_and_align_CF_coordinates(CF, cutout):
    """
    Validate CF coordinates match cutout grid and reindex if necessary.
    
    This is important for filtering the CF with the ISA excluder, which is based 
    on the cutout grid. If they don't match, the filtering will not work and the 
    resulting A_CF and A_CF_ISA will be wrong.
    
    Parameters
    ----------
    CF : xr.DataArray
        Capacity factor data array
    cutout : atlite.Cutout
        Cutout object with grid coordinates
        
    Returns
    -------
    xr.DataArray
        CF data array aligned to cutout grid
    """
    _log_and_print("[validate_and_align_CF_coordinates] Validating coordinate alignment:")
    _log_and_print(f"  Cutout x: [{cutout.data.x.min().values:.2f}, {cutout.data.x.max().values:.2f}], shape: {len(cutout.data.x)}")
    _log_and_print(f"  Cutout y: [{cutout.data.y.min().values:.2f}, {cutout.data.y.max().values:.2f}], shape: {len(cutout.data.y)}")
    _log_and_print(f"  CF x: [{CF.x.min().values:.2f}, {CF.x.max().values:.2f}], shape: {len(CF.x)}")
    _log_and_print(f"  CF y: [{CF.y.min().values:.2f}, {CF.y.max().values:.2f}], shape: {len(CF.y)}")
    
    coords_match = CF.x.equals(cutout.data.x) and CF.y.equals(cutout.data.y)
    
    if not coords_match:
        _log_and_print("Coordinates mismatch detected. Reindexing CF to cutout grid.")
        CF = CF.sel(x=cutout.data.x, y=cutout.data.y, method="nearest")
    else:
        _log_and_print("Coordinates match perfectly.")
    
    return CF



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
    
    ##### Create excluder container, with 25x25 m resolution
    excluder = ExclusionContainer(res=25) # no need to specify CRS here, by default it is 3035, and atlite will handle reprojection internally when adding the raster criterion, as long as we provide the correct input raster CRS.

    with rasterio.open(file_raster_ISA) as raster_ISA:
        log_raster_spatial_info(raster_ISA, source_label=file_raster_ISA)
        raster_crs = raster_ISA.crs

    if raster_crs is None:
        raise ValueError(f"Raster {file_raster_ISA} has no CRS. Cannot build ExclusionContainer mask reliably.")

    _log_and_print(
        f"[get_CAPACITY_matrix] CRS chain -> gdf: {gdf_NUTS_local.crs}, raster_ISA: {raster_crs}, excluder: {excluder.crs}, cutout: {c.crs}"
    )

    ##### Add ISA criterion to excluder
    # codes must go in a list, otherwise, the zero index works wrongly
    # Pass the raster path (not an already-open DatasetReader) so atlite can manage opening/closing internally without ending up with closed handles.
    excluder.add_raster(file_raster_ISA, codes=ISA_list, invert=True, crs=raster_crs) # crs is the input raster crs, not the output excluder crs, because atlite will handle the reprojection internally when building the mask, and it needs to know the input raster CRS to do it correctly. If we pass the excluder CRS (which is 3035), atlite will assume the input raster is already in 3035, which is not the case, and the resulting mask will be wrong.
    
    ### Define shape from the region geometry in cutout CRS
    shape = gdf_NUTS_local.geometry

    ##### Add CF threshold criterion when computing A matrix
    A = c.availabilitymatrix(shape, excluder).where(CF >= CF_threshold, 0)
    

    ##### Generate CAPACITY matrices 
    # Generatearea matrix, AREA, in km2
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
gdf_NUTS_local = load_gdf_nuts_local(file_gdf_NUTS, region)


##### Load and limit cutout
file_cutout = resolve_user_home_path(cutout_params[f"{cutout}_{year}"]["path"])
c = load_and_limit_cutout(file_cutout, gdf_NUTS_local)


##### Load CF
CF = xr.open_dataarray(file_nc_CF)
log_xarray_spatial_info(CF, source_label=file_nc_CF)
# Validate CF coordinates match cutout grid
CF = validate_and_align_CF_coordinates(CF, c)



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







    







