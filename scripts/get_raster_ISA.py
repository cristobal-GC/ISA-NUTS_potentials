import rasterio
from rasterio.mask import mask
from utils import load_gdf_nuts, log_raster_spatial_info, _log_and_print

from typing import Any
snakemake: Any  # This is to avoid my IDE to complain about snakemake variable not being defined, but it is actually defined when running the script with snakemake



############################## Unwrap relevant variables

##### input
file_gdf_NUTS = snakemake.input["gdf_NUTS"]
file_raster_ISA_miteco = snakemake.input["raster_ISA_miteco"]
##### output
file_raster_ISA = snakemake.output["raster_ISA"]
##### wildcards
region = snakemake.wildcards["region"]



############################## Operations

##### Load raster_ISA and apply vectorial mask inside a context manager
with rasterio.open(file_raster_ISA_miteco) as raster_ISA_miteco:

    # Log spatial info of the input raster
    log_raster_spatial_info(raster_ISA_miteco, source_label=file_raster_ISA_miteco)

    # Check if the raster has a CRS defined
    if raster_ISA_miteco.crs is None:
        raise ValueError(f"Input raster has no CRS: {file_raster_ISA_miteco}")

    # Load gdf_NUTS_local and reproject to source raster CRS
    _, gdf_NUTS_local = load_gdf_nuts(file_gdf_NUTS, region)
    gdf_NUTS_local = gdf_NUTS_local.to_crs(raster_ISA_miteco.crs)

    _log_and_print(
        f"[get_raster_ISA] Reprojected gdf_NUTS_local geometry to source raster CRS: {raster_ISA_miteco.crs}."
    )

    # Prepare geometry for masking
    geoms = gdf_NUTS_local.geometry
    
    # Apply mask.
    # This step:
    #   - creates a mask using the geometry
    #   - sets nodata out of the geometry
    #   - if crop=True, cuts the bounding box to polygon bounds, which reduces the output raster size and makes it more manageable for plotting.
    # As a result:
    #   - out_image is a numpy array with the masked raster.
    #   - out_transform is the new affine transform for the masked raster, which accounts for the cropping and shifting
    out_image, out_transform = mask(raster_ISA_miteco, geoms, crop=True)

    # Update metadata for the output raster from a copy of the source raster metadata, but with the new dimensions, transform, and nodata value.
    out_meta = raster_ISA_miteco.meta.copy()
    out_meta.update(
        {
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform,
            "compress": "LZW",
            "nodata": raster_ISA_miteco.nodata,
        }
    )



############################## Create outputs
with rasterio.open(file_raster_ISA, "w", **out_meta) as dest:
    dest.write(out_image)

with rasterio.open(file_raster_ISA) as raster_ISA_out:
    log_raster_spatial_info(raster_ISA_out, source_label=file_raster_ISA)
