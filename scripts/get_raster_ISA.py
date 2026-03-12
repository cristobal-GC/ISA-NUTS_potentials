import json

import rasterio
from rasterio.mask import mask
from utils import load_gdf_nuts_local, log_raster_spatial_info

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
with rasterio.open(file_raster_ISA_miteco) as raster_ISA:
    log_raster_spatial_info(raster_ISA, source_label=file_raster_ISA_miteco)

    if raster_ISA.crs is None:
        raise ValueError(f"Input raster has no CRS: {file_raster_ISA_miteco}")

    ##### Load gdf_NUTS local and reproject to source raster CRS
    gdf_NUTS_local = load_gdf_nuts_local(file_gdf_NUTS, region).to_crs(raster_ISA.crs)

    print(f"[get_raster_ISA] Clipping in source raster CRS: {raster_ISA.crs}.")

    geoms = [json.loads(gdf_NUTS_local.to_json())["features"][0]["geometry"]]
    out_image, out_transform = mask(raster_ISA, geoms, crop=True)

    out_meta = raster_ISA.meta.copy()
    out_meta.update(
        {
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform,
            "dtype": raster_ISA.meta.get("dtype", out_image.dtype),
            "compress": "LZW",
            "nodata": raster_ISA.nodata,
        }
    )



############################## Create outputs
with rasterio.open(file_raster_ISA, "w", **out_meta) as dest:
    dest.write(out_image)

with rasterio.open(file_raster_ISA) as raster_ISA_out:
    log_raster_spatial_info(raster_ISA_out, source_label=file_raster_ISA)
