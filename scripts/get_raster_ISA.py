import json

import rasterio
from rasterio.io import MemoryFile
from rasterio.mask import mask
from rasterio.warp import calculate_default_transform, reproject, Resampling
from utils import load_gdf_nuts_local

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

target_crs = "EPSG:4326"

##### Load raster_ISA and apply vectorial mask inside a context manager
with rasterio.open(file_raster_ISA_miteco) as raster_ISA:
    ##### Reproject source raster to target CRS (EPSG:4326)
    dst_transform, dst_width, dst_height = calculate_default_transform(
        raster_ISA.crs,
        target_crs,
        raster_ISA.width,
        raster_ISA.height,
        *raster_ISA.bounds,
    )

    reproj_meta = raster_ISA.meta.copy()
    reproj_meta.update(
        {
            "driver": "GTiff",
            "crs": target_crs,
            "transform": dst_transform,
            "width": dst_width,
            "height": dst_height,
            "compress": "LZW",
            "nodata": raster_ISA.nodata,
        }
    )

    ##### Load gdf_NUTS local and reproject to target CRS
    gdf_NUTS_local = load_gdf_nuts_local(file_gdf_NUTS, region).to_crs(target_crs)

    ##### Reproject in memory (nearest for categorical ISA), then mask by region geometry
    with MemoryFile() as memfile:
        with memfile.open(**reproj_meta) as raster_ISA_4326:
            for band_idx in range(1, raster_ISA.count + 1):
                reproject(
                    source=rasterio.band(raster_ISA, band_idx),
                    destination=rasterio.band(raster_ISA_4326, band_idx),
                    src_transform=raster_ISA.transform,
                    src_crs=raster_ISA.crs,
                    dst_transform=dst_transform,
                    dst_crs=target_crs,
                    src_nodata=raster_ISA.nodata,
                    dst_nodata=raster_ISA.nodata,
                    resampling=Resampling.nearest,
                )

            geoms = [json.loads(gdf_NUTS_local.to_json())["features"][0]["geometry"]]
            out_image, out_transform = mask(raster_ISA_4326, geoms, crop=True)

            out_meta = raster_ISA_4326.meta.copy()
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
