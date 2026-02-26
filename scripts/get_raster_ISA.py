import geopandas as gpd
import json

import rasterio
from rasterio.mask import mask

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
    raster_crs = raster_ISA.crs

    ##### Load gdf_NUTS and apply operations (reproject to raster CRS)
    gdf_NUTS = (gpd.read_file(file_gdf_NUTS)
                .set_index("NUTS_ID")            # set index 
                .loc[[region]]                   # filter region
                .to_crs(raster_crs)              # change crs to that of the ISA raster
    )

    ##### Filter raster_ISA with vectorial mask from region
    # Put geometry in GeoJSON-like dict
    geoms = [json.loads(gdf_NUTS.to_json())["features"][0]["geometry"]]
    # Apply filter (returns array with shape (bands, height, width))
    out_image, out_transform = mask(raster_ISA, geoms, crop=True)
    # Update metadata based on source
    out_meta = raster_ISA.meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": out_image.shape[1],
        "width": out_image.shape[2],
        "transform": out_transform,
        "dtype": raster_ISA.meta.get("dtype", out_image.dtype),
        "compress": "LZW",
        "nodata": raster_ISA.nodata
    })



############################## Create outputs
with rasterio.open(file_raster_ISA, "w", **out_meta) as dest:
    dest.write(out_image)
