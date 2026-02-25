import geopandas as gpd
import json

import rasterio
from rasterio.mask import mask



############################## Unwrap relevant variables

##### input
file_gdf_NUTS = snakemake.input["gdf_NUTS"]
file_raster_ISA = snakemake.input["raster_ISA"]
##### output
file_raster_ISA_local = snakemake.output["raster_ISA_local"]
##### wildcards
region = snakemake.wildcards["region"]



############################## Operations

##### Load raster_ISA
raster_ISA = rasterio.open(file_raster_ISA)

##### Load gdf_NUTS and apply operations
gdf_NUTS = (gpd.read_file(file_gdf_NUTS)
           .set_index("NUTS_ID")            # set index 
           .loc[[region]]                   # filter region
           .to_crs(raster_ISA.crs)          # change crs to that of the ISA raster
)

##### Filter raster_ISA with vectorial mask from region
# Put geometry in GeoJSON-like dict
geoms = [json.loads(gdf_NUTS.to_json())["features"][0]["geometry"]]
# Apply filter
out_image, out_transform = mask(raster_ISA, geoms, crop=True)
# Update metadata
out_meta = raster_ISA.meta.copy()
out_meta.update({
    "driver": "GTiff",
    "height": out_image.shape[1],
    "width": out_image.shape[2],
    "transform": out_transform,
    "dtype": raster_ISA.meta["dtype"],
    "compress": "LZW",
    "nodata": raster_ISA.nodata
})



############################## Create outputs
with rasterio.open(file_raster_ISA_local, "w", **out_meta) as dest:
    dest.write(out_image)
