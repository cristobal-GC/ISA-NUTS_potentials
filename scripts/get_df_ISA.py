import pandas as pd
import numpy as np
import rasterio
from utils import log_raster_spatial_info

from typing import Any
snakemake: Any



##############################
# This script reads the ISA raster for a given region, counts the number of pixels in each class (0-4), computes the area and percentage for each class, and saves the results in a dataframe as a CSV file. The script also adds a TOTAL row with the sum of counts, area, and percentage.
#
# The area is computed using the pixel resolution from the raster metadata


############################## Unwrap relevant variables

file_raster_ISA = snakemake.input["raster_ISA"]
file_df_ISA = snakemake.output["df_ISA"]

region = snakemake.wildcards["region"]
resource = snakemake.wildcards["resource"]



############################## Operations

with rasterio.open(file_raster_ISA) as raster_ISA:

    log_raster_spatial_info(raster_ISA, source_label=file_raster_ISA)

    band = raster_ISA.read(1)

    ##### Pixel area (km²)
    pixel_width, pixel_height = raster_ISA.res
    pixel_area_km2 = abs(pixel_width * pixel_height) * 1e-6

    ##### Remove nodata
    nodata = raster_ISA.nodata
    if nodata is not None:
        band = band[band != nodata] # 65535



##### Count pixels per class (fast)

# flatten array for bincount
band_flat = band.ravel()

# ensure at least classes 0–4 exist
# bitcount is more efficienty than np.unique for counting occurrences of integer values, but it requires the input to be non-negative integers and will count all integers up to the maximum value in the input. By using minlength=5, we ensure that we get counts for classes 0–4 even if some of them are not present in the data, and we ignore any values above 4.
counts = np.bincount(band_flat, minlength=5)

# keep only 0-4
counts = counts[:5]



##### Build dataframe

df = pd.DataFrame({
    "value": np.arange(5),
    "counts": counts
})


##### Compute areas and percentages

df["area"] = df["counts"] * pixel_area_km2
df["porc"] = 100 * df["counts"] / df["counts"].sum()



##### Add TOTAL row

totals = df.drop(columns=["value"]).sum(numeric_only=True)

row_total = pd.DataFrame(
    {**{"value": "TOTAL"}, **totals.to_dict()},
    index=[0]
)

df = pd.concat([df, row_total], ignore_index=True)

df = df.round({"area": 6, "porc": 6})

df.set_index("value", inplace=True)



############################## Create outputs

df.to_csv(file_df_ISA)