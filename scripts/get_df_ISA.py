import pandas as pd
import numpy as np
import rasterio
from rasterio.warp import calculate_default_transform

from typing import Any
snakemake: Any  # This is to avoid my IDE to complain about snakemake variable not being defined, but it is actually defined when running the script with snakemake



############################## Unwrap relevant variables

##### input
file_raster_ISA = snakemake.input["raster_ISA"]
##### output
file_df_ISA = snakemake.output["df_ISA"]
##### wildcards
region = snakemake.wildcards["region"]
resource = snakemake.wildcards["resource"]



############################## Operations

# NOTE: The ISA raster is in EPSG:25830 (ETRS89 / UTM zone 30N), which is NOT an equal-area projection.
# In UTM, area distortion increases with distance from the central meridian, so pixels don't have uniform area.
# To calculate accurate areas, we need to use the actual pixel area in an equal-area CRS (EPSG:3035 - LAEA Europe),
# which is the same projection used in get_nc_CAPACITY.py for consistency in area calculations.

def calculate_pixel_area_in_equal_area_crs(raster, target_crs='EPSG:3035'):
    """
    Calculate the area of a pixel in square kilometers after reprojecting to an equal-area CRS.
    
    Parameters
    ----------
    raster : rasterio.DatasetReader
        Open raster dataset
    target_crs : str
        Target equal-area CRS (default: EPSG:3035 - LAEA Europe)
    
    Returns
    -------
    float
        Pixel area in km²
    """
    # Get transform and dimensions for reprojection
    transform, width, height = calculate_default_transform(
        raster.crs, 
        target_crs,
        raster.width,
        raster.height,
        *raster.bounds
    )
    
    # Calculate pixel area from transform (in m²)
    pixel_area_m2 = abs(transform.a * transform.e)
    
    # Convert to km²
    pixel_area_km2 = pixel_area_m2 * 1e-6
    
    return pixel_area_km2


##### Load raster_ISA and process
with rasterio.open(file_raster_ISA) as raster_ISA:
    ##### Get ISA band
    band = raster_ISA.read(1)
    
    ##### Calculate pixel area in equal-area projection (km²)
    pixel_area_km2 = calculate_pixel_area_in_equal_area_crs(raster_ISA)

##### Get areas and percentages of each ISA code
# Unique values and counts
unique, counts = np.unique(band, return_counts=True)
df = pd.DataFrame({'value': unique, 'counts': counts})
# Remove 65535
df = df[df['value'] != 65535]
# Make it sure all levels 0-4 exist
all_values = pd.DataFrame({'value': np.arange(5)})
df = all_values.merge(df, on='value', how='left').fillna(0)
# Assign area using actual pixel area in equal-area projection
df['area'] = df['counts'] * pixel_area_km2
# Assign percentage
df['porc'] = 100*df['counts'].div(df['counts'].sum())

##### Add TOTAL row
# Compute sum of all numeric columns
totals = df.drop(columns=['value']).sum(numeric_only=True)
# Add row with totals
row_total = pd.DataFrame({**{'value': 'TOTAL'}, **totals.to_dict()}, index=[0])
# Add at the end
df = pd.concat([df, row_total], ignore_index=True)
# Round
df = df.round({'area': 6, 'porc': 6})

##### Set index: 'value'
df.set_index('value', inplace=True)



############################## Create outputs
df.to_csv(file_df_ISA)



