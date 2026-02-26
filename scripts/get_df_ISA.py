import pandas as pd
import numpy as np
import rasterio

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

##### Load raster_ISA
raster_ISA = rasterio.open(file_raster_ISA)

##### Get ISA band
band = raster_ISA.read(1)

##### Get areas and percentages of each ISA code
# Unique values and counts
unique, counts = np.unique(band, return_counts=True)
df = pd.DataFrame({'value': unique, 'counts': counts})
# Remove 65535
df = df[df['value'] != 65535]
# Make it sure all levels 0-4 exist
all_values = pd.DataFrame({'value': np.arange(5)})
df = all_values.merge(df, on='value', how='left').fillna(0)
# Assign area
df['area'] = df['counts']*0.025*0.025
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
df = df.round({'area': 2, 'porc': 2})

##### Set index: 'value'
df.set_index('value', inplace=True)



############################## Create outputs
df.to_csv(file_df_ISA)



