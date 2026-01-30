

import numpy as np
import geopandas as gpd
import yaml
import json

import rasterio
from rasterio.mask import mask



############################## Unwrap relevant variables

##### input
file_NUTS       = snakemake.input["file_NUTS"]
file_ISA        = snakemake.input["file_ISA"]
##### output
file_ISA_local  = snakemake.output["file_ISA_local"]
##### wildcards
region          = snakemake.wildcards["region"]
resource        = snakemake.wildcards["resource"]



############################## Operations

##### Load gdf with nuts
gdf_NUTS = (gpd.read_file(file_NUTS)
               .set_index("id")
               .loc[[region]]
)


input(f'gdf_NUTS is {gdf_NUTS}')


############################## Create outputs
file_ISA_local.touch()