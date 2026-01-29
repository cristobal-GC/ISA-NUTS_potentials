

import numpy as np
import geopandas as gpd
import yaml
import json

import rasterio
from rasterio.mask import mask



############################## Unwrap relevant variables
ISA_file        = snakemake.input["ISA_file"]
ISA_local_file  = snakemake.output["ISA_local_file"]
region          = snakemake.wildcards["region"]
resource        = snakemake.wildcards["resource"]



############################## Operations

##### Load gdf with nuts

file_gdf_nuts = '../data/nuts/NUTS_RG_01M_2021_4326_LEVL_2.geojson'

gdf = gpd.read_file(file_gdf_nuts)



############################## Create outputs