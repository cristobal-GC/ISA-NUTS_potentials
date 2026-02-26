
#################### Auxiliary functions

RESOURCE_RASTERS = {
    "onwind": "data/ISA/Clas_ISA_eol_pb.tiff",
    "solar": "data/ISA/Clas_ISA_ftv_pb.tiff",
}

def get_file_ISA(wc):
    try:
        return RESOURCE_RASTERS[wc.resource]
    except KeyError:
        raise ValueError(f"Invalid resource: {wc.resource}")



#################### get_raster_ISA
#
# This rule is to generate an ISA raster for a specific region from the Spanis ISA raster
#
# Wildcards:
#   - region    [ES11, ... ]
#   - resource  [onwind, solar]

rule get_raster_ISA:
    message:
        "... Getting raster_ISA for resource: {wildcards.resource} and region: {wildcards.region}."
    input:
        gdf_NUTS ="data/NUTS/NUTS_RG_01M_2021_4326_ES.geojson",
        raster_ISA_miteco=get_file_ISA
    output:
        raster_ISA="results/rasters/ISA/raster_ISA_{resource}_{region}.tiff"
    script:
        "../scripts/get_raster_ISA.py"



#################### get_df_ISA
#
# This rule is to generate a df with the info about surface percentages for ISA classes from raster_ISA
#
# Wildcards:
#   - region    [ES11, ... ]
#   - resource  [onwind, solar]

rule get_df_ISA:
    message:
        "... Getting df_ISA for resource: {wildcards.resource} and region: {wildcards.region}."
    params:
        cutout_params=config["cutout_params"],
    input:        
        raster_ISA="results/rasters/ISA/raster_ISA_{resource}_{region}.tiff"
    output:
        df_ISA="results/dfs/ISA/df_ISA_{resource}_{region}.csv"
    script:
        "../scripts/get_df_ISA.py"
    


#################### get_nc_CF
#
# This rule is to generate a file with the CF for a specific resource, region and cutout
#
# Wildcards:
#   - region    [ES11, ... ]
#   - resource  [onwind, solar]

rule get_nc_CF:
    message:
        "... Getting nc_CF for cutout: {wildcards.cutout}, year: {wildcards.year}, resource: {wildcards.resource} and region: {wildcards.region}."
    params:
        cutout_params=config["cutout_params"],
        CF_params=config["CF_params"]
    input:
        gdf_NUTS ="data/NUTS/NUTS_RG_01M_2021_4326_ES.geojson",
    output:
        nc_CF="results/ncs/CF/CF_{resource}_{region}_{cutout}_{year}.nc"
    script:
        "../scripts/get_nc_CF.py"