
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
        raster_ISA="results/rasters/ISA/raster_ISA_{resource}_{region}.tiff",
    output:
        df_ISA="results/dfs/ISA/df_ISA_{resource}_{region}.csv",
    script:
        "../scripts/get_df_ISA.py"
    


#################### get_nc_CF
#
# This rule is to generate a file with the CF for a specific resource, region and cutout
#
# Wildcards:
#   - region    [ES11, ... ]
#   - resource  [onwind, solar]
#   - cutout    [era5, ...]
#   - year      [2013, ...]

rule get_nc_CF:
    message:
        "... Getting nc_CF for cutout: {wildcards.cutout}, year: {wildcards.year}, resource: {wildcards.resource} and region: {wildcards.region}."
    params:
        cutout_params=config["cutout_params"],
        CF_params=config["CF_params"],
    input:
        gdf_NUTS="data/NUTS/NUTS_RG_01M_2021_4326_ES.geojson",
    output:
        nc_CF="results/ncs/CF/{cutout}/CF_{resource}_{region}_{year}.nc",
    script:
        "../scripts/get_nc_CF.py"



#################### get_nc_CAPACITY_CF_ISA , get_nc_CAPACITY_CF , get_nc_CAPACITY_ISA
#
# This rule is to get the CAPACITY matrix with one or several filters: CF threshold , ISA code
#
# Wildcard 'filters' is constrained to three possibilities. Then, params are conditioned to this filter value. This can be scaled up with more filters
#
# Wildcards:
#   - region    [ES11, ... ]
#   - resource  [onwind, solar]
#   - cutout    [era5, ...]
#   - year      [2013, ...]
#   - isa       [0, 1, ...]


rule get_nc_CAPACITY:
    wildcard_constraints:
        resource="onwind|solar",
        filters="CF|ISA\d+|CF_ISA\d+",

    message:
        "... Getting CAPACITY matrix for cutout: {wildcards.cutout}, year: {wildcards.year}, resource: {wildcards.resource}, region: {wildcards.region}, filters: {wildcards.filters}."

    params:
        cutout_params=config["cutout_params"],

        cap_per_sqkm=lambda w: config["CF_params"][w.resource]["cap_per_sqkm"],

        CF_threshold=lambda w: (
            config["CF_params"][w.resource]["CF_threshold"]
            if "CF" in w.filters
            else 0
        ),

        isa=lambda w: (
            int(w.filters.split("ISA")[1])
            if "ISA" in w.filters
            else None
        )

    input:
        gdf_NUTS="data/NUTS/NUTS_RG_01M_2021_4326_ES.geojson",
        nc_CF="results/ncs/CF/{cutout}/CF_{resource}_{region}_{year}.nc",
        raster_ISA="results/rasters/ISA/raster_ISA_{resource}_{region}.tiff"

    output:
        nc_CAPACITY="results/ncs/CAPACITY/{cutout}/CAPACITY_{filters}_{resource}_{region}_{year}.nc"

    script:
        "../scripts/get_nc_CAPACITY.py"


#################### get_df_CF_CAPACITY
#
# This rule generates a CF-CAPACITY dataframe for a specific ISA code.
#
# Wildcards:
#   - region    [ES11, ... ]
#   - resource  [onwind, solar]
#   - year      [2013, ...]
#   - isa       [0, 1, ...]

rule get_df_CF_CAPACITY:
    message:
        "... Getting df_CF_CAPACITY for cutout: {wildcards.cutout}, year: {wildcards.year}, resource: {wildcards.resource}, region: {wildcards.region}, ISA: {wildcards.isa}."
    input:
        nc_CF="results/ncs/CF/{cutout}/CF_{resource}_{region}_{year}.nc",
        nc_CAPACITY="results/ncs/CAPACITY/{cutout}/CAPACITY_ISA{isa}_{resource}_{region}_{year}.nc",
    output:
        df_CF_CAPACITY="results/dfs/CF_CAPACITY/{cutout}/df_CF_CAPACITY_ISA{isa}_{resource}_{region}_{year}.csv"
    script:
        "../scripts/get_df_CF_CAPACITY.py"





