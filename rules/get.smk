
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


def get_regions_for_nuts(wc):
    regions_cfg = config["regions"]
    if not isinstance(regions_cfg, dict):
        raise ValueError("config['regions'] must be a dictionary with NUTS keys for get_df_summary")
    if wc.nuts not in regions_cfg:
        raise ValueError(f"NUTS level {wc.nuts} not found in config['regions']")
    return regions_cfg[wc.nuts] or []


# Geometry source for a given level. Standard NUTS levels share the official
# NUTS GeoJSON; the custom CIMAS level uses a GeoJSON built from the rectangular
# domains declared in config['CIMAS_domains'] (rule get_geojson_CIMAS); the
# custom ADM3 level uses the geoBoundaries municipalities GeoJSON (identified by
# 'shapeName', without the NUTS columns).
# These are defined in the first included rules file so they are available to all
# rules (get/plot/latex).
NUTS_GEOJSON = "data/NUTS/NUTS_RG_01M_2021_4326_ES.geojson"
CIMAS_GEOJSON = "data/NUTS/CIMAS.geojson"
ADM3_GEOJSON = "data/NUTS/geoBoundaries-ESP-ADM3_simplified.geojson"


def gdf_nuts_input(wc):
    if wc.nuts == "CIMAS":
        return CIMAS_GEOJSON
    if wc.nuts == "ADM3":
        return ADM3_GEOJSON
    return NUTS_GEOJSON


# Source for the thin grey context boundaries drawn behind a region:
#   - NUTS0/2/3 and CIMAS: the official NUTS GeoJSON (same-level regions, or
#     NUTS3 provinces for CIMAS).
#   - ADM3: the municipalities GeoJSON itself (neighbouring municipalities).
# The level filtering happens in utils.load_context_boundaries.
def gdf_context_input(wc):
    if wc.nuts == "ADM3":
        return ADM3_GEOJSON
    return NUTS_GEOJSON



#################### get_geojson_CIMAS
#
# This rule builds data/NUTS/CIMAS.geojson from config['CIMAS_domains']. Each
# domain becomes one feature with a rectangular Polygon derived from its
# lat_min/lat_max/lon_min/lon_max bounds, in EPSG:4326, following the same
# schema (NUTS_ID, LEVL_CODE, geometry, ...) as the official NUTS GeoJSON so the
# rest of the workflow can consume it transparently.

rule get_geojson_CIMAS:
    message:
        "... [get_geojson_CIMAS] Building CIMAS GeoJSON from config['CIMAS_domains']."
    params:
        cimas_domains=config.get("CIMAS_domains", {}),
    input:
        config="config/config.yaml",
    output:
        geojson=CIMAS_GEOJSON,
    script:
        "../scripts/get_geojson_CIMAS.py"



#################### get_raster_ISA
#
# This rule is to generate an ISA raster for a specific region from the Spanish ISA raster
#
# Wildcards:
#   - nuts      [NUTS2, NUTS3]
#   - region    [ES11, ... ]
#   - resource  [onwind, solar]

rule get_raster_ISA:
    wildcard_constraints:
        nuts="NUTS0|NUTS2|NUTS3|CIMAS|ADM3",
    threads:
        lambda w: get_rule_threads(
            "get_raster_ISA",
            f"benchmarks/get_raster_ISA/{w.nuts}/get_raster_ISA_{w.resource}_{w.region}.tsv",
        )
    resources:
        mem_mb=lambda w: get_rule_mem_mb(
            "get_raster_ISA",
            f"benchmarks/get_raster_ISA/{w.nuts}/get_raster_ISA_{w.resource}_{w.region}.tsv",
        )
    benchmark:
        "benchmarks/get_raster_ISA/{nuts}/get_raster_ISA_{resource}_{region}.tsv"
    message:
        "... [get_raster_ISA] Getting raster_ISA for resource: {wildcards.resource} and region: {wildcards.region}."
    input:
        gdf_NUTS=gdf_nuts_input,
        raster_ISA_miteco=get_file_ISA
    output:
        raster_ISA="results/rasters/ISA/{nuts}/raster_ISA_{resource}_{region}.tiff"
    script:
        "../scripts/get_raster_ISA.py"



#################### get_df_ISA
#
# This rule is to generate a df with the info about surface percentages for ISA classes from raster_ISA
#
# Wildcards:
#   - nuts      [NUTS2, NUTS3]
#   - region    [ES11, ... ]
#   - resource  [onwind, solar]

rule get_df_ISA:
    wildcard_constraints:
        nuts="NUTS0|NUTS2|NUTS3|CIMAS|ADM3",
    threads:
        lambda w: get_rule_threads(
            "get_df_ISA",
            f"benchmarks/get_df_ISA/{w.nuts}/get_df_ISA_{w.resource}_{w.region}.tsv",
        )
    resources:
        mem_mb=lambda w: get_rule_mem_mb(
            "get_df_ISA",
            f"benchmarks/get_df_ISA/{w.nuts}/get_df_ISA_{w.resource}_{w.region}.tsv",
        )
    benchmark:
        "benchmarks/get_df_ISA/{nuts}/get_df_ISA_{resource}_{region}.tsv"
    message:
        "... [get_df_ISA] Getting df_ISA for resource: {wildcards.resource} and region: {wildcards.region}."
    params:
        cutout_params=config["cutout_params"],
    input:        
        raster_ISA="results/rasters/ISA/{nuts}/raster_ISA_{resource}_{region}.tiff",
    output:
        df_ISA="results/dfs/ISA/{nuts}/df_ISA_{resource}_{region}.csv",
    script:
        "../scripts/get_df_ISA.py"
    


#################### get_nc_CF
#
# This rule is to generate a file with the CF for a specific resource, region and cutout
#
# Wildcards:
#   - cutout    [era5, ...]
#   - nuts      [NUTS2, NUTS3]
#   - region    [ES11, ... ]
#   - resource  [onwind, solar]
#   - year      [2013, ...]

rule get_nc_CF:
    wildcard_constraints:
        nuts="NUTS0|NUTS2|NUTS3|CIMAS|ADM3",
    threads:
        lambda w: get_rule_threads(
            "get_nc_CF",
            f"benchmarks/get_nc_CF/{w.cutout}/{w.nuts}/get_nc_CF_{w.resource}_{w.region}_{w.year}.tsv",
        )
    resources:
        mem_mb=lambda w: get_rule_mem_mb(
            "get_nc_CF",
            f"benchmarks/get_nc_CF/{w.cutout}/{w.nuts}/get_nc_CF_{w.resource}_{w.region}_{w.year}.tsv",
        )
    benchmark:
        "benchmarks/get_nc_CF/{cutout}/{nuts}/get_nc_CF_{resource}_{region}_{year}.tsv"
    message:
        "... [get_nc_CF] Getting nc_CF for cutout: {wildcards.cutout}, year: {wildcards.year}, resource: {wildcards.resource} and region: {wildcards.region}."
    params:
        cutout_params=config["cutout_params"],
        CF_params=config["CF_params"],
    input:
        gdf_NUTS=gdf_nuts_input,
    output:
        nc_CF="results/ncs/CF/{cutout}/{nuts}/CF_{resource}_{region}_{year}.nc",
    script:
        "../scripts/get_nc_CF.py"



#################### get_nc_CAPACITY_CFth_ISA , get_nc_CAPACITY_CFth , get_nc_CAPACITY_ISA
#
# This rule is to get the CAPACITY matrix with one or several filters: CF threshold , ISA code
#
# Wildcard 'filters' is constrained to three possibilities. Then, params are conditioned to this filter value. This can be scaled up with more filters
#
# Wildcards:
#   - cutout    [era5, ...]
#   - nuts      [NUTS2, NUTS3]
#   - filters   [CFth, ISA0, ... , CFth_ISA0, ...]
#   - region    [ES11, ... ]
#   - resource  [onwind, solar]
#   - year      [2013, ...]


rule get_nc_CAPACITY:
    wildcard_constraints:
        nuts="NUTS0|NUTS2|NUTS3|CIMAS|ADM3",
        resource="onwind|solar",
        filters=r"CFth|ISA\d+|CFth_ISA\d+",
    threads:
        lambda w: get_rule_threads(
            "get_nc_CAPACITY",
            f"benchmarks/get_nc_CAPACITY/{w.cutout}/{w.nuts}/get_nc_CAPACITY_{w.filters}_{w.resource}_{w.region}_{w.year}.tsv",
        )
    resources:
        mem_mb=lambda w: get_rule_mem_mb(
            "get_nc_CAPACITY",
            f"benchmarks/get_nc_CAPACITY/{w.cutout}/{w.nuts}/get_nc_CAPACITY_{w.filters}_{w.resource}_{w.region}_{w.year}.tsv",
        )

    benchmark:
        "benchmarks/get_nc_CAPACITY/{cutout}/{nuts}/get_nc_CAPACITY_{filters}_{resource}_{region}_{year}.tsv"
    message:
        "... [get_nc_CAPACITY] Getting CAPACITY matrix for cutout: {wildcards.cutout}, year: {wildcards.year}, resource: {wildcards.resource}, region: {wildcards.region}, filters: {wildcards.filters}."

    params:
        cutout_params=config["cutout_params"],

        cap_per_sqkm=lambda w: config["CF_params"][w.resource]["cap_per_sqkm"],

        CF_threshold=lambda w: (
            config["CF_params"][w.resource]["CF_threshold"]
            if "CFth" in w.filters
            else 0
        ),

        isa=lambda w: (
            int(w.filters.split("ISA")[1])
            if "ISA" in w.filters
            else None
        )

    input:
        gdf_NUTS=gdf_nuts_input,
        nc_CF="results/ncs/CF/{cutout}/{nuts}/CF_{resource}_{region}_{year}.nc",
        raster_ISA="results/rasters/ISA/{nuts}/raster_ISA_{resource}_{region}.tiff"

    output:
        nc_CAPACITY="results/ncs/CAPACITY/{cutout}/{nuts}/CAPACITY_{filters}_{resource}_{region}_{year}.nc"

    script:
        "../scripts/get_nc_CAPACITY.py"


#################### get_df_CF_CAPACITY
#
# This rule generates a CF-CAPACITY dataframe for a specific ISA code.
#
# Wildcards:
#   - cutout    [era5, ...]
#   - nuts      [NUTS2, NUTS3]
#   - isa       [0, 1, ...]
#   - region    [ES11, ... ]
#   - resource  [onwind, solar]
#   - year      [2013, ...]

rule get_df_CF_CAPACITY:
    wildcard_constraints:
        nuts="NUTS0|NUTS2|NUTS3|CIMAS|ADM3",
    threads:
        lambda w: get_rule_threads(
            "get_df_CF_CAPACITY",
            f"benchmarks/get_df_CF_CAPACITY/{w.cutout}/{w.nuts}/get_df_CF_CAPACITY_ISA{w.isa}_{w.resource}_{w.region}_{w.year}.tsv",
        )
    resources:
        mem_mb=lambda w: get_rule_mem_mb(
            "get_df_CF_CAPACITY",
            f"benchmarks/get_df_CF_CAPACITY/{w.cutout}/{w.nuts}/get_df_CF_CAPACITY_ISA{w.isa}_{w.resource}_{w.region}_{w.year}.tsv",
        )
    benchmark:
        "benchmarks/get_df_CF_CAPACITY/{cutout}/{nuts}/get_df_CF_CAPACITY_ISA{isa}_{resource}_{region}_{year}.tsv"
    message:
        "... [get_df_CF_CAPACITY] Getting df_CF_CAPACITY for cutout: {wildcards.cutout}, year: {wildcards.year}, resource: {wildcards.resource}, region: {wildcards.region}, ISA: {wildcards.isa}."
    input:
        nc_CF="results/ncs/CF/{cutout}/{nuts}/CF_{resource}_{region}_{year}.nc",
        nc_CAPACITY="results/ncs/CAPACITY/{cutout}/{nuts}/CAPACITY_ISA{isa}_{resource}_{region}_{year}.nc",
    output:
        df_CF_CAPACITY="results/dfs/CF_CAPACITY/{cutout}/{nuts}/df_CF_CAPACITY_ISA{isa}_{resource}_{region}_{year}.csv"
    script:
        "../scripts/get_df_CF_CAPACITY.py"


#################### get_df_CAPACITY
#
# This rule generates a CAPACITY dataframe by ISA index for one region.
#
# Wildcards:
#   - cutout    [era5, ...]
#   - nuts      [NUTS2, NUTS3]
#   - region    [ES11, ... ]
#   - resource  [onwind, solar]
#   - year      [2013, ...]

rule get_df_CAPACITY:
    wildcard_constraints:
        nuts="NUTS0|NUTS2|NUTS3|CIMAS|ADM3",
    threads:
        lambda w: get_rule_threads(
            "get_df_CAPACITY",
            f"benchmarks/get_df_CAPACITY/{w.cutout}/{w.nuts}/get_df_CAPACITY_{w.resource}_{w.region}_{w.year}.tsv",
        )
    resources:
        mem_mb=lambda w: get_rule_mem_mb(
            "get_df_CAPACITY",
            f"benchmarks/get_df_CAPACITY/{w.cutout}/{w.nuts}/get_df_CAPACITY_{w.resource}_{w.region}_{w.year}.tsv",
        )
    benchmark:
        "benchmarks/get_df_CAPACITY/{cutout}/{nuts}/get_df_CAPACITY_{resource}_{region}_{year}.tsv"
    message:
        "... [get_df_CAPACITY] Getting df_CAPACITY for cutout: {wildcards.cutout}, year: {wildcards.year}, resource: {wildcards.resource}, region: {wildcards.region}."
    params:
        cap_per_sqkm=lambda w: config["CF_params"][w.resource]["cap_per_sqkm"],
    input:
        ncs_CAPACITY=lambda w: expand(
            "results/ncs/CAPACITY/{cutout}/{nuts}/CAPACITY_ISA{isa}_{resource}_{region}_{year}.nc",
            cutout=w.cutout,
            nuts=w.nuts,
            isa=[0, 1, 2, 3, 4],
            resource=w.resource,
            region=w.region,
            year=w.year,
        ),
    output:
        df_CAPACITY="results/dfs/CAPACITY/{cutout}/{nuts}/df_CAPACITY_{resource}_{region}_{year}.csv"
    script:
        "../scripts/get_df_CAPACITY.py"


#################### get_df_summary
#
# This rule generates a summary dataframe for all regions within a NUTS level.
#
# Wildcards:
#   - cutout    [era5, ...]
#   - nuts      [NUTS2, NUTS3]
#   - resource  [onwind, solar]
#   - year      [2013, ...]

rule get_df_summary:
    wildcard_constraints:
        nuts="NUTS0|NUTS2|NUTS3|CIMAS|ADM3",
    threads:
        lambda w: get_rule_threads(
            "get_df_summary",
            f"benchmarks/get_df_summary/{w.cutout}/{w.nuts}/get_df_summary_{w.resource}_{w.year}.tsv",
        )
    resources:
        mem_mb=lambda w: get_rule_mem_mb(
            "get_df_summary",
            f"benchmarks/get_df_summary/{w.cutout}/{w.nuts}/get_df_summary_{w.resource}_{w.year}.tsv",
        )
    benchmark:
        "benchmarks/get_df_summary/{cutout}/{nuts}/get_df_summary_{resource}_{year}.tsv"
    message:
        "... [get_df_summary] Getting df_summary for cutout: {wildcards.cutout}, nuts: {wildcards.nuts}, year: {wildcards.year}, resource: {wildcards.resource}."
    params:
        regions=get_regions_for_nuts,
        CF_threshold=lambda w: config["CF_params"][w.resource]["CF_threshold"],
        cap_per_sqkm=lambda w: config["CF_params"][w.resource]["cap_per_sqkm"],
    input:
        dfs_CF_CAPACITY=lambda w: expand(
            "results/dfs/CF_CAPACITY/{cutout}/{nuts}/df_CF_CAPACITY_ISA{isa}_{resource}_{region}_{year}.csv",
            cutout=w.cutout,
            nuts=w.nuts,
            isa=[0, 1, 2, 3, 4],
            resource=w.resource,
            region=get_regions_for_nuts(w),
            year=w.year,
        ),
        dfs_CAPACITY=lambda w: expand(
            "results/dfs/CAPACITY/{cutout}/{nuts}/df_CAPACITY_{resource}_{region}_{year}.csv",
            cutout=w.cutout,
            nuts=w.nuts,
            resource=w.resource,
            region=get_regions_for_nuts(w),
            year=w.year,
        )
    output:
        df_summary="results/dfs/summary/{cutout}/{nuts}/df_summary_{resource}_{year}.csv"
    script:
        "../scripts/get_df_summary.py"





