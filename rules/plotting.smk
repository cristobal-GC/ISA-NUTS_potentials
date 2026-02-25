

rule plot_ISA_local:
    message:
        "... Plotting ISA map for resource {wildcards.resource}, region {wildcards.region}, resolution {wildcards.resolution} and format {wildcards.format}."
    params:
        map_params=config["map_params"]
    input:
        gdf_NUTS ="data/NUTS/NUTS_RG_01M_2021_4326_ES.geojson",
        raster_ISA_local="results/rasters/ISA/raster_ISA_local_{resource}_{region}.tiff",
        df_ISA_local="results/dfs/ISA/df_ISA_local_{resource}_{region}.csv"
    output:
        map_ISA_local="results/maps/ISA/{resolution}/{format}/ISA_local_{resource}_{region}_{resolution}.{format}"
    script:
        "../scripts/plot_raster_ISA_local.py"


rule plot_cutout:
    message:
        "... Plotting cutout map for cutout {wildcards.cutout}, year {wildcards.year}, resource {wildcards.resource}, region {wildcards.region}, and format {wildcards.format}."
    params:
        cutout_params=config["cutout_params"],
        map_params=config["map_params"]
    input:
        gdf_NUTS ="data/NUTS/NUTS_RG_01M_2021_4326_ES.geojson",
    output:
        map_cutout="results/maps/cutout/{cutout}/{format}/cutout_{resource}_{region}_{cutout}_{year}.{format}"
    script:
        "../scripts/plot_cutout.py"        