

rule plot_ISA:
    message:
        "... Plotting ISA map for resource: {wildcards.resource}, region: {wildcards.region}, resolution: {wildcards.resolution} and format: {wildcards.format}."
    params:
        map_params=config["map_params"]
    input:
        gdf_NUTS ="data/NUTS/NUTS_RG_01M_2021_4326_ES.geojson",
        raster_ISA="results/rasters/ISA/raster_ISA_{resource}_{region}.tiff",
        df_ISA="results/dfs/ISA/df_ISA_{resource}_{region}.csv"
    output:
        map_ISA="results/maps/ISA/{resolution}/{format}/ISA_{resource}_{region}_{resolution}.{format}"
    script:
        "../scripts/plot_raster_ISA.py"


rule plot_cutout:
    message:
        "... Plotting cutout map for cutout: {wildcards.cutout}, year: {wildcards.year}, resource: {wildcards.resource}, region: {wildcards.region}, and format: {wildcards.format}."
    params:
        cutout_params=config["cutout_params"],
        map_params=config["map_params"]
    input:
        gdf_NUTS ="data/NUTS/NUTS_RG_01M_2021_4326_ES.geojson",
    output:
        map_cutout="results/maps/cutout/{cutout}/{format}/cutout_{resource}_{region}_{cutout}_{year}.{format}"
    script:
        "../scripts/plot_cutout.py"        


rule plot_CF:
    message:
        "... Plotting CF map for cutout: {wildcards.cutout}, year: {wildcards.year}, resource: {wildcards.resource}, region: {wildcards.region}, and format: {wildcards.format}."
    params:
        cutout_params=config["cutout_params"],
        map_params=config["map_params"]
    input:
        gdf_NUTS ="data/NUTS/NUTS_RG_01M_2021_4326_ES.geojson",
        nc_CF="results/ncs/CF/CF_{resource}_{region}_{cutout}_{year}.nc"
    output:
        map_CF="results/maps/CF/{cutout}/{format}/CF_{resource}_{region}_{cutout}_{year}.{format}"
    script:
        "../scripts/plot_nc_CF.py"