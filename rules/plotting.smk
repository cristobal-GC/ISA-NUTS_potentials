

rule plot_ISA_local:
    message:
        "... Plotting ISA map for resource {wildcards.resource}, region {wildcards.region}, resolution {wildcards.resolution} and format {wildcards.format}."
    params:
        map_params=config["map_params"]
    input:
        gdf_NUTS ="data/NUTS/NUTS_RG_01M_2021_4326_ES.geojson",
        raster_ISA_local="results/rasters/ISA/ISA_local_{resource}_{region}.tiff"
    output:
        map_ISA_local="results/maps/ISA/{resolution}/{format}/ISA_local_{resource}_{region}_{resolution}.{format}"
    script:
        "../scripts/plot_ISA_local.py"