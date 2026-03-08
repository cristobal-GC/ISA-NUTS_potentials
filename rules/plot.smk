

rule plot_ISA:
    message:
        "... Plotting ISA map for resource: {wildcards.resource}, region: {wildcards.region}, resolution: {wildcards.resolution} format: {wildcards.format}."
    params:
        fig_params=config["fig_params"]
    input:
        gdf_NUTS ="data/NUTS/NUTS_RG_01M_2021_4326_ES.geojson",
        raster_ISA="results/rasters/ISA/raster_ISA_{resource}_{region}.tiff",
        df_ISA="results/dfs/ISA/df_ISA_{resource}_{region}.csv"
    output:
        map_ISA="results/maps/ISA/{resolution}/ISA_{resource}_{region}_{resolution}.{format}"
    script:
        "../scripts/plot_raster_ISA.py"


rule plot_cutout:
    message:
        "... Plotting cutout map for cutout: {wildcards.cutout}, year: {wildcards.year}, resource: {wildcards.resource}, region: {wildcards.region}, format: {wildcards.format}."
    params:
        cutout_params=config["cutout_params"],
        fig_params=config["fig_params"]
    input:
        gdf_NUTS ="data/NUTS/NUTS_RG_01M_2021_4326_ES.geojson",
    output:
        map_cutout="results/maps/cutout/{cutout}/cutout_{resource}_{region}_{year}.{format}"
    script:
        "../scripts/plot_cutout.py"        


rule plot_CF:
    message:
        "... Plotting CF map for cutout: {wildcards.cutout}, year: {wildcards.year}, resource: {wildcards.resource}, region: {wildcards.region}, format: {wildcards.format}."
    params:
        cutout_params=config["cutout_params"],
        fig_params=config["fig_params"]
    input:
        gdf_NUTS ="data/NUTS/NUTS_RG_01M_2021_4326_ES.geojson",
        nc_CF="results/ncs/CF/{cutout}/CF_{resource}_{region}_{year}.nc"
    output:
        map_CF="results/maps/CF/{cutout}/CF_{resource}_{region}_{year}.{format}"
    script:
        "../scripts/plot_nc_CF.py"


rule plot_CAPACITY_CF_ISA:
    wildcard_constraints:
        resource="onwind|solar",
        filters="CF|ISA\d+|CF_ISA\d+",
    message:
        "... Plotting CAPACITY map for cutout: {wildcards.cutout}, year: {wildcards.year}, resource: {wildcards.resource}, region: {wildcards.region}, filters: {wildcards.filters}, format: {wildcards.format}."
    params:
        fig_params=config["fig_params"]
    input:
        gdf_NUTS ="data/NUTS/NUTS_RG_01M_2021_4326_ES.geojson",
        nc_CAPACITY="results/ncs/CAPACITY/{cutout}/CAPACITY_{filters}_{resource}_{region}_{year}.nc"
    output:
        map_CAPACITY="results/maps/CAPACITY/{cutout}/CAPACITY_{filters}_{resource}_{region}_{year}.{format}"
    script:
        "../scripts/plot_nc_CAPACITY.py"


rule plot_df_CF_CAPACITY:
    message:
        "... Plotting CF vs CUM_CAPACITY curves for cutout: {wildcards.cutout}, year: {wildcards.year}, resource: {wildcards.resource}, region: {wildcards.region}, format: {wildcards.format}."
    params:
        fig_params=config["fig_params"]
    input:
        dfs_CF_CAPACITY=expand(
            "results/dfs/CF_CAPACITY/{{cutout}}/df_CF_CAPACITY_ISA{isa}_{{resource}}_{{region}}_{{year}}.csv",
            isa=[0, 1, 2, 3, 4]
        )
    output:
        plot_CF_CAPACITY="results/figs/CF_CAPACITY/{cutout}/CF_CAPACITY_{resource}_{region}_{year}.{format}"
    script:
        "../scripts/plot_df_CF_CAPACITY.py"