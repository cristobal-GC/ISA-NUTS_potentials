


#################### plot_ISA
#
# Wildcards:
#   - nuts       [NUTS2, NUTS3]
#   - region     [ES11, ... ]
#   - resource   [onwind, solar]
#   - resolution [LR, HR]
#   - format     [png, pdf]

rule plot_ISA:
    wildcard_constraints:
        nuts="NUTS0|NUTS2|NUTS3",
    benchmark:
        "benchmarks/plot_ISA/{nuts}/{resolution}/plot_ISA_{resource}_{region}_{resolution}.{format}.tsv"
    message:
        "... [plot_ISA] Plotting ISA map for resource: {wildcards.resource}, region: {wildcards.region}, resolution: {wildcards.resolution} format: {wildcards.format}."
    params:
        fig_params=config["fig_params"]
    input:
        gdf_NUTS ="data/NUTS/NUTS_RG_01M_2021_4326_ES.geojson",
        raster_ISA="results/rasters/ISA/{nuts}/raster_ISA_{resource}_{region}.tiff",
        df_ISA="results/dfs/ISA/{nuts}/df_ISA_{resource}_{region}.csv"
    output:
        map_ISA="results/maps/ISA/{nuts}/{resolution}/ISA_{resource}_{region}_{resolution}.{format}"
    script:
        "../scripts/plot_raster_ISA.py"



#################### plot_cutout
#
# Wildcards:
#   - cutout     [era5, ...]
#   - nuts       [NUTS2, NUTS3]
#   - region     [ES11, ... ]
#   - resource   [onwind, solar]
#   - year       [2013, ...]
#   - format     [png, pdf]

rule plot_cutout:
    wildcard_constraints:
        nuts="NUTS0|NUTS2|NUTS3",
    benchmark:
        "benchmarks/plot_cutout/{cutout}/{nuts}/plot_cutout_{resource}_{region}_{year}.{format}.tsv"
    message:
        "... [plot_cutout] Plotting cutout map for cutout: {wildcards.cutout}, year: {wildcards.year}, resource: {wildcards.resource}, region: {wildcards.region}, format: {wildcards.format}."
    params:
        cutout_params=config["cutout_params"],
        fig_params=config["fig_params"]
    input:
        gdf_NUTS ="data/NUTS/NUTS_RG_01M_2021_4326_ES.geojson",
    output:
        map_cutout="results/maps/cutout/{cutout}/{nuts}/cutout_{resource}_{region}_{year}.{format}"
    script:
        "../scripts/plot_cutout.py"        



#################### plot_CF
#
# NOTE:
#   - Color scale (vmin/vmax) is computed inside plot_nc_CF.py from all
#     df_CF_CAPACITY files for the same {cutout, nuts, resource, year}
#     and all regions in that NUTS level (ISA0..ISA4).
#   - This enforces a common absolute CF scale across regions.
#
# Wildcards:
#   - cutout     [era5, ...]
#   - nuts       [NUTS2, NUTS3]
#   - region     [ES11, ... ]
#   - resource   [onwind, solar]
#   - year       [2013, ...]
#   - format     [png, pdf]

rule plot_CF:
    wildcard_constraints:
        nuts="NUTS0|NUTS2|NUTS3",
    benchmark:
        "benchmarks/plot_CF/{cutout}/{nuts}/plot_CF_{resource}_{region}_{year}.{format}.tsv"
    message:
        "... [plot_CF] Plotting CF map for cutout: {wildcards.cutout}, year: {wildcards.year}, resource: {wildcards.resource}, region: {wildcards.region}, format: {wildcards.format}."
    params:
        fig_params=config["fig_params"]
    input:
        gdf_NUTS ="data/NUTS/NUTS_RG_01M_2021_4326_ES.geojson",
        nc_CF="results/ncs/CF/{cutout}/{nuts}/CF_{resource}_{region}_{year}.nc",
        dfs_CF_CAPACITY=lambda w: expand(
            "results/dfs/CF_CAPACITY/{cutout}/{nuts}/df_CF_CAPACITY_ISA{isa}_{resource}_{region}_{year}.csv",
            cutout=w.cutout,
            nuts=w.nuts,
            isa=[0, 1, 2, 3, 4],
            resource=w.resource,
            region=get_regions_for_nuts(w),
            year=w.year,
        )
    output:
        map_CF="results/maps/CF/{cutout}/{nuts}/CF_{resource}_{region}_{year}.{format}"
    script:
        "../scripts/plot_nc_CF.py"



#################### plot_CAPACITY
#
# Wildcards:
#   - cutout     [era5, ...]
#   - nuts       [NUTS2, NUTS3]
#   - filters    [CFth, ISA0, ... , CFth_ISA0, ...]
#   - region     [ES11, ... ]
#   - resource   [onwind, solar]
#   - year       [2013, ...]
#   - format     [png, pdf]

rule plot_CAPACITY:
    wildcard_constraints:
        nuts="NUTS0|NUTS2|NUTS3",
        resource="onwind|solar",
        filters=r"CFth|ISA\d+|CFth_ISA\d+",
    benchmark:
        "benchmarks/plot_CAPACITY/{cutout}/{nuts}/plot_CAPACITY_{filters}_{resource}_{region}_{year}.{format}.tsv"
    message:
        "... [plot_CAPACITY] Plotting CAPACITY map for cutout: {wildcards.cutout}, year: {wildcards.year}, resource: {wildcards.resource}, region: {wildcards.region}, filters: {wildcards.filters}, format: {wildcards.format}."
    params:
        fig_params=config["fig_params"]
    input:
        gdf_NUTS ="data/NUTS/NUTS_RG_01M_2021_4326_ES.geojson",
        nc_CAPACITY="results/ncs/CAPACITY/{cutout}/{nuts}/CAPACITY_{filters}_{resource}_{region}_{year}.nc"
    output:
        map_CAPACITY="results/maps/CAPACITY/{cutout}/{nuts}/CAPACITY_{filters}_{resource}_{region}_{year}.{format}"
    script:
        "../scripts/plot_nc_CAPACITY.py"



#################### plot_df_CF_CAPACITY
#
# Wildcards:
#   - cutout     [era5, ...]
#   - nuts       [NUTS2, NUTS3]
#   - region     [ES11, ... ]
#   - resource   [onwind, solar]
#   - year       [2013, ...]
#   - format     [png, pdf]

rule plot_df_CF_CAPACITY:
    wildcard_constraints:
        nuts="NUTS0|NUTS2|NUTS3",
    benchmark:
        "benchmarks/plot_df_CF_CAPACITY/{cutout}/{nuts}/plot_df_CF_CAPACITY_{resource}_{region}_{year}.{format}.tsv"
    message:
        "... [plot_df_CF_CAPACITY] Plotting CF vs CUM_CAPACITY curves for cutout: {wildcards.cutout}, year: {wildcards.year}, resource: {wildcards.resource}, region: {wildcards.region}, format: {wildcards.format}."
    params:
        fig_params=config["fig_params"]
    input:
        dfs_CF_CAPACITY=expand(
            "results/dfs/CF_CAPACITY/{{cutout}}/{{nuts}}/df_CF_CAPACITY_ISA{isa}_{{resource}}_{{region}}_{{year}}.csv",
            isa=[0, 1, 2, 3, 4]
        )
    output:
        plot_CF_CAPACITY="results/figs/CF_CAPACITY/{cutout}/{nuts}/CF_CAPACITY_{resource}_{region}_{year}.{format}"
    script:
        "../scripts/plot_df_CF_CAPACITY.py"



#################### plot_venn_single
#
# Wildcards:
#   - cutout     [era5, ...]
#   - nuts       [NUTS2, NUTS3]
#   - region     [ES11, ... ]
#   - resource   [onwind, solar]
#   - year       [2013, ...]
#   - format     [png, pdf]

rule plot_venn_single:
    wildcard_constraints:
        nuts="NUTS0|NUTS2|NUTS3",
    benchmark:
        "benchmarks/plot_venn_single/{cutout}/{nuts}/plot_venn_single_{resource}_{region}_{year}.{format}.tsv"
    message:
        "... [plot_venn_single] Plotting Venn diagram from df_summary for cutout: {wildcards.cutout}, year: {wildcards.year}, resource: {wildcards.resource}, region: {wildcards.region}, format: {wildcards.format}."
    params:
        fig_params=config["fig_params"]
    input:
        df_summary="results/dfs/summary/{cutout}/{nuts}/df_summary_{resource}_{year}.csv"
    output:
        plot_venn="results/figs/venn/{cutout}/{nuts}/venn_{resource}_{region}_{year}.{format}"
    script:
        "../scripts/plot_venn_single.py"