
rule plot_region:
    input:
        gdf="data/NUTS/NUTS2_ES.geojson"
    output:
        map="results/figs/{region}.pdf"
    script:
        "../scripts/plot_region.py"
