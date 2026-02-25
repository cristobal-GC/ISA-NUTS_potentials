
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



#################### get_ISA_local
#
# This rule is to generate an ISA raster for a specific region from the Spanis ISA raster
#
# Wildcards:
#   - region    [ES11, ... ]
#   - resource  [onwind, solar]

rule get_ISA_local:
    message:
        "... Getting ISA_local raster for resource {wildcards.resource} and region {wildcards.region}."
    input:
        gdf_NUTS ="data/NUTS/NUTS_RG_01M_2021_4326_ES.geojson",
        raster_ISA=get_file_ISA
    output:
        raster_ISA_local="results/rasters/ISA/ISA_local_{resource}_{region}.tiff"
    script:
        "../scripts/get_ISA_local.py"

