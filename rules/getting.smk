

#################### get_ISA_local
#
# This rule is to generate a ISA raster for a specific region
#
# Wildcards:
#   - region    [ES11, ... ]
#   - resource  [onwind, solar]

rule get_ISA_local:
    input:
        file_NUTS = lambda wc: f"data/NUTS/{nuts_from_region(wc['region'])}_ES.geojson",


        file_ISA = lambda wc: (
                                "data/ISA/Clas_ISA_eol_pb.tiff" 
                                if wc['resource']=='onwind'
                                else 
                                "data/ISA/Clas_ISA_ftv_pb.tiff"
                               )
    output:
        file_ISA_local="results/ISA/ISA_local_{resource}_{region}.tiff"
    script:
        "../scripts/get_ISA_local.py"
