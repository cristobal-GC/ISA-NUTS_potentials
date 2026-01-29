

#################### This rule is to generate a ISA raster for a specific region
#
# Wildcards:
#   - region    [ES11, ... ]
#   - resource  [onwind, solar]

rule get_ISA_local:
    input:
        NUTS_file = lambda wc: f"data/{nuts_from_region(wc['region'])}.geojson"


        ISA_file = lambda wc: (
                                "data/ISA/Clas_ISA_eol_pb.tiff" 
                                if wc['resource']=='onwind'
                                else 
                                "data/ISA/Clas_ISA_ftv_pb.tiff"
                               )
    output:
        ISA_local_file="results/ISA/ISA_local_{resource}_{region}.tiff"
    script:
        "../scripts/get_ISA_local.py"
