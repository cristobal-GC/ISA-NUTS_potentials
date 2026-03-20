
from pathlib import Path
import re


########## Set default config file if it exists
if Path("config/config.yaml").exists():
    configfile: "config/config.yaml"



##### Include rules
include: "rules/get.smk",
include: "rules/latex.smk"
include: "rules/plot.smk"
include: "rules/retrieve.smk"



##### Retrieve relevant information from config file
REGIONS_CFG = config["regions"]
RESOURCES = config["resources"]
CUTOUTS = config["cutouts"]
YEARS = config["years"]
FORMATS = config["formats"]
RESOLUTIONS = config["resolutions"]
ISAS = [0, 1, 2, 3, 4]



# This function infers NUTS level from explicit region patterns:
# NUTS0: AB   (A,B uppercase letters only, e.g. 'ES')
# NUTS2: ABxy (A,B uppercase letters; x,y digits)
# NUTS3: ABxyz (A,B uppercase letters; x,y,z digits)
def infer_nuts_level(region):
    if re.fullmatch(r"[A-Z]{2}", region):
        return "NUTS0"
    if re.fullmatch(r"[A-Z]{2}\d{2}", region):
        return "NUTS2"
    if re.fullmatch(r"[A-Z]{2}\d{3}", region):
        return "NUTS3"
    raise ValueError(
        f"Invalid region code '{region}'. Expected NUTS0 pattern AB, NUTS2 pattern ABxy or NUTS3 pattern ABxyz."
    )



# This generates REGION_NUTS_PAIRS = [(NUTS0, region0), (NUTS2, region1), (NUTS2, region2), ..., (NUTS3, regionX), ...]
if isinstance(REGIONS_CFG, dict):
    VALID_NUTS_KEYS = {"NUTS0", "NUTS2", "NUTS3"}
    invalid_keys = [key for key in REGIONS_CFG if key not in VALID_NUTS_KEYS]
    if invalid_keys:
        raise ValueError(
            f"Invalid regions keys {invalid_keys}. Expected only 'NUTS0', 'NUTS2' and/or 'NUTS3'."
        )

    REGION_NUTS_PAIRS = [
        (nuts, region)
        for nuts, regions in REGIONS_CFG.items()
        for region in (regions or [])
    ]

    for nuts, region in REGION_NUTS_PAIRS:
        inferred = infer_nuts_level(region)
        if inferred != nuts:
            raise ValueError(
                f"Region '{region}' does not match key '{nuts}'. It matches '{inferred}'."
            )
else:
    REGION_NUTS_PAIRS = [(infer_nuts_level(region), region) for region in REGIONS_CFG]

if not REGION_NUTS_PAIRS:
    raise ValueError(
        "No regions configured. Define regions in config/config.yaml as a list or under regions.NUTS0/NUTS2/NUTS3."
    )



# This generates FILTERS to apply to the CAPACITY matrix:
#   FILTERS = [CFth, ISA0, ... , ISA4, CFth_ISA0, ... , CFth_ISA4]
FILTERS = (
    ["CFth"]
    + [f"ISA{i}" for i in ISAS]
    + [f"CFth_ISA{i}" for i in ISAS]
)




rule all:
    input:
        #"DAG/dag.pdf",
        "DAG/rulegraph.pdf",
        #"DAG/filegraph.pdf",

        [
            f"results/maps/cutout/{cutout}/{nuts}/cutout_{resource}_{region}_{year}.{fmt}"
            for nuts, region in REGION_NUTS_PAIRS
            for cutout in CUTOUTS
            for year in YEARS
            for resource in RESOURCES
            for fmt in FORMATS
        ],
       
        [
            f"results/figs/venn/{cutout}/{nuts}/venn_{resource}_{region}_{year}.{fmt}"
            for nuts, region in REGION_NUTS_PAIRS
            for cutout in CUTOUTS
            for year in YEARS
            for resource in RESOURCES
            for fmt in FORMATS
        ],

        [
            f"results/LaTex/{cutout}/{year}/{resource}/{nuts}/summary_{region}_{resolution}.pdf"
            for nuts, region in REGION_NUTS_PAIRS
            for cutout in CUTOUTS
            for year in YEARS
            for resource in RESOURCES
            for resolution in RESOLUTIONS
        ]


rule plot_ISAs:
    input:
        [
            f"results/maps/ISA/{nuts}/{resolution}/ISA_{resource}_{region}_{resolution}.{fmt}"
            for nuts, region in REGION_NUTS_PAIRS
            for resource in RESOURCES
            for resolution in RESOLUTIONS
            for fmt in FORMATS
        ]


rule plot_CAPACITYs:
    input:
        [
            f"results/maps/CAPACITY/{cutout}/{nuts}/CAPACITY_{filters}_{resource}_{region}_{year}.{fmt}"
            for nuts, region in REGION_NUTS_PAIRS
            for filters in FILTERS
            for cutout in CUTOUTS
            for year in YEARS
            for resource in RESOURCES
            for fmt in FORMATS
        ]


rule plot_CFs:
    input:
        [
            f"results/maps/CF/{cutout}/{nuts}/CF_{resource}_{region}_{year}.{fmt}"
            for nuts, region in REGION_NUTS_PAIRS
            for cutout in CUTOUTS
            for year in YEARS
            for resource in RESOURCES
            for fmt in FORMATS
        ]



rule dag:
    message:
        "... Generating workflow DAG (PNG, PDF, SVG)"
    output:
        "DAG/dag.png",
        "DAG/dag.pdf",
        "DAG/dag.svg"
    shell:
        (
            "mkdir -p DAG && "
            "snakemake --dag --nolock | dot -Tpng -o {output[0]} && "
            "snakemake --dag --nolock | dot -Tpdf -o {output[1]} && "
            "snakemake --dag --nolock | dot -Tsvg -o {output[2]}"
        )


rule rulegraph:
    message:
        "... Generating workflow rule graph (PNG, PDF, SVG)"
    output:
        "DAG/rulegraph.png",
        "DAG/rulegraph.pdf",
        "DAG/rulegraph.svg"
    shell:
        (
            "mkdir -p DAG && "
            "snakemake --rulegraph --nolock | dot -Tpng -o {output[0]} && "
            "snakemake --rulegraph --nolock | dot -Tpdf -o {output[1]} && "
            "snakemake --rulegraph --nolock | dot -Tsvg -o {output[2]}"
        )


rule filegraph:
    message:
        "... Generating workflow rule graph (PNG, PDF, SVG)"
    output:
        "DAG/filegraph.png",
        "DAG/filegraph.pdf",
        "DAG/filegraph.svg"
    shell:
        (
            "mkdir -p DAG && "
            "snakemake --filegraph --nolock | dot -Tpng -o {output[0]} && "
            "snakemake --filegraph --nolock | dot -Tpdf -o {output[1]} && "
            "snakemake --filegraph --nolock | dot -Tsvg -o {output[2]}"
        )


     