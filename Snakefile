
import csv
import math
from pathlib import Path
import re
from functools import lru_cache


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



##### Resource limits inferred from Snakemake benchmark TSVs
# Snakemake writes one TSV per benchmark path. Re-running the same job rewrites
# that file; different wildcard combinations create different files.
#
# Resource selection in this workflow follows three tiers:
#   1) use the exact benchmark file for the requested job if it exists,
#   2) otherwise use the worst observed benchmark across that rule,
#   3) otherwise fall back to a conservative static default.
#
# Benchmarks are not enforced as hard limits by themselves; we convert observed
# usage into Snakemake resources with a small safety margin.
#
# Important: rule-level resources such as mem_mb only affect scheduling when
# Snakemake is launched with a global resource budget, e.g.
#   snakemake all --cores 32 --resources mem_mb=230000
# Without --resources mem_mb=..., local execution still respects threads, but
# mem_mb is not used to cap overall concurrency.

BENCHMARK_MEM_MARGIN = 1.25
BENCHMARK_THREADS_MARGIN = 1.10

RULE_RESOURCE_DEFAULTS = {
    "get_raster_ISA": {"threads": 1, "mem_mb": 1024},
    "get_df_ISA": {"threads": 1, "mem_mb": 16384},
    "get_nc_CF": {"threads": 1, "mem_mb": 16384},
    "get_nc_CAPACITY": {"threads": 1, "mem_mb": 12288},
    "get_df_CF_CAPACITY": {"threads": 1, "mem_mb": 4096},
    "get_df_CAPACITY": {"threads": 1, "mem_mb": 4096},
    "get_df_summary": {"threads": 1, "mem_mb": 4096},
    "plot_GEBCO": {"threads": 1, "mem_mb": 2048},
    "plot_ISA": {"threads": 1, "mem_mb": 4096},
    "plot_cutout": {"threads": 1, "mem_mb": 4096},
    "plot_CF": {"threads": 1, "mem_mb": 4096},
    "plot_CAPACITY": {"threads": 1, "mem_mb": 4096},
    "plot_df_CF_CAPACITY": {"threads": 1, "mem_mb": 2048},
    "plot_venn_single": {"threads": 1, "mem_mb": 1024},
    "plot_potential_comparison": {"threads": 1, "mem_mb": 2048},
    "plot_NUTS": {"threads": 1, "mem_mb": 1024},
    "get_tex_summary": {"threads": 1, "mem_mb": 1024},
    "compile_tex_summary": {"threads": 1, "mem_mb": 1024},
    "get_tex_NUTS": {"threads": 1, "mem_mb": 1024},
    "compile_tex_NUTS": {"threads": 1, "mem_mb": 1024},
    "get_tex_cover": {"threads": 1, "mem_mb": 1024},
    "compile_tex_cover": {"threads": 1, "mem_mb": 1024},
    "get_tex_DECK": {"threads": 1, "mem_mb": 1024},
}


@lru_cache(maxsize=None)
def _read_benchmark_metrics(benchmark_file):
    # Benchmark TSVs contain a single data row with summary metrics for one job.
    path = Path(benchmark_file)
    if not path.exists():
        return None

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        row = next(reader, None)

    if row is None:
        return None

    def parse_float(field_name):
        try:
            value = row[field_name]
        except KeyError:
            return None

        if value in {None, ""}:
            return None

        try:
            return float(value)
        except ValueError:
            return None

    return {
        "seconds": parse_float("s"),
        "cpu_time": parse_float("cpu_time"),
        "max_rss": parse_float("max_rss"),
    }


@lru_cache(maxsize=None)
def _collect_benchmark_metrics(rule_name):
    # Aggregate all materialized benchmarks for one rule so new jobs can inherit
    # a pessimistic limit even before their exact benchmark exists.
    benchmark_dir = Path("benchmarks") / rule_name
    if not benchmark_dir.exists():
        return tuple()

    metrics = []
    for benchmark_file in sorted(benchmark_dir.rglob("*.tsv")):
        benchmark_metrics = _read_benchmark_metrics(str(benchmark_file))
        if benchmark_metrics is not None:
            metrics.append(benchmark_metrics)

    return tuple(metrics)


def _get_rule_benchmark_metrics(rule_name, benchmark_file):
    # Prefer the exact job benchmark. If it does not exist yet, fall back to the
    # available history for the whole rule.
    benchmark_metrics = _read_benchmark_metrics(benchmark_file)
    if benchmark_metrics is not None:
        return (benchmark_metrics,)

    return _collect_benchmark_metrics(rule_name)


def _estimate_threads_from_metrics(metrics, default_threads):
    # Approximate effective CPU parallelism from cpu_time / wall_time.
    seconds = metrics.get("seconds")
    cpu_time = metrics.get("cpu_time")

    if seconds is None or cpu_time is None or seconds <= 0 or cpu_time <= 0:
        return default_threads

    observed_threads = cpu_time / seconds
    return max(1, math.ceil(observed_threads * BENCHMARK_THREADS_MARGIN))


def _estimate_mem_mb_from_metrics(metrics, default_mem_mb):
    # max_rss is reported in MB by Snakemake benchmark TSVs.
    max_rss = metrics.get("max_rss")
    if max_rss is None or max_rss <= 0:
        return default_mem_mb

    return max(1, math.ceil(max_rss * BENCHMARK_MEM_MARGIN))


def get_rule_threads(rule_name, benchmark_file):
    defaults = RULE_RESOURCE_DEFAULTS[rule_name]
    benchmark_metrics = _get_rule_benchmark_metrics(rule_name, benchmark_file)
    if not benchmark_metrics:
        return defaults["threads"]

    return max(
        _estimate_threads_from_metrics(metrics, defaults["threads"])
        for metrics in benchmark_metrics
    )


def get_rule_mem_mb(rule_name, benchmark_file):
    defaults = RULE_RESOURCE_DEFAULTS[rule_name]
    benchmark_metrics = _get_rule_benchmark_metrics(rule_name, benchmark_file)
    if not benchmark_metrics:
        return defaults["mem_mb"]

    return max(
        _estimate_mem_mb_from_metrics(metrics, defaults["mem_mb"])
        for metrics in benchmark_metrics
    )





#################### RULES

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
            f"results/maps/potential_comparison/{cutout}/{nuts}/potential_comparison_{resource}_{year}.{fmt}"
            for nuts in sorted({n for n, _ in REGION_NUTS_PAIRS if n in {"NUTS2", "NUTS3"}})
            for cutout in CUTOUTS
            for year in YEARS
            for resource in RESOURCES
            for fmt in FORMATS
        ],

        [
            f"results/maps/GEBCO/{nuts}/GEBCO_{region}.{fmt}"
            for nuts, region in REGION_NUTS_PAIRS
            for fmt in FORMATS
        ],

        [
            f"results/maps/NUTS/{nuts}/NUTS_{nuts}.{fmt}"
            for nuts in ["NUTS2", "NUTS3"]
            for fmt in FORMATS
        ],

        [
            f"results/LaTex/{cutout}/{year}/{resource}/{nuts}/NUTS_{nuts}.pdf"
            for cutout in CUTOUTS
            for year in YEARS
            for resource in RESOURCES
            for nuts in ["NUTS2", "NUTS3"]
        ],

        [
            f"results/LaTex/{cutout}/{year}/{resource}/{nuts}/cover_{nuts}.pdf"
            for cutout in CUTOUTS
            for year in YEARS
            for resource in RESOURCES
            for nuts in ["NUTS2", "NUTS3"]
        ],

        [
            f"results/LaTex/{cutout}/{year}/{resource}/{nuts}/DECK_{nuts}.pdf"
            for cutout in CUTOUTS
            for year in YEARS
            for resource in RESOURCES
            for nuts in ["NUTS2", "NUTS3"]
        ],

        [
            f"results/LaTex/{cutout}/{year}/{resource}/{nuts}/summary_{region}_{resolution}.pdf"
            for nuts, region in REGION_NUTS_PAIRS
            for cutout in CUTOUTS
            for year in YEARS
            for resource in RESOURCES
            for resolution in RESOLUTIONS
        ]



rule retrieve:
    input:
        rules.retrieve_isa_onwind.output.tiff_file,
        rules.retrieve_isa_solar.output.tiff_file,
        rules.retrieve_gebco.output.gebco



rule plot_ISAs:
    input:
        [
            f"results/maps/ISA/{nuts}/{resolution}/ISA_{resource}_{region}_{resolution}.{fmt}"
            for nuts, region in REGION_NUTS_PAIRS
            for resource in RESOURCES
            for resolution in RESOLUTIONS
            for fmt in FORMATS
        ]


rule plot_GEBCOs:
    input:
        [
            f"results/maps/GEBCO/{nuts}/GEBCO_{region}.{fmt}"
            for nuts, region in REGION_NUTS_PAIRS
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


rule plot_potential_comparisons:
    input:
        [
            f"results/maps/potential_comparison/{cutout}/{nuts}/potential_comparison_{resource}_{year}.{fmt}"
            for nuts in sorted({n for n, _ in REGION_NUTS_PAIRS if n in {"NUTS2", "NUTS3"}})
            for cutout in CUTOUTS
            for year in YEARS
            for resource in RESOURCES
            for fmt in FORMATS
        ]


rule plot_NUTSs:
    input:
        [
            f"results/maps/NUTS/{nuts}/NUTS_{nuts}.{fmt}"
            for nuts in ["NUTS2", "NUTS3"]
            for fmt in FORMATS
        ]



# rule dag:
#     message:
#         "... Generating workflow DAG (PNG, PDF, SVG)"
#     output:
#         "DAG/dag.png",
#         "DAG/dag.pdf",
#         "DAG/dag.svg"
#     shell:
#         (
#             "mkdir -p DAG && "
#             "snakemake --dag --nolock | dot -Tpng -o {output[0]} && "
#             "snakemake --dag --nolock | dot -Tpdf -o {output[1]} && "
#             "snakemake --dag --nolock | dot -Tsvg -o {output[2]}"
#         )


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


# rule filegraph:
#     message:
#         "... Generating workflow rule graph (PNG, PDF, SVG)"
#     output:
#         "DAG/filegraph.png",
#         "DAG/filegraph.pdf",
#         "DAG/filegraph.svg"
#     shell:
#         (
#             "mkdir -p DAG && "
#             "snakemake --filegraph --nolock | dot -Tpng -o {output[0]} && "
#             "snakemake --filegraph --nolock | dot -Tpdf -o {output[1]} && "
#             "snakemake --filegraph --nolock | dot -Tsvg -o {output[2]}"
#         )


     