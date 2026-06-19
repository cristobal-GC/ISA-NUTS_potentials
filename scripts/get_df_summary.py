from pathlib import Path
import re

import pandas as pd

from typing import Any
snakemake: Any


############################## Unwrap relevant variables

##### input
files_df_CF_CAPACITY = list(snakemake.input["dfs_CF_CAPACITY"])
files_df_CAPACITY = list(snakemake.input["dfs_CAPACITY"])
##### params
regions = list(snakemake.params["regions"])
CF_threshold = float(snakemake.params["CF_threshold"])
cap_per_sqkm = float(snakemake.params["cap_per_sqkm"])
##### wildcards
nuts = snakemake.wildcards["nuts"]
resource = snakemake.wildcards["resource"]
##### output
file_df_summary = snakemake.output["df_summary"]


############################## Operations

# Region codes are NUTS ids (ES, ES11, ES111) or custom domain ids (e.g. DOM1):
# uppercase letters/digits with no underscore.
pattern = re.compile(
    r"df_CF_CAPACITY_ISA(?P<isa>\d+)_[^_]+_(?P<region>[A-Z0-9]+)_\d+\.csv$"
)
pattern_capacity = re.compile(
    r"df_CAPACITY_[^_]+_(?P<region>[A-Z0-9]+)_\d+\.csv$"
)

summary = {
    region: {
        **{f"CAPACITY_ISA{i}": 0.0 for i in range(5)},
        **{f"CAPACITY_CFth_ISA{i}": 0.0 for i in range(5)},
        "GENERATION_CFth_ISA4": 0.0,
    }
    for region in regions
}

for file_path in files_df_CF_CAPACITY:
    match = pattern.search(Path(file_path).name)
    if not match:
        raise ValueError(f"Could not parse ISA/region from file name: {file_path}")

    isa = int(match.group("isa"))
    region = match.group("region")

    if region not in summary:
        continue

    df = pd.read_csv(file_path)

    if df.empty:
        summary[region][f"CAPACITY_ISA{isa}"] = 0.0
        summary[region][f"CAPACITY_CFth_ISA{isa}"] = 0.0
        continue

    cum_capacity_col = "CUM_CAPACITY"
    cf_col = "CF"

    if cum_capacity_col not in df.columns or cf_col not in df.columns:
        raise ValueError(f"Missing required columns in {file_path}. Expected 'CF' and 'CUM_CAPACITY'.")

    cap_total = float(df[cum_capacity_col].iloc[-1])

    df_th = df[df[cf_col] >= CF_threshold]
    cap_th = float(df_th[cum_capacity_col].iloc[-1]) if not df_th.empty else 0.0

    summary[region][f"CAPACITY_ISA{isa}"] = round(cap_total, 6)
    summary[region][f"CAPACITY_CFth_ISA{isa}"] = round(cap_th, 6)

    if isa == 4:
        cap_col = "CAPACITY"
        if cap_col not in df.columns:
            raise ValueError(f"Missing required column '{cap_col}' in {file_path}.")

        GENERATION_th_twh = float((df_th[cap_col] * df_th[cf_col] * 8760 / 1000000).sum()) if not df_th.empty else 0.0
        summary[region]["GENERATION_CFth_ISA4"] = round(GENERATION_th_twh, 6)


def get_total_area_from_df_capacity(file_path):
    df = pd.read_csv(file_path)
    cols_upper = {col.upper(): col for col in df.columns}

    area_col = cols_upper.get("AREA")
    if area_col is None:
        raise ValueError(f"Missing AREA column in {file_path}.")

    value_col = cols_upper.get("VALUE")
    if value_col is None:
        first_col = df.columns[0]
        total_rows = df[df[first_col].astype(str).str.upper() == "TOTAL"]
    else:
        total_rows = df[df[value_col].astype(str).str.upper() == "TOTAL"]

    if total_rows.empty:
        raise ValueError(f"Row TOTAL not found in {file_path}.")

    return float(total_rows.iloc[0][area_col])


total_area_by_region = {}
for file_path in files_df_CAPACITY:
    match = pattern_capacity.search(Path(file_path).name)
    if not match:
        raise ValueError(f"Could not parse region from CAPACITY file name: {file_path}")
    region = match.group("region")
    total_area_by_region[region] = get_total_area_from_df_capacity(file_path)


cols_capacity = [f"CAPACITY_ISA{i}" for i in range(5)]
cols_capacity_th = [f"CAPACITY_CFth_ISA{i}" for i in range(5)]
cols_capacity_all = cols_capacity + cols_capacity_th + ["CAPACITY_CFth"]

summary_df = pd.DataFrame.from_dict(summary, orient="index")
summary_df = summary_df[cols_capacity + cols_capacity_th + ["GENERATION_CFth_ISA4"]]
summary_df["CAPACITY_CFth"] = summary_df[cols_capacity_th].sum(axis=1).round(6)

for region in summary_df.index:
    if region not in total_area_by_region:
        raise ValueError(f"Total area for region {region} not found in ISA inputs.")

for capacity_col in cols_capacity_all:
    suffix = capacity_col.replace("CAPACITY_", "")
    area_col = f"area_{suffix}"
    porc_area_col = f"porc_area_{suffix}"

    summary_df[area_col] = (summary_df[capacity_col] / cap_per_sqkm).round(6)

    summary_df[porc_area_col] = [
        round((summary_df.loc[region, area_col] / total_area_by_region[region]) * 100, 6)
        if total_area_by_region[region] > 0
        else 0.0
        for region in summary_df.index
    ]

cols_area = [f"area_{col.replace('CAPACITY_', '')}" for col in cols_capacity_all]
cols_porc_area = [f"porc_area_{col.replace('CAPACITY_', '')}" for col in cols_capacity_all]

summary_df = summary_df[cols_capacity_all + cols_area + cols_porc_area + ["GENERATION_CFth_ISA4"]]


def load_metric_by_region(file_path, value_col):
    metric = pd.read_csv(file_path)
    if "code" not in metric.columns:
        raise ValueError(f"Missing required column 'code' in {file_path}")
    if value_col not in metric.columns:
        raise ValueError(f"Missing required column '{value_col}' in {file_path}")

    metric["code"] = metric["code"].astype(str).str.strip().str.upper()
    return metric.set_index("code")[value_col].to_dict()


if nuts in {"NUTS0", "NUTS2"}:
    esios_dir = Path("data/esios")
    file_demand_2025 = esios_dir / f"electricity_demand_2025_{nuts}.csv"
    file_installed_2025 = esios_dir / f"esios_onwind_capacity_2025_{nuts}.csv"

    if not file_demand_2025.exists():
        raise FileNotFoundError(
            f"Required demand file for {nuts} not found: {file_demand_2025}"
        )

    demand_by_region = load_metric_by_region(file_demand_2025, "demand_2025")
    summary_df["DEMAND_2025"] = [
        round(float(demand_by_region.get(region, 0.0)), 6) for region in summary_df.index
    ]

    if resource == "onwind":
        if not file_installed_2025.exists():
            raise FileNotFoundError(
                f"Required installed capacity file for {nuts} not found: {file_installed_2025}"
            )
        installed_by_region = load_metric_by_region(file_installed_2025, "wind_capacity")
        summary_df["installed_2025"] = [
            round(float(installed_by_region.get(region, 0.0)), 6) for region in summary_df.index
        ]
    else:
        summary_df["installed_2025"] = 0.0

    summary_df["GENERATION_CFth_ISA4_perc_DEMAND_2025"] = [
        round((summary_df.loc[region, "GENERATION_CFth_ISA4"] / summary_df.loc[region, "DEMAND_2025"]) * 100, 6)
        if summary_df.loc[region, "DEMAND_2025"] > 0
        else 0.0
        for region in summary_df.index
    ]

summary_df.index.name = "region"


############################## Create outputs

Path(file_df_summary).parent.mkdir(parents=True, exist_ok=True)
summary_df.to_csv(file_df_summary)
