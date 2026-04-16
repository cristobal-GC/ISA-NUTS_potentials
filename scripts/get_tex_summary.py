from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import yaml


def normalize_resolution(value: str) -> str:
    mapping = {
        "hr": "HR",
        "lr": "LR",
        "highres": "HR",
        "lowres": "LR",
    }
    key = value.strip().lower()
    if key not in mapping:
        raise ValueError("resolution must be one of: HR, LR, HighRes, LowRes")
    return mapping[key]


def resolve_fig_ext(resolution: str, fig_ext: str | None) -> str:
    if fig_ext:
        return fig_ext
    return "pdf" if resolution == "HR" else "png"


def read_summary(project_root: Path, cutout: str, nuts: str, resource: str, year: int) -> pd.DataFrame:
    csv_path = project_root / "results" / "dfs" / "summary" / cutout / nuts / f"df_summary_{resource}_{year}.csv"
    df_path = project_root / "results" / "dfs" / "summary" / cutout / nuts / f"df_summary_{resource}_{year}.df"

    if csv_path.exists():
        summary_file = csv_path
    elif df_path.exists():
        summary_file = df_path
    else:
        raise FileNotFoundError(
            f"Summary file not found. Expected one of: {csv_path} or {df_path}"
        )

    df = pd.read_csv(summary_file)
    if "region" in df.columns:
        df = df.set_index("region")
    else:
        first_col = df.columns[0]
        if first_col.lower() in {"region", "unnamed: 0"}:
            df = df.set_index(first_col)

    return df


def read_region_name(nuts_geojson: Path, region: str) -> str:
    with nuts_geojson.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    for feature in data.get("features", []):
        props = feature.get("properties", {})
        if props.get("NUTS_ID") == region:
            return props.get("NUTS_NAME") or props.get("NAME_LATN") or region

    return region


def read_threshold(config_file: Path, resource: str, cli_threshold: float | None) -> float:
    if cli_threshold is not None:
        return cli_threshold

    with config_file.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    try:
        return float(config["CF_params"][resource]["CF_threshold"])
    except KeyError as exc:
        raise KeyError(
            f"CF_threshold for resource='{resource}' not found in {config_file}"
        ) from exc


def format_number(value: object, decimals: int = 2, default: str = "N/A") -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default

    if pd.isna(number):
        return default

    return f"{number:.{decimals}f}"


def get_row_value(row: pd.Series, candidates: list[str]) -> object:
    for col in candidates:
        if col in row.index:
            return row[col]
    return None


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate a LaTeX summary file for one region/cutout/year/resource."
    )
    parser.add_argument("--region", required=True)
    parser.add_argument("--cutout", required=True)
    parser.add_argument("--year", required=True, type=int)
    parser.add_argument("--resource", required=True)
    parser.add_argument("--nuts", required=True, help="e.g. NUTS2")
    parser.add_argument("--resolution", required=True, help="HR/LR or HighRes/LowRes")
    parser.add_argument("--fig-ext", default=None, help="Override figure extension (png, pdf, ...) ")
    parser.add_argument("--threshold", type=float, default=None, help="Override CF threshold")

    parser.add_argument(
        "--project-root",
        default=str(Path(__file__).resolve().parents[1]),
        help="Repository root path",
    )
    parser.add_argument(
        "--template",
        default=None,
        help="Path to LaTeX template. Defaults to <project-root>/LaTex/template.tex",
    )
    parser.add_argument(
        "--config",
        default=None,
        help="Path to config YAML. Defaults to <project-root>/config/config.yaml",
    )
    parser.add_argument(
        "--nuts-geojson",
        default=None,
        help="Path to NUTS geojson. Defaults to <project-root>/data/NUTS/NUTS_RG_01M_2021_4326_ES.geojson",
    )

    args = parser.parse_args()

    project_root = Path(args.project_root).expanduser().resolve()
    template_file = (
        Path(args.template).expanduser().resolve()
        if args.template
        else project_root / "LaTex" / "template.tex"
    )
    config_file = (
        Path(args.config).expanduser().resolve()
        if args.config
        else project_root / "config" / "config.yaml"
    )
    nuts_geojson = (
        Path(args.nuts_geojson).expanduser().resolve()
        if args.nuts_geojson
        else project_root / "data" / "NUTS" / "NUTS_RG_01M_2021_4326_ES.geojson"
    )

    resolution = normalize_resolution(args.resolution)
    fig_ext = resolve_fig_ext(resolution, args.fig_ext)
    threshold = read_threshold(config_file, args.resource, args.threshold)

    df_summary = read_summary(project_root, args.cutout, args.nuts, args.resource, args.year)
    if args.region not in df_summary.index:
        raise KeyError(
            f"Region '{args.region}' not found in summary index: {list(df_summary.index)}"
        )

    row = df_summary.loc[args.region]

    porc_land = row.get("porc_area_CFth_ISA4")
    potcap_mw = row.get("CAPACITY_CFth_ISA4")
    potcap_gw = float(potcap_mw) / 1000 if pd.notna(potcap_mw) else None

    pot_actual_mw = get_row_value(row, ["installed_2025"])
    pot_actual = float(pot_actual_mw) / 1000.0 if pd.notna(pot_actual_mw) else None
    potene_twh = get_row_value(row, ["GENERATION_CFth_ISA4", "gen_CFth_ISA4", "ENERGY_potential_TWh"])
    perc_demand = get_row_value(
        row,
        ["GENERATION_CFth_ISA4_perc_DEMAND_2025", "gen_CFth_ISA4_perc_demand_2025", "perc_demand"],
    )
    demand_2025 = get_row_value(row, ["DEMAND_2025", "demand_2025"])

    region_name = read_region_name(nuts_geojson, args.region)

    template = template_file.read_text(encoding="utf-8")
    tex_content = template

    replacements = {
        "CODIGONUTS": args.region,
        "RESOLUTION": resolution,
        "FIGEXT": fig_ext,
        "THRESHOLD": format_number(threshold, decimals=2),
        "CUTOUT": args.cutout,
        "ROOTPATH": str(project_root),
        "RESOURCE": args.resource,
        "PORCLAND": format_number(porc_land, decimals=2),
        "POTCAP": format_number(potcap_gw, decimals=2),
        "POTACTUAL": format_number(pot_actual, decimals=2),
        "POTENE": format_number(potene_twh, decimals=2),
        "PERCDEMAND": format_number(perc_demand, decimals=0),
        "REGION": region_name,
        "NUTSLEVEL": args.nuts,
        "YEAR": str(args.year),
    }

    include_demand_block = args.nuts in {"NUTS0", "NUTS2"}
    if include_demand_block:
        pot_actual_str = format_number(pot_actual, decimals=2)
        potene_str = format_number(potene_twh, decimals=2)
        perc_demand_str = format_number(perc_demand, decimals=0)
        demand_2025_str = format_number(demand_2025, decimals=2)

        installed_note = (
            f"{{\\color{{gris}} Installed wind power capacity in {region_name} in 2025 was "
            f"\\textbf{{{pot_actual_str} GW}}.}} \\vspace{{0.4cm}}"
        )
        energy_note = (
            f"\\item \\texttt{{\\textbf{{ENERGY:}}}} {{\\color{{low}}\\textbf{{{potene_str} TWh}}}} could be generated with this potential capacity.\n"
            f"\n"
            f"        {{\\color{{gris}} This represents around \\textbf{{{perc_demand_str}\\%}} of the actual electricity demand in the region ({demand_2025_str} TWh).}}"
        )
    else:
        installed_note = "% Installed capacity note omitted for this NUTS level"
        energy_note = "% Energy/demand note omitted for this NUTS level"

    replacements["INSTALLED_NOTE"] = installed_note
    replacements["ENERGY_NOTE"] = energy_note

    for key, value in replacements.items():
        tex_content = tex_content.replace(key, value)

    output_dir = (
        project_root
        / "results"
        / "LaTex"
        / args.cutout
        / str(args.year)
        / args.resource
        / args.nuts
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"summary_{args.region}_{resolution}.tex"
    output_file.write_text(tex_content, encoding="utf-8")

    print(f"LaTeX file generated: {output_file}")


if __name__ == "__main__":
    main()
