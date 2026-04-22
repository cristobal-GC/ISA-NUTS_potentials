from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def _latex_escape(text: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(ch, ch) for ch in text)


def _read_region_names(nuts_geojson: Path) -> dict[str, str]:
    with nuts_geojson.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    region_names: dict[str, str] = {}
    for feature in data.get("features", []):
        props = feature.get("properties", {})
        nuts_id = props.get("NUTS_ID")
        if not nuts_id:
            continue

        name = (
            props.get("NUTS_NAME")
            or props.get("NAME_LATN")
            or props.get("NAME_ENGL")
            or nuts_id
        )
        region_names[str(nuts_id)] = str(name)

    return region_names


def _format_number(value: object, decimals: int = 2, scale: float = 1.0, default: str = "N/A") -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default

    if pd.isna(number):
        return default

    return f"{number / scale:.{decimals}f}"


def _load_summary(df_summary_file: Path) -> pd.DataFrame:
    df = pd.read_csv(df_summary_file)
    if "region" in df.columns:
        df = df.set_index("region")
    else:
        first_col = df.columns[0]
        if first_col.lower() in {"region", "unnamed: 0"}:
            df = df.set_index(first_col)

    df.index = df.index.map(str)
    return df.sort_index()


def _load_nuts0_total_row(project_root: Path, cutout: str, resource: str, year: int) -> pd.Series:
    df_summary_file = (
        project_root
        / "results"
        / "dfs"
        / "summary"
        / cutout
        / "NUTS0"
        / f"df_summary_{resource}_{year}.csv"
    )
    if not df_summary_file.exists():
        raise FileNotFoundError(f"NUTS0 summary file not found: {df_summary_file}")

    df_summary = _load_summary(df_summary_file)
    if "ES" not in df_summary.index:
        raise KeyError(f"Region 'ES' not found in NUTS0 summary file: {df_summary_file}")

    return df_summary.loc["ES"]


def _build_row_cells(nuts: str, region: str, row: pd.Series, region_names: dict[str, str]) -> list[str]:
    cells = [
        _latex_escape(region),
        _latex_escape(region_names.get(region, region)),
        _format_number(row.get("CAPACITY_ISA4"), scale=1000.0),
        _format_number(row.get("CAPACITY_CFth"), scale=1000.0),
        _format_number(row.get("CAPACITY_CFth_ISA4"), scale=1000.0),
        _format_number(row.get("porc_area_CFth_ISA4")),
    ]

    if nuts == "NUTS2":
        cells.append(_format_number(row.get("installed_2025"), scale=1000.0))

    cells.append(_format_number(row.get("GENERATION_CFth_ISA4")))

    if nuts == "NUTS2":
        cells.append(_format_number(row.get("GENERATION_CFth_ISA4_perc_DEMAND_2025")))

    return cells


def _render_table(
    nuts: str,
    df_summary: pd.DataFrame,
    region_names: dict[str, str],
    nuts0_total_row: pd.Series,
) -> str:
    nuts_level = nuts.removeprefix("NUTS")

    header_cells = [
        r"\shortstack{NUTS \\ code}",
        r"Name",
        r"\shortstack{$P_{_{\mathrm{ISA}4}}$ \\ {[GW]}}",
        r"\shortstack{$P^{^{\mathit{CF}^*}}$ \\ {[GW]}}",
        r"\shortstack{$P_{_{\mathrm{ISA}4}}^{^{\mathit{CF}^*}}$ \\ {[GW]}}",
        r"\shortstack{$A_{_{\mathrm{ISA}4}}^{^{\mathit{CF}^*}}$ \\ {$(\%)$}}",
    ]

    if nuts == "NUTS2":
        header_cells.append(r"\shortstack{$P_{_{2025}}$ \\ {[GW]}}")

    header_cells.append(r"\shortstack{$E_{_{\mathrm{ISA}4}}^{^{\mathit{CF}^*}}$ \\ {[TWh]}}")

    if nuts == "NUTS2":
        header_cells.append(r"\shortstack{$E_{_{\mathrm{ISA}4}}^{^{\mathit{CF}^*}}$ \\ {$(\%_{2025})$}}")

    col_spec = "ll" + "r" * (len(header_cells) - 2)

    lines = [
        r"\begin{table}[t]",
        r"\centering",
        r"\small",
        rf"\caption{{Summary of results for NUTS {nuts_level}.}} \label{{table_results_{nuts_level}}}",
        rf"\begin{{tabular}}{{{col_spec}}}",
        r"\hline",
        " & ".join(header_cells) + r" \\",
        r"\hline",
    ]

    for region, row in df_summary.iterrows():
        lines.append(" & ".join(_build_row_cells(nuts, region, row, region_names)) + r" \\")

    lines.extend([
        r"\hline",
        " & ".join(_build_row_cells(nuts, "ES", nuts0_total_row, region_names)) + r" \\",
        r"\hline",
        r"\end{tabular}",
        r"\end{table}",
        "",
    ])

    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate LaTeX summary table for one cutout/year/resource/NUTS level."
    )
    parser.add_argument("--cutout", required=True)
    parser.add_argument("--year", required=True, type=int)
    parser.add_argument("--nuts", required=True, choices=["NUTS2", "NUTS3"])
    parser.add_argument("--resource", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--project-root",
        default=str(Path(__file__).resolve().parents[1]),
        help="Repository root path",
    )
    parser.add_argument(
        "--nuts-geojson",
        default=None,
        help="Path to NUTS geojson. Defaults to <project-root>/data/NUTS/NUTS_RG_01M_2021_4326_ES.geojson",
    )
    args = parser.parse_args()

    project_root = Path(args.project_root).expanduser().resolve()
    output_file = Path(args.output).expanduser()
    if not output_file.is_absolute():
        output_file = project_root / output_file

    df_summary_file = (
        project_root
        / "results"
        / "dfs"
        / "summary"
        / args.cutout
        / args.nuts
        / f"df_summary_{args.resource}_{args.year}.csv"
    )
    if not df_summary_file.exists():
        raise FileNotFoundError(f"Summary file not found: {df_summary_file}")

    nuts_geojson = (
        Path(args.nuts_geojson).expanduser().resolve()
        if args.nuts_geojson
        else project_root / "data" / "NUTS" / "NUTS_RG_01M_2021_4326_ES.geojson"
    )
    if not nuts_geojson.exists():
        raise FileNotFoundError(f"NUTS geojson not found: {nuts_geojson}")

    df_summary = _load_summary(df_summary_file)
    nuts0_total_row = _load_nuts0_total_row(project_root, args.cutout, args.resource, args.year)
    region_names = _read_region_names(nuts_geojson)
    tex_content = _render_table(args.nuts, df_summary, region_names, nuts0_total_row)

    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text(tex_content, encoding="utf-8")

    print(f"LaTeX summary table generated: {output_file}")


if __name__ == "__main__":
    main()
