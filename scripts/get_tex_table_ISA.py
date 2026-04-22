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


def _get_row_value(df: pd.DataFrame, index_value: str, column_name: str) -> float | None:
    if column_name not in df.columns:
        return None

    if index_value in df.index:
        value = df.loc[index_value, column_name]
    else:
        return None

    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _extract_region_code(csv_file: Path, resource: str) -> str:
    prefix = f"df_ISA_{resource}_"
    stem = csv_file.stem
    if not stem.startswith(prefix):
        raise ValueError(f"Unexpected file name format: {csv_file.name}")
    return stem[len(prefix):]


def _build_table_rows(df_files: list[Path], resource: str, region_names: dict[str, str]) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []

    for csv_file in sorted(df_files):
        region_code = _extract_region_code(csv_file, resource)

        df = pd.read_csv(csv_file, index_col=0)
        df.index = df.index.map(str)

        area_km2 = _get_row_value(df, "TOTAL", "area")
        if area_km2 is None:
            area_km2 = 0.0
            for isa in range(5):
                isa_area = _get_row_value(df, str(isa), "area")
                if isa_area is not None:
                    area_km2 += isa_area

        isa_perc = []
        for isa in range(5):
            perc = _get_row_value(df, str(isa), "porc")
            isa_perc.append(perc if perc is not None else 0.0)

        rows.append(
            {
                "code": region_code,
                "name": region_names.get(region_code, region_code),
                "area_1e3_km2": area_km2 / 1000.0,
                "A_ISA0": isa_perc[0],
                "A_ISA1": isa_perc[1],
                "A_ISA2": isa_perc[2],
                "A_ISA3": isa_perc[3],
                "A_ISA4": isa_perc[4],
            }
        )

    return sorted(rows, key=lambda item: str(item["code"]))


def _render_table(nuts: str, rows: list[dict[str, float | str]]) -> str:
    nuts_level = nuts.removeprefix("NUTS")

    lines = [
        r"\begin{table}[t]",
        r"\centering",
        r"\small",
        rf"\caption{{Spain's regions surface for NUTS {nuts_level}, and percentage per ISA level}} \label{{table_isa_{nuts_level}}}",
        r"\begin{tabular}{llrrrrrr}",
        r"\hline",
        r"NUTS code & Name & \shortstack{Area \\ ($\times 10^3$ km$^2$)} & \shortstack{$A_{\mathrm{ISA0}}$ \\ (\%)} & \shortstack{$A_{\mathrm{ISA1}}$ \\ (\%)} & \shortstack{$A_{\mathrm{ISA2}}$ \\ (\%)} & \shortstack{$A_{\mathrm{ISA3}}$ \\ (\%)} & \shortstack{$A_{\mathrm{ISA4}}$ \\ (\%)} \\",
        r"\hline",
    ]

    for row in rows:
        code = _latex_escape(str(row["code"]))
        name = _latex_escape(str(row["name"]))
        area = float(row["area_1e3_km2"])
        a0 = float(row["A_ISA0"])
        a1 = float(row["A_ISA1"])
        a2 = float(row["A_ISA2"])
        a3 = float(row["A_ISA3"])
        a4 = float(row["A_ISA4"])

        lines.append(
            f"{code} & {name} & {area:.2f} & {a0:.2f} & {a1:.2f} & {a2:.2f} & {a3:.2f} & {a4:.2f} \\\\"  # noqa: E501
        )

    lines.extend([
        r"\hline",
        r"\end{tabular}",
        r"\end{table}",
        "",
    ])

    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate LaTeX table for ISA area distribution by NUTS level and resource."
    )
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

    df_dir = project_root / "results" / "dfs" / "ISA" / args.nuts
    pattern = f"df_ISA_{args.resource}_*.csv"
    df_files = sorted(df_dir.glob(pattern))
    if not df_files:
        raise FileNotFoundError(
            f"No ISA CSV files found for nuts={args.nuts}, resource={args.resource} in {df_dir} with pattern {pattern}"
        )

    nuts_geojson = (
        Path(args.nuts_geojson).expanduser().resolve()
        if args.nuts_geojson
        else project_root / "data" / "NUTS" / "NUTS_RG_01M_2021_4326_ES.geojson"
    )
    if not nuts_geojson.exists():
        raise FileNotFoundError(f"NUTS geojson not found: {nuts_geojson}")

    region_names = _read_region_names(nuts_geojson)
    rows = _build_table_rows(df_files, args.resource, region_names)
    tex_content = _render_table(args.nuts, rows)

    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text(tex_content, encoding="utf-8")

    print(f"LaTeX ISA table generated: {output_file}")


if __name__ == "__main__":
    main()
