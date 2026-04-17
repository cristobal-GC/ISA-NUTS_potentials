from __future__ import annotations

import argparse
from pathlib import Path

import yaml


def format_number(value: object, decimals: int = 2, default: str = "N/A") -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default

    return f"{number:.{decimals}f}"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate NUTS LaTeX overview file from a resource-specific template."
    )
    parser.add_argument("--nuts", required=True, choices=["NUTS2", "NUTS3"])
    parser.add_argument("--cutout", required=True)
    parser.add_argument("--cutout-label", required=True)
    parser.add_argument("--year", required=True)
    parser.add_argument("--resource", required=True)
    parser.add_argument(
        "--project-root",
        default=str(Path(__file__).resolve().parents[1]),
        help="Repository root path",
    )
    args = parser.parse_args()

    project_root = Path(args.project_root).expanduser().resolve()
    template_dir = project_root / "LaTex"
    if args.resource == "onwind":
        template_file = template_dir / "template_NUTS_onwind.tex"
    else:
        template_file = template_dir / f"template_NUTS_{args.resource}.tex"
    config_file = project_root / "config" / "config.yaml"

    if not template_file.exists():
        if args.resource == "onwind":
            raise FileNotFoundError(f"Template not found: {template_file}")
        raise FileNotFoundError(
            f"Resource-specific template not found for '{args.resource}'. "
            f"Expected: {template_file}"
        )

    if not config_file.exists():
        raise FileNotFoundError(f"Config not found: {config_file}")

    with config_file.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    cf_params = config.get("CF_params", {}).get(args.resource, {})

    threshold = format_number(cf_params.get("CF_threshold"), decimals=3)
    cap_per_sqkm = format_number(cf_params.get("cap_per_sqkm"), decimals=0)
    correction_factor = format_number(cf_params.get("correction_factor"), decimals=2)
    turbine_model = str(cf_params.get("turbine", "N/A")).replace("_", " ")

    try:
        equiv_hours_raw = float(cf_params.get("CF_threshold")) * 8760
        equiv_hours = str(round(equiv_hours_raw / 100) * 100)
    except (TypeError, ValueError):
        equiv_hours = "N/A"

    tex_content = template_file.read_text(encoding="utf-8")
    replacements = {
        "ROOTPATH": str(project_root),
        "NUTSLEVEL": args.nuts,
        "CUTOUT": args.cutout_label,
        "YEAR": str(args.year),
        "THRESHOLD": threshold,
        "CAPPERSQKM": cap_per_sqkm,
        "CORRECTIONFACTOR": correction_factor,
        "TURBINEMODEL": turbine_model,
        "EQUIVHOURS": equiv_hours,
    }

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

    output_file = output_dir / f"NUTS_{args.nuts}.tex"
    output_file.write_text(tex_content, encoding="utf-8")

    print(f"LaTeX file generated: {output_file}")


if __name__ == "__main__":
    main()
