from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate Venn LaTeX file from a resource-specific template."
    )
    parser.add_argument("--nuts", required=True, choices=["NUTS2", "NUTS3"])
    parser.add_argument("--cutout", required=True)
    parser.add_argument("--year", required=True)
    parser.add_argument("--resource", required=True)
    parser.add_argument("--template", required=True)
    parser.add_argument(
        "--project-root",
        default=str(Path(__file__).resolve().parents[1]),
        help="Repository root path",
    )
    args = parser.parse_args()

    project_root = Path(args.project_root).expanduser().resolve()
    template_file = Path(args.template).expanduser()
    if not template_file.is_absolute():
        template_file = project_root / template_file

    if not template_file.exists():
        raise FileNotFoundError(f"Template not found: {template_file}")

    tex_content = template_file.read_text(encoding="utf-8")
    replacements = {
        "ROOTPATH": str(project_root),
        "NUTSLEVEL": args.nuts,
        "CUTOUTNAME": args.cutout,
        "RESOURCE": args.resource,
        "YEAR": str(args.year),
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

    output_file = output_dir / f"Venn_{args.nuts}.tex"
    output_file.write_text(tex_content, encoding="utf-8")

    print(f"LaTeX file generated: {output_file}")


if __name__ == "__main__":
    main()
