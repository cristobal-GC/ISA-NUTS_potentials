from __future__ import annotations

import argparse
from pathlib import Path


def _get_pdf_merger():
    try:
        from pypdf import PdfWriter  # type: ignore

        return PdfWriter()
    except Exception:
        try:
            from PyPDF2 import PdfWriter  # type: ignore

            return PdfWriter()
        except Exception as exc:
            raise RuntimeError(
                "No PDF merger library found. Install 'pypdf' (recommended) or 'PyPDF2'."
            ) from exc


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Concatenate cover + NUTS + summary PDFs into a DECK PDF in the provided order."
        )
    )
    parser.add_argument("--output", required=True, help="Output DECK PDF path")
    parser.add_argument("inputs", nargs="+", help="Input PDFs in desired order")
    args = parser.parse_args()

    output_path = Path(args.output).expanduser().resolve()
    input_paths = [Path(p).expanduser().resolve() for p in args.inputs]

    missing = [str(p) for p in input_paths if not p.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing input PDFs for DECK concatenation:\n" + "\n".join(missing)
        )

    writer = _get_pdf_merger()

    for pdf_path in input_paths:
        with pdf_path.open("rb") as handle:
            writer.append(handle)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wb") as handle:
        writer.write(handle)

    print(f"DECK PDF generated: {output_path}")


if __name__ == "__main__":
    main()
