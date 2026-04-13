import argparse
import importlib
from pathlib import Path

import pandas as pd
import tabula as tb

PRESETS = {
    'cytokine_xl': {
        'pages': [17, 18, 19],
        'input_name': 'cytoXL array kit - protocol.pdf',
        'output_name': 'cytoXL array kit - protocol.xlsx',
    },
    'angiogenesis': {
        'pages': [15, 16],
        'input_name': 'Mouse Angiogenesis Array Kit - Protocol.pdf',
        'output_name': 'Mouse Angiogenesis Array Kit - Protocol.xlsx',
    },
}


def parse_args():
    parser = argparse.ArgumentParser(
        description='Extract analyte tables from a proteome profiler protocol PDF.',
    )
    parser.add_argument(
        '--preset',
        choices=sorted(PRESETS),
        default='cytokine_xl',
        help='Use one of the built-in protocol presets.',
    )
    parser.add_argument(
        '--input',
        help='Path to the protocol PDF. Defaults to the preset filename searched in the current directory and ./protocols/.',
    )
    parser.add_argument(
        '--output',
        help='Path to the output .xlsx file. Defaults to output/<preset workbook name>.',
    )
    parser.add_argument(
        '--pages',
        help='Comma-separated PDF pages to extract, overriding the preset.',
    )
    return parser.parse_args()


def resolve_input_path(input_arg: str | None, preset_input_name: str) -> Path:
    candidates = []
    if input_arg:
        candidates.append(Path(input_arg))
    else:
        candidates.extend([
            Path(preset_input_name),
            Path('protocols') / preset_input_name,
        ])

    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()

    raise FileNotFoundError(
        'Could not find the protocol PDF. Checked:\n'
        + '\n'.join(f'- {candidate}' for candidate in candidates),
    )


def read_protocol_tables(input_path: Path, pages: list[int]) -> list[pd.DataFrame]:
    """Read protocol tables using whichever tabula API variant is installed."""
    if hasattr(tb, 'read_pdf'):
        return tb.read_pdf(str(input_path), pages=pages, multiple_tables=True)

    io_module = getattr(tb, 'io', None)
    if io_module is None:
        try:
            io_module = importlib.import_module('tabula.io')
        except ImportError:
            io_module = None

    if io_module is not None and hasattr(io_module, 'read_pdf'):
        return io_module.read_pdf(str(input_path), pages=pages, multiple_tables=True)

    raise AttributeError(
        "Installed tabula module does not expose read_pdf. "
        "Expected tabula-py with tabula.read_pdf or tabula.io.read_pdf."
    )


def main():
    args = parse_args()
    preset = PRESETS[args.preset]
    input_path = resolve_input_path(args.input, preset['input_name'])

    if args.pages:
        pages = [int(page.strip()) for page in args.pages.split(',') if page.strip()]
    else:
        pages = preset['pages']

    if args.output:
        output_path = Path(args.output)
    else:
        output_path = Path('output') / preset['output_name']

    output_path.parent.mkdir(parents=True, exist_ok=True)

    tables = read_protocol_tables(input_path, pages)
    if not tables:
        raise RuntimeError(f'No tables were extracted from {input_path}')

    df = pd.concat(tables, ignore_index=True)
    df.to_excel(output_path, index=False)
    print(f'Wrote {output_path}')


if __name__ == '__main__':
    main()
