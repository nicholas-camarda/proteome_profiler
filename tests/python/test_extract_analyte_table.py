from pathlib import Path

import pytest

import importlib.util


MODULE_PATH = Path("scripts/setup/extract_analyte_table.py")
SPEC = importlib.util.spec_from_file_location("extract_analyte_table", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_resolve_input_path_prefers_protocols_dir(tmp_path: Path) -> None:
    protocols_dir = tmp_path / "protocols"
    protocols_dir.mkdir()
    pdf_path = protocols_dir / "mock.pdf"
    pdf_path.write_text("x")

    cwd = Path.cwd()
    try:
        import os

        os.chdir(tmp_path)
        resolved = MODULE.resolve_input_path(None, "mock.pdf")
        assert resolved == pdf_path.resolve()
    finally:
        os.chdir(cwd)


def test_resolve_input_path_errors_when_missing(tmp_path: Path) -> None:
    cwd = Path.cwd()
    try:
        import os

        os.chdir(tmp_path)
        with pytest.raises(FileNotFoundError):
            MODULE.resolve_input_path(None, "missing.pdf")
    finally:
        os.chdir(cwd)


def test_read_protocol_tables_uses_tabula_io_when_top_level_missing(monkeypatch, tmp_path: Path) -> None:
    pdf_path = tmp_path / "mock.pdf"
    pdf_path.write_text("x")

    calls = []

    class FakeIO:
        @staticmethod
        def read_pdf(path, pages, multiple_tables):
            calls.append((path, tuple(pages), multiple_tables))
            return ["ok"]

    class FakeTabula:
        io = FakeIO()

    monkeypatch.setattr(MODULE, "tb", FakeTabula())

    result = MODULE.read_protocol_tables(pdf_path, [1, 2])

    assert result == ["ok"]
    assert calls == [(str(pdf_path), (1, 2), True)]
