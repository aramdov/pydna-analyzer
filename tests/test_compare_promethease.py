"""Tests for compare_promethease helpers."""

import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


def _load_compare_module():
    module_path = Path(__file__).resolve().parent.parent / "compare_promethease.py"
    spec = spec_from_file_location("compare_promethease", module_path)
    assert spec is not None
    assert spec.loader is not None
    module = module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_load_pydna_from_dna_file_uses_live_analysis(tmp_path):
    """The comparison script should analyze the provided DNA file directly."""
    compare_promethease = _load_compare_module()

    content = """#AncestryDNA raw data download
#rsid\tchromosome\tposition\tallele1\tallele2
rs1801133\t1\t11856378\tG\tA
rs6025\t1\t169519049\tC\tC
"""
    dna_file = tmp_path / "sample.txt"
    dna_file.write_text(content)

    variants = compare_promethease.load_pydna_from_dna_file(dna_file)
    by_rsid = {variant.rsid: variant for variant in variants}

    assert by_rsid["rs1801133"].risk_level == "moderate"
    assert by_rsid["rs6025"].risk_level == "normal"
