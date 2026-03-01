"""
PyDNA Analyzer: A comprehensive personal genomics toolkit.

Analyze your DTC DNA data for clinical variants, pharmacogenomics,
polygenic risk scores, and ancestry estimation.
"""

__version__ = "0.1.0"
__author__ = "Aram Dovlatyan"

from pydna_analyzer.clinical.analyzer import ClinicalAnalyzer
from pydna_analyzer.clinical.apoe import determine_apoe_status
from pydna_analyzer.core.data_loader import load_dna_data
from pydna_analyzer.ancestry import AncestryAnalyzer, AncestryResult
from pydna_analyzer.pharmacogenomics import MetabolizerPhenotype, PGxAnalyzer, PGxResult

__all__ = [
    "load_dna_data",
    "ClinicalAnalyzer",
    "determine_apoe_status",
    "PGxAnalyzer",
    "PGxResult",
    "MetabolizerPhenotype",
    "AncestryAnalyzer",
    "AncestryResult",
    "__version__",
]
