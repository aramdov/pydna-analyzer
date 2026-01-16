"""
GenomeInsight: A comprehensive personal genomics toolkit.

Analyze your DTC DNA data for clinical variants, pharmacogenomics,
polygenic risk scores, and ancestry estimation.
"""

__version__ = "0.1.0"
__author__ = "Aram Dovlatyan"

from genomeinsight.core.data_loader import load_dna_data
from genomeinsight.clinical.analyzer import ClinicalAnalyzer
from genomeinsight.clinical.apoe import determine_apoe_status

__all__ = [
    "load_dna_data",
    "ClinicalAnalyzer",
    "determine_apoe_status",
    "__version__",
]
