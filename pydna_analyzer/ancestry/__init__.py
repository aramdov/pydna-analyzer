"""Ancestry estimation module.

Pipeline:
    DNADataset -> match AIMs -> build likelihood matrix -> MLE optimize
               -> bootstrap CI -> AncestryResult
"""

from __future__ import annotations

from pydna_analyzer.ancestry.estimator import AncestryEstimator, AncestryResult, PopulationResult
from pydna_analyzer.ancestry.reference_data import AIMDatabase
from pydna_analyzer.core.data_loader import DNADataset


class AncestryAnalyzer:
    """High-level ancestry estimation interface."""

    def __init__(self, n_bootstrap: int = 100) -> None:
        self.n_bootstrap = n_bootstrap
        self._db = AIMDatabase.load()
        self._estimator = AncestryEstimator()

    def analyze(self, dataset: DNADataset) -> AncestryResult:
        """Estimate ancestry composition from a DNA dataset."""
        aim_rsids = self._db.get_aim_rsids()
        user_rsids = dataset.get_rsids()
        overlap = aim_rsids & user_rsids

        genotypes: dict[str, str] = {}
        for rsid in overlap:
            gt = dataset.get_genotype(rsid)
            if gt is not None and gt not in ("--", "00", "NC", ""):
                genotypes[rsid] = gt

        return self._estimator.estimate(
            genotypes, self._db, n_bootstrap=self.n_bootstrap
        )


__all__ = [
    "AncestryAnalyzer",
    "AncestryResult",
    "PopulationResult",
    "AIMDatabase",
]
