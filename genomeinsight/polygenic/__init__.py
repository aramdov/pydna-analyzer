"""
Polygenic Risk Score (PRS) Calculator.

Calculates polygenic risk scores from user genotypes and published
SNP weight files (e.g., from PGS Catalog).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from genomeinsight.core.data_loader import DNADataset


class PRSNormalization(str, Enum):
    """Normalization methods for PRS."""

    RAW = "raw"  # Raw weighted sum
    ZSCORE = "zscore"  # Z-score normalized
    PERCENTILE = "percentile"  # Percentile rank (requires population stats)


@dataclass
class SNPWeight:
    """A single SNP weight entry."""

    rsid: str
    effect_allele: str
    weight: float
    other_allele: Optional[str] = None
    chromosome: Optional[str] = None
    position: Optional[int] = None


@dataclass
class PRSResult:
    """Result from PRS calculation."""

    score_name: str
    raw_score: float
    normalized_score: Optional[float] = None
    percentile: Optional[float] = None
    risk_category: str = "unknown"
    snps_used: int = 0
    snps_available: int = 0
    coverage: float = 0.0
    missing_snps: list[str] = field(default_factory=list)
    interpretation: str = ""

    @property
    def coverage_percent(self) -> str:
        """Return coverage as formatted percentage."""
        return f"{self.coverage * 100:.1f}%"


class PRSWeightLoader:
    """Load PRS weights from various formats."""

    @staticmethod
    def load_csv(
        filepath: Path | str,
        rsid_col: str = "rsid",
        effect_allele_col: str = "effect_allele",
        weight_col: str = "weight",
        other_allele_col: Optional[str] = None,
    ) -> list[SNPWeight]:
        """
        Load weights from CSV file.

        Expected format:
            rsid,effect_allele,weight
            rs123456,A,0.05
            rs789012,G,-0.03

        Args:
            filepath: Path to CSV file
            rsid_col: Column name for rsID
            effect_allele_col: Column name for effect allele
            weight_col: Column name for weight/beta
            other_allele_col: Optional column for other allele

        Returns:
            List of SNPWeight objects
        """
        filepath = Path(filepath)
        df = pd.read_csv(filepath)

        # Normalize column names
        df.columns = df.columns.str.lower().str.strip()
        rsid_col = rsid_col.lower()
        effect_allele_col = effect_allele_col.lower()
        weight_col = weight_col.lower()

        if rsid_col not in df.columns:
            # Try common alternatives
            for alt in ["snp", "variant", "variant_id", "snp_id"]:
                if alt in df.columns:
                    rsid_col = alt
                    break

        if effect_allele_col not in df.columns:
            for alt in ["allele", "effect", "alt", "a1", "allele1"]:
                if alt in df.columns:
                    effect_allele_col = alt
                    break

        if weight_col not in df.columns:
            for alt in ["beta", "effect_weight", "effect_size", "or", "odds_ratio"]:
                if alt in df.columns:
                    weight_col = alt
                    break

        weights = []
        for _, row in df.iterrows():
            other_allele = None
            if other_allele_col and other_allele_col.lower() in df.columns:
                other_allele = str(row[other_allele_col.lower()])

            weights.append(
                SNPWeight(
                    rsid=str(row[rsid_col]).strip(),
                    effect_allele=str(row[effect_allele_col]).strip().upper(),
                    weight=float(row[weight_col]),
                    other_allele=other_allele,
                )
            )

        return weights

    @staticmethod
    def load_pgs_catalog_format(filepath: Path | str) -> list[SNPWeight]:
        """
        Load PGS Catalog scoring file format.

        PGS Catalog files typically have columns:
        rsID, chr_name, chr_position, effect_allele, other_allele, effect_weight
        """
        filepath = Path(filepath)

        # Skip header lines starting with #
        with open(filepath) as f:
            skip_rows = 0
            for line in f:
                if line.startswith("#"):
                    skip_rows += 1
                else:
                    break

        df = pd.read_csv(filepath, skiprows=skip_rows, sep="\t")
        df.columns = df.columns.str.lower().str.strip()

        weights = []
        for _, row in df.iterrows():
            rsid = row.get("rsid", row.get("snp", ""))
            if pd.isna(rsid) or not str(rsid).startswith("rs"):
                continue

            weights.append(
                SNPWeight(
                    rsid=str(rsid).strip(),
                    effect_allele=str(row.get("effect_allele", row.get("allele1", ""))).upper(),
                    weight=float(row.get("effect_weight", row.get("weight", row.get("beta", 0)))),
                    other_allele=str(row.get("other_allele", row.get("allele2", ""))).upper()
                    if pd.notna(row.get("other_allele", row.get("allele2")))
                    else None,
                    chromosome=str(row.get("chr_name", row.get("chr", "")))
                    if pd.notna(row.get("chr_name", row.get("chr")))
                    else None,
                    position=int(row.get("chr_position", row.get("pos", 0)))
                    if pd.notna(row.get("chr_position", row.get("pos")))
                    else None,
                )
            )

        return weights


class PRSCalculator:
    """
    Calculate polygenic risk scores.

    Example:
        >>> from genomeinsight import load_dna_data
        >>> from genomeinsight.polygenic import PRSCalculator
        >>> 
        >>> dataset = load_dna_data("my_dna.txt")
        >>> calculator = PRSCalculator()
        >>> 
        >>> # Load weights and calculate
        >>> result = calculator.calculate_from_file(
        ...     dataset,
        ...     "heart_disease_prs.csv",
        ...     score_name="Coronary Artery Disease"
        ... )
        >>> print(f"Score: {result.raw_score:.4f}")
        >>> print(f"Coverage: {result.coverage_percent}")
    """

    def __init__(
        self,
        population_mean: Optional[float] = None,
        population_std: Optional[float] = None,
    ):
        """
        Initialize calculator.

        Args:
            population_mean: Mean score in reference population (for z-score)
            population_std: Standard deviation in reference population
        """
        self.population_mean = population_mean
        self.population_std = population_std

    def calculate(
        self,
        dataset: DNADataset,
        weights: list[SNPWeight],
        score_name: str = "PRS",
    ) -> PRSResult:
        """
        Calculate PRS from dataset and weights.

        The score is calculated as:
            PRS = Σ (dosage × weight)

        Where dosage is:
            - 0 if genotype has 0 copies of effect allele
            - 1 if genotype has 1 copy (heterozygous)
            - 2 if genotype has 2 copies (homozygous effect)

        Args:
            dataset: DNADataset with genotypes
            weights: List of SNPWeight objects
            score_name: Name for this score

        Returns:
            PRSResult with calculated score and metadata
        """
        raw_score = 0.0
        snps_used = 0
        missing_snps = []

        for snp_weight in weights:
            genotype = dataset.get_genotype(snp_weight.rsid)

            if genotype is None or genotype in ("--", "00", "NC", ""):
                missing_snps.append(snp_weight.rsid)
                continue

            # Calculate dosage (count of effect allele)
            dosage = self._calculate_dosage(genotype, snp_weight.effect_allele)
            raw_score += dosage * snp_weight.weight
            snps_used += 1

        # Calculate coverage
        coverage = snps_used / len(weights) if weights else 0.0

        # Normalize if population stats available
        normalized_score = None
        percentile = None
        if self.population_mean is not None and self.population_std is not None:
            normalized_score = (raw_score - self.population_mean) / self.population_std
            # Approximate percentile from z-score (assuming normal distribution)
            from scipy import stats

            percentile = stats.norm.cdf(normalized_score) * 100

        # Determine risk category
        risk_category = self._categorize_risk(normalized_score, percentile)

        # Generate interpretation
        interpretation = self._generate_interpretation(
            score_name, raw_score, normalized_score, percentile, coverage
        )

        return PRSResult(
            score_name=score_name,
            raw_score=raw_score,
            normalized_score=normalized_score,
            percentile=percentile,
            risk_category=risk_category,
            snps_used=snps_used,
            snps_available=len(weights),
            coverage=coverage,
            missing_snps=missing_snps[:20],  # Limit to first 20
            interpretation=interpretation,
        )

    def calculate_from_file(
        self,
        dataset: DNADataset,
        weights_file: Path | str,
        score_name: Optional[str] = None,
        file_format: str = "auto",
    ) -> PRSResult:
        """
        Calculate PRS from a weights file.

        Args:
            dataset: DNADataset with genotypes
            weights_file: Path to weights file (CSV or PGS Catalog format)
            score_name: Name for this score (defaults to filename)
            file_format: "csv", "pgs", or "auto"

        Returns:
            PRSResult with calculated score
        """
        weights_file = Path(weights_file)

        if score_name is None:
            score_name = weights_file.stem.replace("_", " ").title()

        # Detect format
        if file_format == "auto":
            # Check if it's PGS Catalog format (has # header)
            with open(weights_file) as f:
                first_line = f.readline()
                if first_line.startswith("#"):
                    file_format = "pgs"
                else:
                    file_format = "csv"

        # Load weights
        if file_format == "pgs":
            weights = PRSWeightLoader.load_pgs_catalog_format(weights_file)
        else:
            weights = PRSWeightLoader.load_csv(weights_file)

        return self.calculate(dataset, weights, score_name)

    def _calculate_dosage(self, genotype: str, effect_allele: str) -> int:
        """
        Calculate dosage (0, 1, or 2) of effect allele in genotype.

        Args:
            genotype: Two-character genotype (e.g., "AG", "CC")
            effect_allele: The allele to count

        Returns:
            Count of effect allele (0, 1, or 2)
        """
        genotype = genotype.upper()
        effect_allele = effect_allele.upper()

        count = 0
        for allele in genotype:
            if allele == effect_allele:
                count += 1

        return count

    def _categorize_risk(
        self,
        normalized_score: Optional[float],
        percentile: Optional[float],
    ) -> str:
        """Categorize risk based on score or percentile."""
        if percentile is not None:
            if percentile >= 95:
                return "very_high"
            elif percentile >= 80:
                return "high"
            elif percentile >= 60:
                return "moderate"
            elif percentile >= 40:
                return "average"
            elif percentile >= 20:
                return "low"
            else:
                return "very_low"
        elif normalized_score is not None:
            if normalized_score >= 2:
                return "very_high"
            elif normalized_score >= 1:
                return "high"
            elif normalized_score >= 0:
                return "moderate"
            elif normalized_score >= -1:
                return "average"
            else:
                return "low"
        else:
            return "unknown"

    def _generate_interpretation(
        self,
        score_name: str,
        raw_score: float,
        normalized_score: Optional[float],
        percentile: Optional[float],
        coverage: float,
    ) -> str:
        """Generate human-readable interpretation."""
        parts = []

        if percentile is not None:
            parts.append(
                f"Your {score_name} score places you in the {percentile:.0f}th percentile."
            )
        else:
            parts.append(f"Your raw {score_name} score is {raw_score:.4f}.")

        if coverage < 0.5:
            parts.append(
                f"Note: Only {coverage * 100:.0f}% of SNPs were available in your data. "
                "This may affect accuracy."
            )
        elif coverage < 0.8:
            parts.append(
                f"Coverage: {coverage * 100:.0f}% of SNPs found in your data."
            )

        return " ".join(parts)


# Convenience function
def calculate_prs(
    dataset: DNADataset,
    weights_file: Path | str,
    score_name: Optional[str] = None,
) -> PRSResult:
    """
    Calculate a polygenic risk score.

    Args:
        dataset: DNADataset with genotypes
        weights_file: Path to weights CSV file
        score_name: Optional name for the score

    Returns:
        PRSResult with calculated score

    Example:
        >>> from genomeinsight import load_dna_data
        >>> from genomeinsight.polygenic import calculate_prs
        >>> 
        >>> data = load_dna_data("my_dna.txt")
        >>> result = calculate_prs(data, "cad_weights.csv")
        >>> print(f"Score: {result.raw_score:.3f}, Coverage: {result.coverage_percent}")
    """
    calculator = PRSCalculator()
    return calculator.calculate_from_file(dataset, weights_file, score_name)
