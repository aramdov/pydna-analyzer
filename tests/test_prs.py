"""Tests for Polygenic Risk Score calculator."""

import pytest
from pathlib import Path

from pydna_analyzer.core.data_loader import load_dna_data
from pydna_analyzer.polygenic import (
    PRSCalculator,
    PRSResult,
    PRSWeightLoader,
    SNPWeight,
    calculate_prs,
)


class TestSNPWeight:
    """Tests for SNPWeight data class."""

    def test_snp_weight_creation(self):
        """Should create SNPWeight with required fields."""
        weight = SNPWeight(
            rsid="rs123456",
            effect_allele="A",
            weight=0.15,
        )
        assert weight.rsid == "rs123456"
        assert weight.effect_allele == "A"
        assert weight.weight == 0.15


class TestPRSWeightLoader:
    """Tests for weight file loading."""

    def test_load_csv(self, tmp_path):
        """Should load weights from CSV file."""
        csv_content = """rsid,effect_allele,weight
rs123456,A,0.15
rs789012,G,-0.10
rs345678,C,0.25
"""
        csv_file = tmp_path / "weights.csv"
        csv_file.write_text(csv_content)

        weights = PRSWeightLoader.load_csv(csv_file)

        assert len(weights) == 3
        assert weights[0].rsid == "rs123456"
        assert weights[0].effect_allele == "A"
        assert weights[0].weight == 0.15
        assert weights[1].weight == -0.10

    def test_load_csv_alternative_columns(self, tmp_path):
        """Should handle alternative column names."""
        csv_content = """snp,allele,beta
rs123456,A,0.15
rs789012,G,-0.10
"""
        csv_file = tmp_path / "weights.csv"
        csv_file.write_text(csv_content)

        weights = PRSWeightLoader.load_csv(csv_file)

        assert len(weights) == 2
        assert weights[0].rsid == "rs123456"


class TestPRSCalculator:
    """Tests for PRS calculation."""

    @pytest.fixture
    def sample_weights(self):
        """Sample weights for testing."""
        return [
            SNPWeight(rsid="rs1801133", effect_allele="T", weight=0.20),
            SNPWeight(rsid="rs1801131", effect_allele="C", weight=0.15),
            SNPWeight(rsid="rs4680", effect_allele="A", weight=0.10),
            SNPWeight(rsid="rs999999", effect_allele="G", weight=0.05),  # Missing
        ]

    def test_calculate_returns_result(self, sample_ancestrydna_file, sample_weights):
        """Should return PRSResult."""
        dataset = load_dna_data(sample_ancestrydna_file)
        calculator = PRSCalculator()
        result = calculator.calculate(dataset, sample_weights, "Test Score")

        assert isinstance(result, PRSResult)
        assert result.score_name == "Test Score"

    def test_dosage_calculation(self, sample_ancestrydna_file, sample_weights):
        """Should correctly calculate dosage from genotypes."""
        dataset = load_dna_data(sample_ancestrydna_file)
        calculator = PRSCalculator()
        result = calculator.calculate(dataset, sample_weights, "Test")

        # Should use 3 SNPs (rs999999 is missing)
        assert result.snps_used == 3
        assert result.snps_available == 4
        assert "rs999999" in result.missing_snps

    def test_coverage_calculation(self, sample_ancestrydna_file, sample_weights):
        """Should calculate coverage correctly."""
        dataset = load_dna_data(sample_ancestrydna_file)
        calculator = PRSCalculator()
        result = calculator.calculate(dataset, sample_weights, "Test")

        # 3/4 = 75%
        assert result.coverage == 0.75
        assert result.coverage_percent == "75.0%"

    def test_raw_score_calculation(self, tmp_path):
        """Should calculate raw score as weighted sum of dosages."""
        # Create DNA data with known genotypes
        content = """#AncestryDNA raw data download
#rsid	chromosome	position	allele1	allele2
rs111111	1	1000	A	A
rs222222	1	2000	A	G
rs333333	1	3000	G	G
"""
        file_path = tmp_path / "test.txt"
        file_path.write_text(content)

        weights = [
            SNPWeight(rsid="rs111111", effect_allele="A", weight=0.10),  # dosage=2
            SNPWeight(rsid="rs222222", effect_allele="A", weight=0.20),  # dosage=1
            SNPWeight(rsid="rs333333", effect_allele="A", weight=0.30),  # dosage=0
        ]

        dataset = load_dna_data(file_path)
        calculator = PRSCalculator()
        result = calculator.calculate(dataset, weights, "Test")

        # Score = (2 * 0.10) + (1 * 0.20) + (0 * 0.30) = 0.40
        assert result.raw_score == pytest.approx(0.40)

    def test_calculate_from_file(self, sample_ancestrydna_file, tmp_path):
        """Should load weights from file and calculate."""
        weights_content = """rsid,effect_allele,weight
rs1801133,T,0.20
rs4680,A,0.10
"""
        weights_file = tmp_path / "weights.csv"
        weights_file.write_text(weights_content)

        dataset = load_dna_data(sample_ancestrydna_file)
        calculator = PRSCalculator()
        result = calculator.calculate_from_file(dataset, weights_file)

        assert result.snps_used > 0
        assert result.raw_score != 0


class TestCalculatePRSConvenience:
    """Tests for the convenience function."""

    def test_calculate_prs_function(self, sample_ancestrydna_file, tmp_path):
        """Should work as a simple convenience function."""
        weights_content = """rsid,effect_allele,weight
rs1801133,T,0.20
"""
        weights_file = tmp_path / "weights.csv"
        weights_file.write_text(weights_content)

        dataset = load_dna_data(sample_ancestrydna_file)
        result = calculate_prs(dataset, weights_file, "Simple Test")

        assert isinstance(result, PRSResult)
        assert result.score_name == "Simple Test"
