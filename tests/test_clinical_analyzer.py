"""Tests for clinical variant analyzer."""

import pytest

from pydna_analyzer.core.data_loader import load_dna_data
from pydna_analyzer.clinical.analyzer import ClinicalAnalyzer, AnalysisResult
from pydna_analyzer.clinical.variants import RiskLevel, Category


class TestClinicalAnalyzer:
    """Tests for the ClinicalAnalyzer."""

    def test_analyze_returns_result(self, sample_ancestrydna_file):
        """Should return an AnalysisResult."""
        dataset = load_dna_data(sample_ancestrydna_file)
        analyzer = ClinicalAnalyzer()
        result = analyzer.analyze(dataset)

        assert isinstance(result, AnalysisResult)
        assert result.snp_count == dataset.snp_count

    def test_variants_found(self, sample_ancestrydna_file):
        """Should find matching clinical variants."""
        dataset = load_dna_data(sample_ancestrydna_file)
        analyzer = ClinicalAnalyzer()
        result = analyzer.analyze(dataset)

        # Sample data contains known clinical variants
        assert result.variants_found > 0
        assert len(result.variant_results) > 0

    def test_mthfr_detection(self, sample_ancestrydna_file):
        """Should detect MTHFR C677T variant."""
        dataset = load_dna_data(sample_ancestrydna_file)
        analyzer = ClinicalAnalyzer()
        result = analyzer.analyze(dataset)

        # Find rs1801133 (MTHFR C677T)
        mthfr = next(
            (v for v in result.variant_results if v.rsid == "rs1801133"),
            None
        )
        assert mthfr is not None
        assert mthfr.gene == "MTHFR"
        assert mthfr.genotype == "CT"  # Heterozygous in sample

    def test_priority_classification(self, sample_ancestrydna_file):
        """Should classify variants by priority."""
        dataset = load_dna_data(sample_ancestrydna_file)
        analyzer = ClinicalAnalyzer()
        result = analyzer.analyze(dataset)

        # All variants should be in one priority list
        total = (
            len(result.high_priority) +
            len(result.moderate_priority) +
            len(result.low_priority)
        )
        assert total == len(result.variant_results)

    def test_category_grouping(self, sample_ancestrydna_file):
        """Should group variants by category."""
        dataset = load_dna_data(sample_ancestrydna_file)
        analyzer = ClinicalAnalyzer()
        result = analyzer.analyze(dataset)

        # Check that categories are populated
        assert isinstance(result.by_category, dict)

        # All variants should be in a category
        categorized = sum(len(v) for v in result.by_category.values())
        assert categorized == len(result.variant_results)

    def test_reverse_complement_strand_genotypes_are_recognized(self, tmp_path):
        """Consumer reverse-strand calls should map to the curated genotype definitions."""
        content = """#AncestryDNA raw data download
#rsid\tchromosome\tposition\tallele1\tallele2
rs1801133\t1\t11856378\tG\tA
rs6025\t1\t169519049\tC\tC
rs9923231\t16\t31107689\tC\tC
rs1800497\t11\t113400106\tA\tG
rs1800566\t16\t69745134\tA\tG
rs1801260\t4\t56357331\tA\tA
rs1048943\t15\t75011876\tT\tT
"""
        file_path = tmp_path / "reverse_strand_calls.txt"
        file_path.write_text(content)

        dataset = load_dna_data(file_path)
        analyzer = ClinicalAnalyzer()
        result = analyzer.analyze(dataset)
        by_rsid = {variant.rsid: variant for variant in result.variant_results}

        assert by_rsid["rs1801133"].matched_genotype in {"CT", "TC"}
        assert by_rsid["rs1801133"].risk_level == RiskLevel.MODERATE
        assert by_rsid["rs6025"].matched_genotype == "GG"
        assert by_rsid["rs6025"].risk_level == RiskLevel.NORMAL
        assert by_rsid["rs9923231"].matched_genotype == "GG"
        assert by_rsid["rs9923231"].risk_level == RiskLevel.NORMAL
        assert by_rsid["rs1800497"].matched_genotype in {"CT", "TC"}
        assert by_rsid["rs1800497"].risk_level == RiskLevel.MODERATE
        assert by_rsid["rs1800566"].matched_genotype in {"CT", "TC"}
        assert by_rsid["rs1800566"].risk_level == RiskLevel.MODERATE
        assert by_rsid["rs1801260"].matched_genotype == "TT"
        assert by_rsid["rs1801260"].risk_level == RiskLevel.NORMAL
        assert by_rsid["rs1048943"].matched_genotype == "AA"
        assert by_rsid["rs1048943"].risk_level == RiskLevel.NORMAL


class TestVariantResult:
    """Tests for individual variant results."""

    def test_variant_has_required_fields(self, sample_ancestrydna_file):
        """Each variant should have all required fields."""
        dataset = load_dna_data(sample_ancestrydna_file)
        analyzer = ClinicalAnalyzer()
        result = analyzer.analyze(dataset)

        for variant in result.variant_results:
            assert variant.rsid is not None
            assert variant.gene is not None
            assert variant.genotype is not None
            assert isinstance(variant.risk_level, RiskLevel)
            assert isinstance(variant.category, Category)
            assert variant.description is not None


class TestGeneInteractions:
    """Tests for gene-gene interaction detection."""

    def test_mthfr_compound_heterozygosity(self, tmp_path):
        """Should detect MTHFR compound heterozygosity."""
        # Create data with both MTHFR variants heterozygous
        content = """#AncestryDNA raw data download
#rsid	chromosome	position	allele1	allele2
rs1801133	1	11856378	C	T
rs1801131	1	11854476	A	C
"""
        file_path = tmp_path / "mthfr_compound.txt"
        file_path.write_text(content)

        dataset = load_dna_data(file_path)
        analyzer = ClinicalAnalyzer()
        result = analyzer.analyze(dataset)

        # Should detect the interaction
        assert len(result.gene_interactions) > 0

        # Find the MTHFR interaction
        mthfr_interaction = next(
            (i for i in result.gene_interactions if "MTHFR" in str(i)),
            None
        )
        assert mthfr_interaction is not None
