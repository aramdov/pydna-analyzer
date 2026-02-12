"""Tests for the reports module (HTML + JSON export)."""

import json

import pytest

from genomeinsight.clinical.analyzer import AnalysisResult, VariantResult
from genomeinsight.clinical.apoe import APOEAllele, APOEResult
from genomeinsight.clinical.variants import Category, EvidenceLevel, RiskLevel
from genomeinsight.reports.html_report import generate_html_report
from genomeinsight.reports.json_export import export_to_json


@pytest.fixture
def sample_apoe_result():
    """Build a realistic APOEResult for e3/e4."""
    return APOEResult(
        genotype="e3/e4",
        allele1=APOEAllele.E3,
        allele2=APOEAllele.E4,
        rs429358_genotype="TC",
        rs7412_genotype="CC",
        risk_category="elevated",
        interpretation=(
            "One e4 allele increases Alzheimer's risk ~2-3x compared to e3/e3. "
            "Also elevated cardiovascular risk."
        ),
        recommendations=[
            "Regular cognitive screening after age 50",
            "Mediterranean/MIND diet recommended",
        ],
    )


@pytest.fixture
def sample_variant_high():
    """A high-risk cardiovascular variant (MTHFR C677T TT)."""
    return VariantResult(
        rsid="rs1801133",
        gene="MTHFR",
        name="MTHFR C677T",
        category=Category.CARDIOVASCULAR,
        genotype="TT",
        matched_genotype="TT",
        risk_level=RiskLevel.ELEVATED,
        risk_allele="T",
        description="~70% reduced enzyme activity, elevated homocysteine risk",
        effect="Reduced enzyme activity affecting folate metabolism and homocysteine levels",
        evidence=EvidenceLevel.STRONG,
        interventions=["Methylfolate supplementation", "B12 monitoring"],
        chromosome="1",
        position=11856378,
        has_interactions=True,
        interaction_notes=["Compound heterozygosity with A1298C reduces MTHFR function"],
    )


@pytest.fixture
def sample_variant_moderate():
    """A moderate-risk pharmacogenomics variant (CYP2C19*2 GA)."""
    return VariantResult(
        rsid="rs4244285",
        gene="CYP2C19",
        name="CYP2C19*2",
        category=Category.PHARMACOGENOMICS,
        genotype="GA",
        matched_genotype="GA",
        risk_level=RiskLevel.MODERATE,
        risk_allele="A",
        description="Intermediate metabolizer - reduced clopidogrel effect",
        effect="Loss of function - affects clopidogrel (Plavix), PPIs, antidepressants",
        evidence=EvidenceLevel.STRONG,
        interventions=["Alternative to clopidogrel (prasugrel, ticagrelor)"],
        chromosome="10",
        position=96521657,
    )


@pytest.fixture
def sample_variant_low():
    """A low/normal-risk nutrient variant (HFE H63D CC)."""
    return VariantResult(
        rsid="rs1799945",
        gene="HFE",
        name="HFE H63D",
        category=Category.NUTRIENT,
        genotype="CC",
        matched_genotype="CC",
        risk_level=RiskLevel.NORMAL,
        risk_allele="G",
        description="Normal iron absorption",
        effect="Affects iron absorption - hereditary hemochromatosis risk",
        evidence=EvidenceLevel.STRONG,
        interventions=["Regular ferritin monitoring"],
        chromosome="6",
        position=26091179,
    )


@pytest.fixture
def sample_gene_interaction():
    """A sample gene interaction dict."""
    return {
        "type": "compound_heterozygosity",
        "genes": ["MTHFR"],
        "rsids": ["rs1801133", "rs1801131"],
        "severity": "elevated",
        "note": (
            "Compound heterozygosity detected: MTHFR C677T (TT) + "
            "MTHFR A1298C (AC). This combination significantly reduces "
            "MTHFR enzyme activity."
        ),
        "recommendations": [
            "Methylfolate (L-methylfolate) supplementation",
            "Monitor homocysteine levels",
        ],
    }


@pytest.fixture
def sample_analysis_result(
    sample_variant_high,
    sample_variant_moderate,
    sample_variant_low,
    sample_apoe_result,
    sample_gene_interaction,
):
    """Build a realistic AnalysisResult with mixed priorities, APOE, and interactions."""
    all_variants = [sample_variant_high, sample_variant_moderate, sample_variant_low]
    return AnalysisResult(
        source_file="test_dna_data.txt",
        snp_count=650000,
        variants_analyzed=35,
        variants_found=3,
        apoe_status=sample_apoe_result,
        variant_results=all_variants,
        high_priority=[sample_variant_high],
        moderate_priority=[sample_variant_moderate],
        low_priority=[sample_variant_low],
        by_category={
            Category.CARDIOVASCULAR: [sample_variant_high],
            Category.PHARMACOGENOMICS: [sample_variant_moderate],
            Category.NUTRIENT: [sample_variant_low],
        },
        gene_interactions=[sample_gene_interaction],
    )


@pytest.fixture
def empty_analysis_result():
    """An AnalysisResult with no variants found."""
    return AnalysisResult(
        source_file="empty_test.txt",
        snp_count=500000,
        variants_analyzed=35,
        variants_found=0,
        apoe_status=None,
        variant_results=[],
        high_priority=[],
        moderate_priority=[],
        low_priority=[],
        by_category={},
        gene_interactions=[],
    )


# =============================================================================
# HTML Report Tests
# =============================================================================


class TestGenerateHtmlReport:
    """Tests for generate_html_report()."""

    def test_creates_html_file(self, tmp_path, sample_analysis_result):
        """Should create an HTML file at the given path."""
        output = tmp_path / "report.html"
        generate_html_report(sample_analysis_result, output)

        assert output.exists()
        assert output.stat().st_size > 0

    def test_html_contains_structural_elements(self, tmp_path, sample_analysis_result):
        """HTML should contain key structural elements: title, disclaimer, summary stats."""
        output = tmp_path / "report.html"
        generate_html_report(sample_analysis_result, output)
        html = output.read_text()

        assert "<title>GenomeInsight Report</title>" in html
        assert "Important Disclaimer" in html
        # Summary stat cards
        assert "SNPs Analyzed" in html
        assert "Clinical Variants Found" in html
        assert "High Priority" in html
        assert "Moderate Priority" in html

    def test_variant_data_in_output(self, tmp_path, sample_analysis_result):
        """Variant data should appear in the HTML output."""
        output = tmp_path / "report.html"
        generate_html_report(sample_analysis_result, output)
        html = output.read_text()

        # Gene names from our fixture variants
        assert "MTHFR" in html
        assert "CYP2C19" in html
        assert "HFE" in html
        # rsIDs
        assert "rs1801133" in html
        assert "rs4244285" in html
        assert "rs1799945" in html

    def test_snp_count_formatted(self, tmp_path, sample_analysis_result):
        """SNP count should be formatted with comma separators."""
        output = tmp_path / "report.html"
        generate_html_report(sample_analysis_result, output)
        html = output.read_text()

        assert "650,000" in html

    def test_works_with_empty_variant_results(self, tmp_path, empty_analysis_result):
        """Should generate valid HTML even with no variants."""
        output = tmp_path / "report.html"
        generate_html_report(empty_analysis_result, output)
        html = output.read_text()

        assert output.exists()
        assert "<title>GenomeInsight Report</title>" in html
        # The "no high priority findings" placeholder should appear
        assert "No high priority findings" in html

    def test_works_without_apoe_status(self, tmp_path, sample_analysis_result):
        """Should work when apoe_status is None (APOE section omitted)."""
        sample_analysis_result.apoe_status = None
        output = tmp_path / "report.html"
        generate_html_report(sample_analysis_result, output)
        html = output.read_text()

        assert output.exists()
        # The APOE genotype string should NOT appear
        assert "e3/e4" not in html

    def test_apoe_section_present_when_set(self, tmp_path, sample_analysis_result):
        """APOE section should render when apoe_status is provided."""
        output = tmp_path / "report.html"
        generate_html_report(sample_analysis_result, output)
        html = output.read_text()

        assert "APOE Status" in html
        assert "e3/e4" in html
        assert "elevated" in html.lower()

    def test_gene_interactions_present(self, tmp_path, sample_analysis_result):
        """Gene interaction alerts should appear in the HTML."""
        output = tmp_path / "report.html"
        generate_html_report(sample_analysis_result, output)
        html = output.read_text()

        assert "Gene Interactions Detected" in html
        assert "Compound heterozygosity detected" in html

    def test_no_gene_interactions_section_when_empty(self, tmp_path, empty_analysis_result):
        """Gene interactions section should not appear when there are none."""
        output = tmp_path / "report.html"
        generate_html_report(empty_analysis_result, output)
        html = output.read_text()

        assert "Gene Interactions Detected" not in html

    def test_risk_level_badges(self, tmp_path, sample_analysis_result):
        """Risk level values should appear as badge classes in the HTML."""
        output = tmp_path / "report.html"
        generate_html_report(sample_analysis_result, output)
        html = output.read_text()

        assert "risk-elevated" in html
        assert "risk-moderate" in html
        assert "risk-normal" in html

    def test_version_in_footer(self, tmp_path, sample_analysis_result):
        """Footer should contain the GenomeInsight version."""
        output = tmp_path / "report.html"
        generate_html_report(sample_analysis_result, output)
        html = output.read_text()

        from genomeinsight import __version__

        assert f"GenomeInsight v{__version__}" in html


# =============================================================================
# JSON Export Tests
# =============================================================================


class TestExportToJson:
    """Tests for export_to_json()."""

    def test_creates_json_file(self, tmp_path, sample_analysis_result):
        """Should create a valid JSON file at the given path."""
        output = tmp_path / "report.json"
        export_to_json(sample_analysis_result, output)

        assert output.exists()
        # Should be parseable JSON
        data = json.loads(output.read_text())
        assert isinstance(data, dict)

    def test_top_level_keys(self, tmp_path, sample_analysis_result):
        """JSON should have expected top-level keys."""
        output = tmp_path / "report.json"
        export_to_json(sample_analysis_result, output)
        data = json.loads(output.read_text())

        expected_keys = {
            "metadata",
            "summary",
            "apoe_status",
            "variants",
            "by_category",
            "gene_interactions",
        }
        assert expected_keys.issubset(set(data.keys()))

    def test_metadata_fields(self, tmp_path, sample_analysis_result):
        """Metadata should contain version, generated_at, and source_file."""
        output = tmp_path / "report.json"
        export_to_json(sample_analysis_result, output)
        data = json.loads(output.read_text())

        meta = data["metadata"]
        assert "version" in meta
        assert "generated_at" in meta
        assert meta["source_file"] == "test_dna_data.txt"

    def test_variant_serialization(self, tmp_path, sample_analysis_result):
        """Variant data should be correctly serialized (enums to strings)."""
        output = tmp_path / "report.json"
        export_to_json(sample_analysis_result, output)
        data = json.loads(output.read_text())

        high = data["variants"]["high_priority"]
        assert len(high) == 1
        variant = high[0]
        assert variant["rsid"] == "rs1801133"
        assert variant["gene"] == "MTHFR"
        assert variant["risk_level"] == "elevated"
        assert variant["category"] == "cardiovascular"
        assert variant["evidence"] == "strong"
        assert isinstance(variant["interventions"], list)

    def test_summary_counts_match(self, tmp_path, sample_analysis_result):
        """Summary counts should match the input data."""
        output = tmp_path / "report.json"
        export_to_json(sample_analysis_result, output)
        data = json.loads(output.read_text())

        summary = data["summary"]
        assert summary["snp_count"] == 650000
        assert summary["variants_analyzed"] == 35
        assert summary["variants_found"] == 3
        assert summary["high_priority_count"] == 1
        assert summary["moderate_priority_count"] == 1
        assert summary["low_priority_count"] == 1
        assert summary["gene_interactions_count"] == 1

    def test_apoe_status_serialized_when_present(self, tmp_path, sample_analysis_result):
        """APOE status should be correctly serialized when present."""
        output = tmp_path / "report.json"
        export_to_json(sample_analysis_result, output)
        data = json.loads(output.read_text())

        apoe = data["apoe_status"]
        assert apoe is not None
        assert apoe["genotype"] == "e3/e4"
        assert apoe["allele1"] == "e3"
        assert apoe["allele2"] == "e4"
        assert apoe["rs429358_genotype"] == "TC"
        assert apoe["rs7412_genotype"] == "CC"
        assert apoe["risk_category"] == "elevated"
        assert isinstance(apoe["recommendations"], list)
        assert len(apoe["recommendations"]) == 2
        assert apoe["has_e4"] is True
        assert apoe["is_e4_homozygous"] is False
        assert apoe["has_e2"] is False

    def test_apoe_status_null_when_absent(self, tmp_path, empty_analysis_result):
        """APOE status should be null when not present."""
        output = tmp_path / "report.json"
        export_to_json(empty_analysis_result, output)
        data = json.loads(output.read_text())

        assert data["apoe_status"] is None

    def test_gene_interactions_included(self, tmp_path, sample_analysis_result):
        """Gene interactions should be included in the JSON."""
        output = tmp_path / "report.json"
        export_to_json(sample_analysis_result, output)
        data = json.loads(output.read_text())

        interactions = data["gene_interactions"]
        assert len(interactions) == 1
        interaction = interactions[0]
        assert interaction["type"] == "compound_heterozygosity"
        assert interaction["genes"] == ["MTHFR"]
        assert "rs1801133" in interaction["rsids"]
        assert "rs1801131" in interaction["rsids"]
        assert isinstance(interaction["recommendations"], list)

    def test_by_category_keys(self, tmp_path, sample_analysis_result):
        """by_category should use category string values as keys."""
        output = tmp_path / "report.json"
        export_to_json(sample_analysis_result, output)
        data = json.loads(output.read_text())

        by_cat = data["by_category"]
        assert "cardiovascular" in by_cat
        assert "pharmacogenomics" in by_cat
        assert "nutrient" in by_cat
        assert len(by_cat["cardiovascular"]) == 1
        assert by_cat["cardiovascular"][0]["gene"] == "MTHFR"

    def test_empty_variants_produces_empty_lists(self, tmp_path, empty_analysis_result):
        """Empty variant results should produce empty lists in JSON."""
        output = tmp_path / "report.json"
        export_to_json(empty_analysis_result, output)
        data = json.loads(output.read_text())

        assert data["variants"]["high_priority"] == []
        assert data["variants"]["moderate_priority"] == []
        assert data["variants"]["low_priority"] == []
        assert data["by_category"] == {}
        assert data["gene_interactions"] == []

    def test_variant_interaction_fields(self, tmp_path, sample_analysis_result):
        """Variant serialization should include has_interactions and interaction_notes."""
        output = tmp_path / "report.json"
        export_to_json(sample_analysis_result, output)
        data = json.loads(output.read_text())

        high_variant = data["variants"]["high_priority"][0]
        assert high_variant["has_interactions"] is True
        assert isinstance(high_variant["interaction_notes"], list)
        assert len(high_variant["interaction_notes"]) == 1

        moderate_variant = data["variants"]["moderate_priority"][0]
        assert moderate_variant["has_interactions"] is False
        assert moderate_variant["interaction_notes"] == []
