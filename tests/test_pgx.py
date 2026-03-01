"""Tests for the pharmacogenomics module."""

import pytest

from pydna_analyzer.pharmacogenomics import (
    DRUG_RECOMMENDATIONS,
    GENE_DEFINITIONS,
    CPICLevel,
    MetabolizerPhenotype,
    PGxAnalyzer,
)

# =============================================================================
# TestStarAlleleDefinitions — data integrity
# =============================================================================


class TestStarAlleleDefinitions:
    """Verify the gene/star allele database is well-formed."""

    def test_all_genes_have_relevant_rsids(self):
        for gene_name, gene_def in GENE_DEFINITIONS.items():
            assert len(gene_def.relevant_rsids) > 0, f"{gene_name} has no relevant rsIDs"

    def test_star_allele_snps_match_gene_rsids(self):
        """Every SNP referenced in a star allele must be in the gene's relevant_rsids."""
        for gene_name, gene_def in GENE_DEFINITIONS.items():
            for star in gene_def.star_alleles:
                for rsid in star.defining_snps:
                    assert rsid in gene_def.relevant_rsids, (
                        f"{gene_name} {star.name}: {rsid} not in relevant_rsids"
                    )

    def test_activity_scores_are_valid(self):
        """Activity scores must be between 0.0 and 2.0."""
        for gene_name, gene_def in GENE_DEFINITIONS.items():
            for star in gene_def.star_alleles:
                assert 0.0 <= star.activity_score <= 2.0, (
                    f"{gene_name} {star.name}: invalid activity_score {star.activity_score}"
                )

    def test_single_snp_genes_have_no_star_alleles(self):
        """VKORC1 and COMT should be flagged as single-SNP with no star alleles."""
        for gene_name, gene_def in GENE_DEFINITIONS.items():
            if gene_def.is_single_snp_gene:
                assert len(gene_def.star_alleles) == 0, (
                    f"{gene_name} is single-SNP but has star alleles"
                )
                assert len(gene_def.relevant_rsids) == 1, (
                    f"{gene_name} is single-SNP but has {len(gene_def.relevant_rsids)} rsIDs"
                )


# =============================================================================
# TestDiplotypeCaller — star allele calling logic
# =============================================================================


class TestDiplotypeCaller:
    """Test diplotype calling from genotype data."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.analyzer = PGxAnalyzer()

    def test_wildtype_both_alleles(self, pgx_all_wildtype_dataset):
        """No variant alleles -> *1/*1."""
        result = self.analyzer.analyze_gene(pgx_all_wildtype_dataset, "CYP2C19")
        assert result is not None
        assert result.diplotype == "*1/*1"

    def test_heterozygous_single_variant(self, pgx_missing_snps_dataset):
        """One copy of *2 -> *1/*2."""
        result = self.analyzer.analyze_gene(pgx_missing_snps_dataset, "CYP2C19")
        assert result is not None
        assert result.diplotype == "*1/*2"

    def test_homozygous_variant(self, pgx_cyp2c19_pm_dataset):
        """Two copies of *2 -> *2/*2."""
        result = self.analyzer.analyze_gene(pgx_cyp2c19_pm_dataset, "CYP2C19")
        assert result is not None
        assert result.diplotype == "*2/*2"

    def test_compound_heterozygote(self, pgx_cyp2c9_compound_het_dataset):
        """CYP2C9 *2 het + *3 het -> *2/*3."""
        result = self.analyzer.analyze_gene(pgx_cyp2c9_compound_het_dataset, "CYP2C9")
        assert result is not None
        assert result.diplotype == "*2/*3"

    def test_cyp2c19_ultrarapid_star17(self, pgx_cyp2c19_um_dataset):
        """CYP2C19 *17/*17 ultrarapid."""
        result = self.analyzer.analyze_gene(pgx_cyp2c19_um_dataset, "CYP2C19")
        assert result is not None
        assert result.diplotype == "*17/*17"

    def test_missing_all_snps_defaults_wildtype(self, pgx_vkorc1_high_dataset):
        """When all SNPs for a multi-SNP gene are missing, default to *1/*1."""
        # This dataset only has VKORC1, so CYP2D6 SNPs are all missing
        result = self.analyzer.analyze_gene(pgx_vkorc1_high_dataset, "CYP2D6")
        assert result is not None
        assert result.diplotype == "*1/*1"
        assert result.confidence == "low"

    def test_unknown_gene_returns_none(self, pgx_all_wildtype_dataset):
        """Requesting an unknown gene returns None."""
        result = self.analyzer.analyze_gene(pgx_all_wildtype_dataset, "FAKE_GENE")
        assert result is None

    def test_partial_data_confidence(self, pgx_missing_snps_dataset):
        """Partial SNP coverage should reduce confidence."""
        result = self.analyzer.analyze_gene(pgx_missing_snps_dataset, "CYP2C19")
        assert result is not None
        assert result.confidence == "moderate"
        assert result.snps_missing > 0


# =============================================================================
# TestPhenotypePredictor — activity score to phenotype mapping
# =============================================================================


class TestPhenotypePredictor:
    """Test phenotype prediction from activity scores and genotypes."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.analyzer = PGxAnalyzer()

    def test_normal_metabolizer(self, pgx_all_wildtype_dataset):
        """*1/*1 (score 2.0) -> Normal Metabolizer."""
        result = self.analyzer.analyze_gene(pgx_all_wildtype_dataset, "CYP2C19")
        assert result.phenotype == MetabolizerPhenotype.NORMAL

    def test_intermediate_metabolizer(self, pgx_missing_snps_dataset):
        """CYP2C19 *1/*2 (score 1.0) -> Intermediate Metabolizer."""
        result = self.analyzer.analyze_gene(pgx_missing_snps_dataset, "CYP2C19")
        assert result.phenotype == MetabolizerPhenotype.INTERMEDIATE

    def test_poor_metabolizer(self, pgx_cyp2c19_pm_dataset):
        """CYP2C19 *2/*2 (score 0.0) -> Poor Metabolizer."""
        result = self.analyzer.analyze_gene(pgx_cyp2c19_pm_dataset, "CYP2C19")
        assert result.phenotype == MetabolizerPhenotype.POOR
        assert result.activity_score == 0.0

    def test_ultrarapid_metabolizer(self, pgx_cyp2c19_um_dataset):
        """CYP2C19 *17/*17 (score 3.0) -> Ultrarapid Metabolizer."""
        result = self.analyzer.analyze_gene(pgx_cyp2c19_um_dataset, "CYP2C19")
        assert result.phenotype == MetabolizerPhenotype.ULTRARAPID
        assert result.activity_score == 3.0

    def test_vkorc1_high_sensitivity(self, pgx_vkorc1_high_dataset):
        """VKORC1 AA -> High Sensitivity."""
        result = self.analyzer.analyze_gene(pgx_vkorc1_high_dataset, "VKORC1")
        assert result.phenotype == MetabolizerPhenotype.HIGH_SENSITIVITY
        assert result.diplotype == "AA"

    def test_vkorc1_normal_sensitivity(self, pgx_all_wildtype_dataset):
        """VKORC1 GG -> Normal Sensitivity."""
        result = self.analyzer.analyze_gene(pgx_all_wildtype_dataset, "VKORC1")
        assert result.phenotype == MetabolizerPhenotype.NORMAL_SENSITIVITY

    def test_comt_low_activity(self, pgx_comt_low_dataset):
        """COMT AA (Met/Met) -> Low Activity."""
        result = self.analyzer.analyze_gene(pgx_comt_low_dataset, "COMT")
        assert result.phenotype == MetabolizerPhenotype.LOW_ACTIVITY

    def test_comt_high_activity(self, pgx_all_wildtype_dataset):
        """COMT GG (Val/Val) -> High Activity."""
        result = self.analyzer.analyze_gene(pgx_all_wildtype_dataset, "COMT")
        assert result.phenotype == MetabolizerPhenotype.HIGH_ACTIVITY


# =============================================================================
# TestDrugRecommendations — correct recs for phenotype
# =============================================================================


class TestDrugRecommendations:
    """Test drug recommendation lookup."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.analyzer = PGxAnalyzer()

    def test_clopidogrel_poor_metabolizer(self, pgx_cyp2c19_pm_dataset):
        """CYP2C19 PM should get clopidogrel avoidance recommendation."""
        result = self.analyzer.analyze_gene(pgx_cyp2c19_pm_dataset, "CYP2C19")
        drug_names = [r.drug_name for r in result.drug_recommendations]
        assert "Clopidogrel" in drug_names
        clopi_rec = next(r for r in result.drug_recommendations if r.drug_name == "Clopidogrel")
        assert clopi_rec.cpic_level == CPICLevel.STRONG
        assert len(clopi_rec.alternatives) > 0

    def test_warfarin_vkorc1_sensitivity(self, pgx_vkorc1_high_dataset):
        """VKORC1 high sensitivity should recommend lower warfarin dose."""
        result = self.analyzer.analyze_gene(pgx_vkorc1_high_dataset, "VKORC1")
        drug_names = [r.drug_name for r in result.drug_recommendations]
        assert "Warfarin" in drug_names

    def test_codeine_cyp2d6_poor(self, pgx_all_wildtype_dataset, tmp_path):
        """CYP2D6 PM should get codeine avoidance recommendation."""
        # Create a dataset with CYP2D6 *4/*4
        content = (
            "#AncestryDNA raw data download\n"
            "#rsid\tchromosome\tposition\tallele1\tallele2\n"
            "rs3892097\t22\t42524947\tA\tA\n"
            "rs1065852\t22\t42522613\tG\tG\n"
        )
        from pydna_analyzer.core.data_loader import load_dna_data

        filepath = tmp_path / "cyp2d6_pm.txt"
        filepath.write_text(content)
        dataset = load_dna_data(filepath)
        result = self.analyzer.analyze_gene(dataset, "CYP2D6")
        assert result.phenotype == MetabolizerPhenotype.POOR
        drug_names = [r.drug_name for r in result.drug_recommendations]
        assert "Codeine" in drug_names

    def test_normal_metabolizer_no_actionable_recs(self, pgx_all_wildtype_dataset):
        """Normal metabolizer should have no actionable recommendations."""
        result = self.analyzer.analyze_gene(pgx_all_wildtype_dataset, "CYP2C19")
        assert result.phenotype == MetabolizerPhenotype.NORMAL
        assert not result.is_actionable


# =============================================================================
# TestPGxResult — structure and helpers
# =============================================================================


class TestPGxResult:
    """Test PGxResult data structure and computed properties."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.analyzer = PGxAnalyzer()

    def test_full_analysis_returns_all_genes(self, pgx_all_wildtype_dataset):
        """Full analysis should return results for all 6 genes."""
        result = self.analyzer.analyze(pgx_all_wildtype_dataset)
        assert result.total_genes_tested == 6
        gene_names = {r.gene for r in result.gene_results}
        assert gene_names == {"CYP2C9", "CYP2C19", "CYP2D6", "TPMT", "VKORC1", "COMT"}

    def test_all_wildtype_zero_actionable(self, pgx_all_wildtype_dataset):
        """All wild-type should have 0 actionable results."""
        result = self.analyzer.analyze(pgx_all_wildtype_dataset)
        assert result.actionable_count == 0

    def test_poor_metabolizer_actionable(self, pgx_cyp2c19_pm_dataset):
        """CYP2C19 PM should be counted as actionable."""
        result = self.analyzer.analyze(pgx_cyp2c19_pm_dataset)
        assert result.actionable_count >= 1

    def test_all_drug_recommendations_aggregated(self, pgx_cyp2c19_pm_dataset):
        """All drug recs should be aggregated across genes."""
        result = self.analyzer.analyze(pgx_cyp2c19_pm_dataset)
        all_recs = result.all_drug_recommendations
        assert len(all_recs) >= 1
        # Should include clopidogrel rec from CYP2C19 PM
        drug_names = [r.drug_name for r in all_recs]
        assert "Clopidogrel" in drug_names

    def test_summary_message(self, pgx_all_wildtype_dataset):
        """Summary should mention number of genes analyzed."""
        result = self.analyzer.analyze(pgx_all_wildtype_dataset)
        assert "6" in result.summary
        assert "pharmacogene" in result.summary.lower()

    def test_snp_counts(self, pgx_all_wildtype_dataset):
        """Total SNP counts should add up correctly."""
        result = self.analyzer.analyze(pgx_all_wildtype_dataset)
        expected_total = sum(
            len(g.relevant_rsids) for g in GENE_DEFINITIONS.values()
        )
        assert result.total_snps_tested + result.total_snps_missing == expected_total


# =============================================================================
# TestDrugRecommendationDatabase — database integrity
# =============================================================================


class TestDrugRecommendationDatabase:
    """Verify the drug recommendation database is well-formed."""

    def test_all_recs_reference_valid_genes(self):
        """Every drug recommendation must reference a gene in GENE_DEFINITIONS."""
        for rec in DRUG_RECOMMENDATIONS:
            assert rec.gene in GENE_DEFINITIONS, (
                f"Drug rec for {rec.drug_name} references unknown gene {rec.gene}"
            )

    def test_all_recs_have_nonempty_recommendation(self):
        for rec in DRUG_RECOMMENDATIONS:
            assert len(rec.recommendation) > 0

    def test_cpic_levels_are_valid(self):
        for rec in DRUG_RECOMMENDATIONS:
            assert isinstance(rec.cpic_level, CPICLevel)
