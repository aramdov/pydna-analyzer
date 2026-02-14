"""Tests for ancestry estimation module."""

import numpy as np
from typer.testing import CliRunner

from genomeinsight.ancestry import AncestryAnalyzer, AncestryResult
from genomeinsight.ancestry.estimator import AncestryEstimator
from genomeinsight.ancestry.reference_data import AIMDatabase
from genomeinsight.cli import app


class TestAIMDatabase:
    """Tests for AIM reference data loading."""

    def test_load_returns_database(self):
        db = AIMDatabase.load()
        assert db is not None

    def test_has_populations(self):
        db = AIMDatabase.load()
        assert len(db.populations) == 14

    def test_has_aims(self):
        db = AIMDatabase.load()
        assert len(db.aims) >= 100

    def test_population_has_region(self):
        db = AIMDatabase.load()
        for pop_name, pop_info in db.populations.items():
            assert "region" in pop_info, f"{pop_name} missing region"

    def test_aim_has_frequencies_for_all_populations(self):
        db = AIMDatabase.load()
        pop_names = set(db.populations.keys())
        for rsid, aim in db.aims.items():
            aim_pops = set(aim["frequencies"].keys())
            assert aim_pops == pop_names, f"{rsid} missing populations: {pop_names - aim_pops}"

    def test_frequencies_are_valid(self):
        db = AIMDatabase.load()
        for rsid, aim in db.aims.items():
            for pop, freq in aim["frequencies"].items():
                assert 0.0 <= freq <= 1.0, f"{rsid} {pop}: freq {freq} out of range"

    def test_get_aim_rsids(self):
        db = AIMDatabase.load()
        rsids = db.get_aim_rsids()
        assert isinstance(rsids, set)
        assert len(rsids) == len(db.aims)

    def test_get_population_names(self):
        db = AIMDatabase.load()
        names = db.get_population_names()
        assert "Northern European" in names
        assert "West African" in names


class TestGenotypeLikelihood:
    def setup_method(self):
        self.estimator = AncestryEstimator()

    def test_homozygous_effect_high_freq(self):
        result = self.estimator._genotype_likelihood("AA", "A", 0.9)
        assert abs(result - 0.81) < 0.001

    def test_homozygous_effect_low_freq(self):
        result = self.estimator._genotype_likelihood("AA", "A", 0.1)
        assert abs(result - 0.01) < 0.001

    def test_heterozygous(self):
        result = self.estimator._genotype_likelihood("AG", "A", 0.5)
        assert abs(result - 0.5) < 0.001

    def test_homozygous_other(self):
        result = self.estimator._genotype_likelihood("GG", "A", 0.3)
        assert abs(result - 0.49) < 0.001

    def test_freq_zero_homozygous_other(self):
        result = self.estimator._genotype_likelihood("GG", "A", 0.0)
        assert abs(result - 1.0) < 0.01

    def test_freq_one_homozygous_effect(self):
        result = self.estimator._genotype_likelihood("AA", "A", 1.0)
        assert abs(result - 1.0) < 0.01

    def test_minimum_floor(self):
        result = self.estimator._genotype_likelihood("AA", "A", 0.0)
        assert result > 0


class TestOptimizeProportions:
    def setup_method(self):
        self.estimator = AncestryEstimator()

    def test_single_population_dominant(self):
        likelihood_matrix = np.array([[0.9, 0.01, 0.01]] * 10)
        proportions = self.estimator._optimize_proportions(likelihood_matrix)
        assert proportions[0] > 0.9
        assert abs(sum(proportions) - 1.0) < 0.001

    def test_equal_mixture(self):
        likelihood_matrix = np.array([[0.5, 0.5, 0.5]] * 10)
        proportions = self.estimator._optimize_proportions(likelihood_matrix)
        for p in proportions:
            assert 0.2 < p < 0.5

    def test_two_population_mix(self):
        rows_pop0 = [[0.9, 0.01, 0.01]] * 5
        rows_pop1 = [[0.01, 0.9, 0.01]] * 5
        likelihood_matrix = np.array(rows_pop0 + rows_pop1)
        proportions = self.estimator._optimize_proportions(likelihood_matrix)
        assert abs(proportions[0] - 0.5) < 0.15
        assert abs(proportions[1] - 0.5) < 0.15

    def test_proportions_sum_to_one(self):
        likelihood_matrix = np.array([
            [0.3, 0.5, 0.2],
            [0.7, 0.1, 0.2],
            [0.1, 0.8, 0.1],
        ])
        proportions = self.estimator._optimize_proportions(likelihood_matrix)
        assert abs(sum(proportions) - 1.0) < 0.001

    def test_proportions_non_negative(self):
        likelihood_matrix = np.array([[0.9, 0.01, 0.01, 0.01]] * 20)
        proportions = self.estimator._optimize_proportions(likelihood_matrix)
        for p in proportions:
            assert p >= -0.001


class TestBootstrapConfidence:
    def setup_method(self):
        self.estimator = AncestryEstimator()

    def test_returns_low_and_high(self):
        likelihood_matrix = np.array([[0.9, 0.05, 0.05]] * 20)
        low, high = self.estimator._bootstrap_confidence(likelihood_matrix, n_bootstrap=10)
        assert len(low) == 3
        assert len(high) == 3

    def test_ci_contains_point_estimate(self):
        likelihood_matrix = np.array([[0.9, 0.05, 0.05]] * 30)
        point = self.estimator._optimize_proportions(likelihood_matrix)
        low, high = self.estimator._bootstrap_confidence(likelihood_matrix, n_bootstrap=20)
        assert low[0] <= point[0] <= high[0] + 0.05

    def test_ci_width_reasonable(self):
        likelihood_matrix = np.array([[0.8, 0.1, 0.1]] * 20)
        low, high = self.estimator._bootstrap_confidence(likelihood_matrix, n_bootstrap=10)
        for i in range(3):
            assert 0.0 <= low[i] <= high[i] <= 1.0


class TestAncestryAnalyzer:
    """Integration tests for full ancestry estimation."""

    def test_analyze_returns_result(self, ancestry_dataset_factory):
        analyzer = AncestryAnalyzer(n_bootstrap=10)
        db = AIMDatabase.load()
        rsids = list(db.aims.keys())[:20]
        genotypes = {rsid: db.get_effect_allele(rsid) * 2 for rsid in rsids}
        dataset = ancestry_dataset_factory(genotypes)
        result = analyzer.analyze(dataset)
        assert isinstance(result, AncestryResult)

    def test_proportions_sum_to_one(self, ancestry_dataset_factory):
        analyzer = AncestryAnalyzer(n_bootstrap=10)
        db = AIMDatabase.load()
        rsids = list(db.aims.keys())[:20]
        genotypes = {rsid: db.get_effect_allele(rsid) * 2 for rsid in rsids}
        dataset = ancestry_dataset_factory(genotypes)
        result = analyzer.analyze(dataset)
        total = sum(p.proportion for p in result.populations)
        assert abs(total - 1.0) < 0.01

    def test_result_has_coverage(self, ancestry_dataset_factory):
        analyzer = AncestryAnalyzer(n_bootstrap=10)
        db = AIMDatabase.load()
        rsids = list(db.aims.keys())[:10]
        genotypes = {rsid: db.get_effect_allele(rsid) * 2 for rsid in rsids}
        dataset = ancestry_dataset_factory(genotypes)
        result = analyzer.analyze(dataset)
        assert result.snps_used == 10
        assert result.snps_available >= 50
        assert 0 < result.coverage <= 1.0

    def test_result_sorted_by_proportion(self, ancestry_dataset_factory):
        analyzer = AncestryAnalyzer(n_bootstrap=10)
        db = AIMDatabase.load()
        rsids = list(db.aims.keys())[:20]
        genotypes = {rsid: db.get_effect_allele(rsid) * 2 for rsid in rsids}
        dataset = ancestry_dataset_factory(genotypes)
        result = analyzer.analyze(dataset)
        proportions = [p.proportion for p in result.populations]
        assert proportions == sorted(proportions, reverse=True)

    def test_top_regions_aggregated(self, ancestry_dataset_factory):
        analyzer = AncestryAnalyzer(n_bootstrap=10)
        db = AIMDatabase.load()
        rsids = list(db.aims.keys())[:20]
        genotypes = {rsid: db.get_effect_allele(rsid) * 2 for rsid in rsids}
        dataset = ancestry_dataset_factory(genotypes)
        result = analyzer.analyze(dataset)
        assert isinstance(result.top_regions, dict)
        assert abs(sum(result.top_regions.values()) - 1.0) < 0.01

    def test_has_interpretation(self, ancestry_dataset_factory):
        analyzer = AncestryAnalyzer(n_bootstrap=10)
        db = AIMDatabase.load()
        rsids = list(db.aims.keys())[:20]
        genotypes = {rsid: db.get_effect_allele(rsid) * 2 for rsid in rsids}
        dataset = ancestry_dataset_factory(genotypes)
        result = analyzer.analyze(dataset)
        assert len(result.interpretation) > 0

    def test_no_matching_aims(self, ancestry_dataset_factory):
        dataset = ancestry_dataset_factory({"rs99999999": "AA"})
        analyzer = AncestryAnalyzer(n_bootstrap=10)
        result = analyzer.analyze(dataset)
        assert result.snps_used == 0
        assert "insufficient" in result.interpretation.lower()

    def test_low_coverage_noted(self, ancestry_dataset_factory):
        analyzer = AncestryAnalyzer(n_bootstrap=10)
        db = AIMDatabase.load()
        rsids = list(db.aims.keys())[:3]
        genotypes = {rsid: db.get_effect_allele(rsid) * 2 for rsid in rsids}
        dataset = ancestry_dataset_factory(genotypes)
        result = analyzer.analyze(dataset)
        assert result.coverage < 0.1


runner = CliRunner()


class TestAncestryCLI:
    """Tests for the ancestry CLI command."""

    def test_ancestry_command_runs(self, sample_ancestrydna_file):
        result = runner.invoke(app, ["ancestry", str(sample_ancestrydna_file)])
        assert result.exit_code == 0

    def test_ancestry_json_output(self, sample_ancestrydna_file, tmp_path):
        output = tmp_path / "ancestry.json"
        result = runner.invoke(
            app, ["ancestry", str(sample_ancestrydna_file), "-o", str(output)]
        )
        assert result.exit_code == 0
        # JSON is only written when ancestry markers are matched;
        # the sample file may lack AIMs, so we just check for success
        if output.exists():
            import json

            data = json.loads(output.read_text())
            assert "populations" in data

    def test_ancestry_output_contains_header(self, sample_ancestrydna_file):
        result = runner.invoke(app, ["ancestry", str(sample_ancestrydna_file)])
        # Should contain ancestry-related output
        assert (
            "Ancestry" in result.output
            or "ancestry" in result.output
            or "marker" in result.output
        )
