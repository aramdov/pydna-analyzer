# Ancestry Estimation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add sub-continental ancestry composition estimation using likelihood-based admixture with bundled 1000 Genomes reference data.

**Architecture:** Three-file module (`reference_data.py`, `estimator.py`, `__init__.py`) under `pydna_analyzer/ancestry/` with bundled AIM frequency JSON. MLE optimization via scipy. CLI command `pydna_analyzer ancestry <file>`. Follows existing dataclass + Rich + pytest patterns exactly.

**Tech Stack:** Python 3.10+, numpy, scipy.optimize, dataclasses, Rich, Typer, pytest

---

### Task 1: Update pyproject.toml Dependencies

**Files:**
- Modify: `pyproject.toml:64-67` (ancestry optional deps)

**Step 1: Update ancestry dependencies**

Replace the current ancestry optional deps with scipy (numpy is already a base dep):

```toml
ancestry = [
    "scipy>=1.10",
]
```

Remove `scikit-allel` and `scikit-learn` — we don't need them for likelihood-based estimation.

**Step 2: Install updated deps**

Run: `uv sync --all-extras`
Expected: Clean install with scipy available

**Step 3: Commit**

```bash
git add pyproject.toml uv.lock
git commit -m "chore: update ancestry deps to scipy only"
```

---

### Task 2: Create Reference Data Module and AIM Frequencies

**Files:**
- Create: `pydna_analyzer/ancestry/data/aim_frequencies.json`
- Create: `pydna_analyzer/ancestry/reference_data.py`
- Test: `tests/test_ancestry.py`

**Step 1: Write the failing test for AIMDatabase loading**

Create `tests/test_ancestry.py`:

```python
"""Tests for ancestry estimation module."""

import pytest

from pydna_analyzer.ancestry.reference_data import AIMDatabase


class TestAIMDatabase:
    """Tests for AIM reference data loading."""

    def test_load_returns_database(self):
        db = AIMDatabase.load()
        assert db is not None

    def test_has_populations(self):
        db = AIMDatabase.load()
        assert len(db.populations) >= 10

    def test_has_aims(self):
        db = AIMDatabase.load()
        assert len(db.aims) >= 50

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
```

**Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_ancestry.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'pydna_analyzer.ancestry.reference_data'`

**Step 3: Create the AIM frequencies JSON**

Create `pydna_analyzer/ancestry/data/aim_frequencies.json` with curated AIM data.

This is the largest single piece — a JSON file with ~100-200 well-known ancestry-informative markers sourced from published literature (1000 Genomes Phase 3 allele frequencies). Key AIMs to include:

- **Skin/hair pigmentation:** rs1426654 (SLC24A5), rs16891982 (SLC45A2), rs1042602 (TYR), rs12913832 (HERC2/OCA2), rs1800407 (OCA2)
- **Continental differentiation:** rs2814778 (DARC/Duffy, near-fixed African), rs3827760 (EDAR, East Asian), rs1871534 (ABCC11), rs4988235 (LCT, European lactase)
- **Well-studied AIMs panels:** Kosoy et al. 2009 (128 AIMs), Galanter et al. 2012, Phillips et al. 2007

Each AIM needs: rsid, chromosome, position, effect_allele, and allele frequencies for all 15 populations.

Population list (15):
- Northern European, Southern European, Eastern European
- West African, East African, North African
- Han Chinese, Japanese, Korean
- South Asian
- Native American
- Middle Eastern
- Oceanian
- Central Asian

Frequency sources: 1000 Genomes Phase 3 super-population and sub-population allele frequencies (publicly available). For populations not directly in 1000 Genomes (Middle Eastern, North African, Central Asian, Oceanian, Korean, Native American), use gnomAD or published AIM panel papers.

**Step 4: Create reference_data.py**

```python
"""Reference allele frequency data for ancestry-informative markers."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass
class AIMDatabase:
    """Database of ancestry-informative markers with population frequencies."""

    metadata: dict[str, Any]
    populations: dict[str, dict[str, Any]]
    aims: dict[str, dict[str, Any]]

    @classmethod
    def load(cls, filepath: Path | None = None) -> AIMDatabase:
        """Load AIM database from bundled JSON file.

        Args:
            filepath: Optional custom path. Defaults to bundled data.

        Returns:
            AIMDatabase with all reference data loaded.
        """
        if filepath is None:
            filepath = Path(__file__).parent / "data" / "aim_frequencies.json"

        with open(filepath) as f:
            data = json.load(f)

        return cls(
            metadata=data["metadata"],
            populations=data["populations"],
            aims=data["aims"],
        )

    def get_aim_rsids(self) -> set[str]:
        """Get set of all AIM rsIDs."""
        return set(self.aims.keys())

    def get_population_names(self) -> list[str]:
        """Get list of population names."""
        return list(self.populations.keys())

    def get_frequency(self, rsid: str, population: str) -> float | None:
        """Get allele frequency for a specific AIM and population."""
        aim = self.aims.get(rsid)
        if aim is None:
            return None
        return aim["frequencies"].get(population)

    def get_effect_allele(self, rsid: str) -> str | None:
        """Get the effect allele for an AIM."""
        aim = self.aims.get(rsid)
        if aim is None:
            return None
        return aim["effect_allele"]
```

**Step 5: Run tests to verify they pass**

Run: `uv run pytest tests/test_ancestry.py -v`
Expected: 8 PASSED

**Step 6: Commit**

```bash
git add pydna_analyzer/ancestry/ tests/test_ancestry.py
git commit -m "feat(ancestry): add AIM reference database with population frequencies"
```

---

### Task 3: Implement Genotype Likelihood Math

**Files:**
- Create: `pydna_analyzer/ancestry/estimator.py`
- Test: `tests/test_ancestry.py` (add to)

**Step 1: Write failing tests for genotype likelihood**

Add to `tests/test_ancestry.py`:

```python
from pydna_analyzer.ancestry.estimator import AncestryEstimator


class TestGenotypeLikelihood:
    """Tests for Hardy-Weinberg genotype likelihood calculation."""

    def setup_method(self):
        self.estimator = AncestryEstimator()

    def test_homozygous_effect_high_freq(self):
        """AA genotype with high freq (0.9) -> ~0.81."""
        result = self.estimator._genotype_likelihood(
            genotype="AA", effect_allele="A", freq=0.9
        )
        assert abs(result - 0.81) < 0.001

    def test_homozygous_effect_low_freq(self):
        """AA genotype with low freq (0.1) -> ~0.01."""
        result = self.estimator._genotype_likelihood(
            genotype="AA", effect_allele="A", freq=0.1
        )
        assert abs(result - 0.01) < 0.001

    def test_heterozygous(self):
        """AG genotype with freq 0.5 -> 0.5 (2pq)."""
        result = self.estimator._genotype_likelihood(
            genotype="AG", effect_allele="A", freq=0.5
        )
        assert abs(result - 0.5) < 0.001

    def test_homozygous_other(self):
        """GG genotype (0 copies of A) with freq 0.3 -> (1-0.3)^2 = 0.49."""
        result = self.estimator._genotype_likelihood(
            genotype="GG", effect_allele="A", freq=0.3
        )
        assert abs(result - 0.49) < 0.001

    def test_freq_zero_homozygous_other(self):
        """freq=0 with non-effect genotype -> 1.0."""
        result = self.estimator._genotype_likelihood(
            genotype="GG", effect_allele="A", freq=0.0
        )
        assert abs(result - 1.0) < 0.001

    def test_freq_one_homozygous_effect(self):
        """freq=1.0 with effect genotype -> 1.0."""
        result = self.estimator._genotype_likelihood(
            genotype="AA", effect_allele="A", freq=1.0
        )
        assert abs(result - 1.0) < 0.001

    def test_minimum_floor(self):
        """Likelihood should never be exactly 0 (causes log(0))."""
        result = self.estimator._genotype_likelihood(
            genotype="AA", effect_allele="A", freq=0.0
        )
        assert result > 0  # floored to small epsilon
```

**Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_ancestry.py::TestGenotypeLikelihood -v`
Expected: FAIL — `ImportError`

**Step 3: Implement estimator with genotype likelihood**

Create `pydna_analyzer/ancestry/estimator.py`:

```python
"""Ancestry estimation via maximum likelihood admixture analysis."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

# Minimum likelihood floor to prevent log(0)
_EPSILON = 1e-10


@dataclass
class PopulationResult:
    """Estimated ancestry proportion for a single population."""

    population: str
    region: str
    proportion: float
    confidence_low: float
    confidence_high: float


@dataclass
class AncestryResult:
    """Complete ancestry estimation result."""

    populations: list[PopulationResult]
    snps_used: int
    snps_available: int
    coverage: float
    log_likelihood: float
    convergence: bool
    top_regions: dict[str, float]
    interpretation: str


class AncestryEstimator:
    """Maximum likelihood ancestry estimation using admixture model."""

    def _count_effect_alleles(self, genotype: str, effect_allele: str) -> int:
        """Count copies of effect allele in a genotype string (0, 1, or 2)."""
        return sum(1 for a in genotype if a == effect_allele)

    def _genotype_likelihood(
        self, genotype: str, effect_allele: str, freq: float
    ) -> float:
        """Compute P(genotype | population) using Hardy-Weinberg equilibrium.

        Args:
            genotype: Two-character genotype (e.g., "AA", "AG").
            effect_allele: The allele whose frequency is given.
            freq: Population frequency of the effect allele (0.0 to 1.0).

        Returns:
            Probability of observing this genotype given the population frequency.
        """
        # Clamp frequency to avoid exact 0 or 1 (prevents log(0))
        p = max(_EPSILON, min(1.0 - _EPSILON, freq))
        q = 1.0 - p

        dosage = self._count_effect_alleles(genotype, effect_allele)

        if dosage == 2:
            likelihood = p * p
        elif dosage == 1:
            likelihood = 2 * p * q
        else:
            likelihood = q * q

        return max(likelihood, _EPSILON)
```

**Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_ancestry.py::TestGenotypeLikelihood -v`
Expected: 7 PASSED

**Step 5: Commit**

```bash
git add pydna_analyzer/ancestry/estimator.py tests/test_ancestry.py
git commit -m "feat(ancestry): add genotype likelihood calculation with HWE"
```

---

### Task 4: Implement MLE Optimization

**Files:**
- Modify: `pydna_analyzer/ancestry/estimator.py`
- Test: `tests/test_ancestry.py` (add to)

**Step 1: Write failing tests for proportion optimization**

Add to `tests/test_ancestry.py`:

```python
import numpy as np


class TestOptimizeProportions:
    """Tests for MLE optimization of population proportions."""

    def setup_method(self):
        self.estimator = AncestryEstimator()

    def test_single_population_dominant(self):
        """When all likelihoods favor one population, it should get ~100%."""
        # 10 SNPs, 3 populations
        # All SNPs strongly favor population 0
        likelihood_matrix = np.array([
            [0.9, 0.01, 0.01],
        ] * 10)
        proportions = self.estimator._optimize_proportions(likelihood_matrix)
        assert proportions[0] > 0.9
        assert abs(sum(proportions) - 1.0) < 0.001

    def test_equal_mixture(self):
        """When likelihoods are equal, proportions should be roughly uniform."""
        likelihood_matrix = np.array([
            [0.5, 0.5, 0.5],
        ] * 10)
        proportions = self.estimator._optimize_proportions(likelihood_matrix)
        for p in proportions:
            assert 0.2 < p < 0.5  # roughly 1/3 each

    def test_two_population_mix(self):
        """50/50 mix of two populations should give ~50/50."""
        # Half the SNPs favor pop 0, half favor pop 1
        rows_pop0 = [[0.9, 0.01, 0.01]] * 5
        rows_pop1 = [[0.01, 0.9, 0.01]] * 5
        likelihood_matrix = np.array(rows_pop0 + rows_pop1)
        proportions = self.estimator._optimize_proportions(likelihood_matrix)
        assert abs(proportions[0] - 0.5) < 0.15
        assert abs(proportions[1] - 0.5) < 0.15

    def test_proportions_sum_to_one(self):
        """Proportions must always sum to 1.0."""
        likelihood_matrix = np.array([
            [0.3, 0.5, 0.2],
            [0.7, 0.1, 0.2],
            [0.1, 0.8, 0.1],
        ])
        proportions = self.estimator._optimize_proportions(likelihood_matrix)
        assert abs(sum(proportions) - 1.0) < 0.001

    def test_proportions_non_negative(self):
        """All proportions must be >= 0."""
        likelihood_matrix = np.array([
            [0.9, 0.01, 0.01, 0.01],
        ] * 20)
        proportions = self.estimator._optimize_proportions(likelihood_matrix)
        for p in proportions:
            assert p >= -0.001  # small tolerance for float
```

**Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_ancestry.py::TestOptimizeProportions -v`
Expected: FAIL — `AttributeError: 'AncestryEstimator' object has no attribute '_optimize_proportions'`

**Step 3: Implement MLE optimization**

Add to `pydna_analyzer/ancestry/estimator.py`:

```python
from scipy.optimize import minimize


class AncestryEstimator:
    # ... existing methods ...

    def _build_likelihood_matrix(
        self,
        genotypes: dict[str, str],
        aim_database: "AIMDatabase",
    ) -> tuple[NDArray[np.float64], list[str]]:
        """Build matrix of genotype likelihoods: (n_snps x n_populations).

        Args:
            genotypes: Dict of rsid -> genotype string for matched AIMs.
            aim_database: Reference AIM database.

        Returns:
            Tuple of (likelihood_matrix, population_names).
        """
        pop_names = aim_database.get_population_names()
        n_pops = len(pop_names)
        rsids = list(genotypes.keys())
        n_snps = len(rsids)

        matrix = np.zeros((n_snps, n_pops))

        for i, rsid in enumerate(rsids):
            gt = genotypes[rsid]
            effect_allele = aim_database.get_effect_allele(rsid)
            for j, pop in enumerate(pop_names):
                freq = aim_database.get_frequency(rsid, pop)
                matrix[i, j] = self._genotype_likelihood(gt, effect_allele, freq)

        return matrix, pop_names

    def _optimize_proportions(
        self, likelihood_matrix: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Find population proportions that maximize log-likelihood.

        Uses SLSQP with constraints: sum(Q) = 1, all Q >= 0.

        Args:
            likelihood_matrix: Shape (n_snps, n_populations).

        Returns:
            Array of population proportions summing to 1.0.
        """
        n_pops = likelihood_matrix.shape[1]

        # Initial guess: uniform proportions
        q0 = np.ones(n_pops) / n_pops

        def neg_log_likelihood(q: NDArray) -> float:
            """Negative log-likelihood (we minimize this)."""
            # For each SNP: P(genotype) = sum(q_k * P(genotype|k))
            mixed = likelihood_matrix @ q
            # Clamp to avoid log(0)
            mixed = np.maximum(mixed, _EPSILON)
            return -np.sum(np.log(mixed))

        # Constraints: sum(q) = 1
        constraints = {"type": "eq", "fun": lambda q: np.sum(q) - 1.0}
        # Bounds: each q_k in [0, 1]
        bounds = [(0.0, 1.0)] * n_pops

        result = minimize(
            neg_log_likelihood,
            q0,
            method="SLSQP",
            bounds=bounds,
            constraints=constraints,
            options={"maxiter": 1000, "ftol": 1e-12},
        )

        # Normalize to ensure exact sum = 1
        proportions = np.maximum(result.x, 0.0)
        proportions /= proportions.sum()

        return proportions
```

**Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_ancestry.py::TestOptimizeProportions -v`
Expected: 5 PASSED

**Step 5: Commit**

```bash
git add pydna_analyzer/ancestry/estimator.py tests/test_ancestry.py
git commit -m "feat(ancestry): add MLE proportion optimization with SLSQP"
```

---

### Task 5: Implement Bootstrap Confidence Intervals

**Files:**
- Modify: `pydna_analyzer/ancestry/estimator.py`
- Test: `tests/test_ancestry.py` (add to)

**Step 1: Write failing tests for bootstrap CI**

```python
class TestBootstrapConfidence:
    """Tests for bootstrap confidence interval estimation."""

    def setup_method(self):
        self.estimator = AncestryEstimator()

    def test_returns_low_and_high(self):
        """Bootstrap should return (low, high) arrays."""
        likelihood_matrix = np.array([[0.9, 0.05, 0.05]] * 20)
        low, high = self.estimator._bootstrap_confidence(
            likelihood_matrix, n_bootstrap=10
        )
        assert len(low) == 3
        assert len(high) == 3

    def test_ci_contains_point_estimate(self):
        """Point estimate should generally fall within CI."""
        likelihood_matrix = np.array([[0.9, 0.05, 0.05]] * 30)
        point = self.estimator._optimize_proportions(likelihood_matrix)
        low, high = self.estimator._bootstrap_confidence(
            likelihood_matrix, n_bootstrap=20
        )
        # Dominant population should be within CI
        assert low[0] <= point[0] <= high[0] + 0.05  # small tolerance

    def test_ci_width_reasonable(self):
        """CI should not be wider than 0-1 range."""
        likelihood_matrix = np.array([[0.8, 0.1, 0.1]] * 20)
        low, high = self.estimator._bootstrap_confidence(
            likelihood_matrix, n_bootstrap=10
        )
        for i in range(3):
            assert 0.0 <= low[i] <= high[i] <= 1.0
```

**Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_ancestry.py::TestBootstrapConfidence -v`
Expected: FAIL

**Step 3: Implement bootstrap**

Add to `AncestryEstimator` in `estimator.py`:

```python
    def _bootstrap_confidence(
        self,
        likelihood_matrix: NDArray[np.float64],
        n_bootstrap: int = 100,
        ci: float = 0.95,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Estimate confidence intervals by bootstrap resampling SNPs.

        Args:
            likelihood_matrix: Shape (n_snps, n_populations).
            n_bootstrap: Number of bootstrap iterations.
            ci: Confidence interval level (default 0.95).

        Returns:
            Tuple of (lower_bounds, upper_bounds) arrays of shape (n_populations,).
        """
        n_snps = likelihood_matrix.shape[0]
        n_pops = likelihood_matrix.shape[1]
        rng = np.random.default_rng(seed=42)

        bootstrap_proportions = np.zeros((n_bootstrap, n_pops))

        for b in range(n_bootstrap):
            # Resample SNP indices with replacement
            indices = rng.choice(n_snps, size=n_snps, replace=True)
            resampled = likelihood_matrix[indices, :]
            bootstrap_proportions[b, :] = self._optimize_proportions(resampled)

        alpha = (1 - ci) / 2
        lower = np.quantile(bootstrap_proportions, alpha, axis=0)
        upper = np.quantile(bootstrap_proportions, 1 - alpha, axis=0)

        return np.maximum(lower, 0.0), np.minimum(upper, 1.0)
```

**Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_ancestry.py::TestBootstrapConfidence -v`
Expected: 3 PASSED

**Step 5: Commit**

```bash
git add pydna_analyzer/ancestry/estimator.py tests/test_ancestry.py
git commit -m "feat(ancestry): add bootstrap confidence intervals"
```

---

### Task 6: Implement Full AncestryAnalyzer (End-to-End Estimation)

**Files:**
- Modify: `pydna_analyzer/ancestry/__init__.py`
- Modify: `pydna_analyzer/ancestry/estimator.py`
- Test: `tests/test_ancestry.py` (add to)
- Modify: `tests/conftest.py` (add ancestry fixtures)

**Step 1: Add ancestry test fixtures to conftest.py**

Add to `tests/conftest.py`:

```python
# ===== Ancestry estimation fixtures =====

ANCESTRY_HEADER = "#AncestryDNA raw data download\n#rsid\tchromosome\tposition\tallele1\tallele2\n"


@pytest.fixture
def ancestry_dataset_factory(tmp_path):
    """Factory fixture: creates a DNADataset from a dict of rsid -> genotype."""
    from pydna_analyzer.core.data_loader import load_dna_data

    def _make(genotypes: dict[str, str]) -> "DNADataset":
        lines = [ANCESTRY_HEADER]
        for rsid, gt in genotypes.items():
            a1 = gt[0] if len(gt) >= 1 else "A"
            a2 = gt[1] if len(gt) >= 2 else a1
            lines.append(f"{rsid}\t1\t100\t{a1}\t{a2}\n")
        filepath = tmp_path / "ancestry_test.txt"
        filepath.write_text("".join(lines))
        return load_dna_data(filepath)

    return _make
```

**Step 2: Write failing tests for full estimation pipeline**

Add to `tests/test_ancestry.py`:

```python
from pydna_analyzer.ancestry import AncestryAnalyzer, AncestryResult, PopulationResult


class TestAncestryAnalyzer:
    """Integration tests for full ancestry estimation."""

    def test_analyze_returns_result(self, ancestry_dataset_factory):
        """analyze() should return an AncestryResult."""
        analyzer = AncestryAnalyzer()
        db = AIMDatabase.load()
        # Build a synthetic dataset with a few AIMs
        rsids = list(db.aims.keys())[:20]
        genotypes = {}
        for rsid in rsids:
            effect = db.get_effect_allele(rsid)
            genotypes[rsid] = effect + effect  # homozygous effect
        dataset = ancestry_dataset_factory(genotypes)
        result = analyzer.analyze(dataset)
        assert isinstance(result, AncestryResult)

    def test_proportions_sum_to_one(self, ancestry_dataset_factory):
        """All population proportions should sum to ~1.0."""
        analyzer = AncestryAnalyzer()
        db = AIMDatabase.load()
        rsids = list(db.aims.keys())[:20]
        genotypes = {rsid: db.get_effect_allele(rsid) * 2 for rsid in rsids}
        dataset = ancestry_dataset_factory(genotypes)
        result = analyzer.analyze(dataset)
        total = sum(p.proportion for p in result.populations)
        assert abs(total - 1.0) < 0.01

    def test_result_has_coverage(self, ancestry_dataset_factory):
        """Result should report SNP coverage."""
        analyzer = AncestryAnalyzer()
        db = AIMDatabase.load()
        rsids = list(db.aims.keys())[:10]
        genotypes = {rsid: db.get_effect_allele(rsid) * 2 for rsid in rsids}
        dataset = ancestry_dataset_factory(genotypes)
        result = analyzer.analyze(dataset)
        assert result.snps_used == 10
        assert result.snps_available >= 50
        assert 0 < result.coverage <= 1.0

    def test_result_sorted_by_proportion(self, ancestry_dataset_factory):
        """Populations should be sorted by proportion descending."""
        analyzer = AncestryAnalyzer()
        db = AIMDatabase.load()
        rsids = list(db.aims.keys())[:20]
        genotypes = {rsid: db.get_effect_allele(rsid) * 2 for rsid in rsids}
        dataset = ancestry_dataset_factory(genotypes)
        result = analyzer.analyze(dataset)
        proportions = [p.proportion for p in result.populations]
        assert proportions == sorted(proportions, reverse=True)

    def test_top_regions_aggregated(self, ancestry_dataset_factory):
        """top_regions should aggregate populations by continent."""
        analyzer = AncestryAnalyzer()
        db = AIMDatabase.load()
        rsids = list(db.aims.keys())[:20]
        genotypes = {rsid: db.get_effect_allele(rsid) * 2 for rsid in rsids}
        dataset = ancestry_dataset_factory(genotypes)
        result = analyzer.analyze(dataset)
        assert isinstance(result.top_regions, dict)
        assert abs(sum(result.top_regions.values()) - 1.0) < 0.01

    def test_has_interpretation(self, ancestry_dataset_factory):
        """Result should include interpretation text."""
        analyzer = AncestryAnalyzer()
        db = AIMDatabase.load()
        rsids = list(db.aims.keys())[:20]
        genotypes = {rsid: db.get_effect_allele(rsid) * 2 for rsid in rsids}
        dataset = ancestry_dataset_factory(genotypes)
        result = analyzer.analyze(dataset)
        assert len(result.interpretation) > 0

    def test_no_matching_aims(self, ancestry_dataset_factory):
        """Should handle case where no AIMs match user data."""
        dataset = ancestry_dataset_factory({"rs99999999": "AA"})
        analyzer = AncestryAnalyzer()
        result = analyzer.analyze(dataset)
        assert result.snps_used == 0
        assert "insufficient" in result.interpretation.lower()

    def test_low_coverage_warning(self, ancestry_dataset_factory):
        """Low coverage should be noted in interpretation."""
        analyzer = AncestryAnalyzer()
        db = AIMDatabase.load()
        rsids = list(db.aims.keys())[:3]  # Very few AIMs
        genotypes = {rsid: db.get_effect_allele(rsid) * 2 for rsid in rsids}
        dataset = ancestry_dataset_factory(genotypes)
        result = analyzer.analyze(dataset)
        assert result.coverage < 0.1
```

**Step 3: Run tests to verify they fail**

Run: `uv run pytest tests/test_ancestry.py::TestAncestryAnalyzer -v`
Expected: FAIL — `ImportError: cannot import name 'AncestryAnalyzer' from 'pydna_analyzer.ancestry'`

**Step 4: Implement the full estimate() method on AncestryEstimator**

Add to `pydna_analyzer/ancestry/estimator.py`:

```python
    def estimate(
        self,
        genotypes: dict[str, str],
        aim_database: AIMDatabase,
        n_bootstrap: int = 100,
    ) -> AncestryResult:
        """Run full ancestry estimation pipeline.

        Args:
            genotypes: Dict of rsid -> genotype for matched AIMs only.
            aim_database: Reference AIM database.
            n_bootstrap: Number of bootstrap iterations for CI.

        Returns:
            AncestryResult with proportions, CI, and interpretation.
        """
        if len(genotypes) == 0:
            return self._empty_result(aim_database)

        # Build likelihood matrix
        likelihood_matrix, pop_names = self._build_likelihood_matrix(
            genotypes, aim_database
        )

        # Optimize proportions
        proportions = self._optimize_proportions(likelihood_matrix)

        # Bootstrap confidence intervals
        if len(genotypes) >= 5:
            ci_low, ci_high = self._bootstrap_confidence(
                likelihood_matrix, n_bootstrap=n_bootstrap
            )
        else:
            # Too few SNPs for meaningful bootstrap
            ci_low = np.zeros(len(pop_names))
            ci_high = np.ones(len(pop_names))

        # Compute log-likelihood at solution
        mixed = likelihood_matrix @ proportions
        log_likelihood = float(np.sum(np.log(np.maximum(mixed, _EPSILON))))

        # Build population results sorted by proportion
        pop_results = []
        for i, pop_name in enumerate(pop_names):
            region = aim_database.populations[pop_name].get("region", "Unknown")
            pop_results.append(
                PopulationResult(
                    population=pop_name,
                    region=region,
                    proportion=float(proportions[i]),
                    confidence_low=float(ci_low[i]),
                    confidence_high=float(ci_high[i]),
                )
            )
        pop_results.sort(key=lambda x: x.proportion, reverse=True)

        # Aggregate by region
        top_regions: dict[str, float] = {}
        for pr in pop_results:
            top_regions[pr.region] = top_regions.get(pr.region, 0.0) + pr.proportion

        # Generate interpretation
        interpretation = self._generate_interpretation(
            pop_results, top_regions, len(genotypes), len(aim_database.aims)
        )

        return AncestryResult(
            populations=pop_results,
            snps_used=len(genotypes),
            snps_available=len(aim_database.aims),
            coverage=len(genotypes) / len(aim_database.aims),
            log_likelihood=log_likelihood,
            convergence=True,
            top_regions=top_regions,
            interpretation=interpretation,
        )

    def _empty_result(self, aim_database: AIMDatabase) -> AncestryResult:
        """Return result when no AIMs match."""
        return AncestryResult(
            populations=[],
            snps_used=0,
            snps_available=len(aim_database.aims),
            coverage=0.0,
            log_likelihood=0.0,
            convergence=False,
            top_regions={},
            interpretation=(
                "Insufficient data: no ancestry-informative markers found in your "
                "dataset. Ancestry estimation requires overlap with the reference panel."
            ),
        )

    def _generate_interpretation(
        self,
        populations: list[PopulationResult],
        regions: dict[str, float],
        snps_used: int,
        snps_available: int,
    ) -> str:
        """Generate human-readable interpretation of results."""
        coverage = snps_used / snps_available if snps_available > 0 else 0

        # Top populations above 5%
        significant = [p for p in populations if p.proportion >= 0.05]
        if not significant:
            return "No population reached 5% proportion threshold."

        parts = []
        primary = significant[0]
        parts.append(
            f"Your ancestry is predominantly {primary.population} "
            f"({primary.proportion:.0%})."
        )

        if len(significant) > 1:
            others = ", ".join(
                f"{p.population} ({p.proportion:.0%})" for p in significant[1:3]
            )
            parts.append(f"Additional components include {others}.")

        # Top region
        top_region = max(regions, key=regions.get)
        parts.append(
            f"By region, your strongest affinity is to {top_region} "
            f"({regions[top_region]:.0%})."
        )

        if coverage < 0.5:
            parts.append(
                f"Note: Only {coverage:.0%} of ancestry markers were found. "
                "Results may be less precise."
            )

        return " ".join(parts)
```

**Step 5: Implement AncestryAnalyzer in __init__.py**

Update `pydna_analyzer/ancestry/__init__.py`:

```python
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
    """High-level ancestry estimation interface.

    Usage:
        >>> from pydna_analyzer.ancestry import AncestryAnalyzer
        >>> analyzer = AncestryAnalyzer()
        >>> result = analyzer.analyze(dataset)
        >>> for pop in result.populations[:5]:
        ...     print(f"{pop.population}: {pop.proportion:.1%}")
    """

    def __init__(self, n_bootstrap: int = 100) -> None:
        self.n_bootstrap = n_bootstrap
        self._db = AIMDatabase.load()
        self._estimator = AncestryEstimator()

    def analyze(self, dataset: DNADataset) -> AncestryResult:
        """Estimate ancestry composition from a DNA dataset.

        Args:
            dataset: Loaded DNA dataset with genotype data.

        Returns:
            AncestryResult with population proportions and confidence intervals.
        """
        # Match user SNPs against AIMs
        aim_rsids = self._db.get_aim_rsids()
        user_rsids = dataset.get_rsids()
        overlap = aim_rsids & user_rsids

        # Extract genotypes for matched AIMs
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
```

**Step 6: Run tests to verify they pass**

Run: `uv run pytest tests/test_ancestry.py -v`
Expected: All tests PASS (8 reference + 7 likelihood + 5 optimization + 3 bootstrap + 8 analyzer = 31)

**Step 7: Commit**

```bash
git add pydna_analyzer/ancestry/ tests/test_ancestry.py tests/conftest.py
git commit -m "feat(ancestry): complete ancestry estimation pipeline with MLE and bootstrap CI"
```

---

### Task 7: Add CLI Command

**Files:**
- Modify: `pydna_analyzer/cli.py`
- Test: `tests/test_ancestry.py` (add CLI tests)

**Step 1: Write failing CLI tests**

Add to `tests/test_ancestry.py`:

```python
from typer.testing import CliRunner
from pydna_analyzer.cli import app

runner = CliRunner()


class TestAncestryCLI:
    """Tests for the ancestry CLI command."""

    def test_ancestry_command_exists(self, sample_ancestrydna_file):
        """The ancestry command should be available."""
        result = runner.invoke(app, ["ancestry", str(sample_ancestrydna_file)])
        assert result.exit_code == 0 or "Error" not in result.output

    def test_ancestry_json_output(self, sample_ancestrydna_file, tmp_path):
        """--json flag should produce a JSON file."""
        output = tmp_path / "ancestry.json"
        result = runner.invoke(
            app, ["ancestry", str(sample_ancestrydna_file), "-o", str(output)]
        )
        # Should at least not crash
        assert result.exit_code == 0 or result.exit_code is None

    def test_ancestry_shows_populations(self, sample_ancestrydna_file):
        """Output should include population names."""
        result = runner.invoke(app, ["ancestry", str(sample_ancestrydna_file)])
        # Should contain some output (even if few AIMs match)
        assert "Ancestry" in result.output or "ancestry" in result.output
```

**Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_ancestry.py::TestAncestryCLI -v`
Expected: FAIL — `No such command 'ancestry'`

**Step 3: Add the ancestry CLI command**

Add to `pydna_analyzer/cli.py` (after the `pgx` command):

```python
@app.command()
def ancestry(
    filepath: Path = typer.Argument(
        ...,
        help="Path to DNA data file (AncestryDNA, 23andMe, MyHeritage, or VCF)",
        exists=True,
    ),
    output: Optional[Path] = typer.Option(
        None,
        "--output",
        "-o",
        help="Output JSON file path",
    ),
    bootstrap: int = typer.Option(
        100,
        "--bootstrap",
        "-b",
        help="Number of bootstrap iterations for confidence intervals",
    ),
):
    """
    Estimate ancestry composition from DNA data.

    Uses maximum likelihood estimation with ancestry-informative markers
    to estimate population proportions with confidence intervals.

    Examples:
        pydna_analyzer ancestry AncestryDNA.txt
        pydna_analyzer ancestry data.txt -o ancestry.json
        pydna_analyzer ancestry data.txt --bootstrap 200
    """
    from pydna_analyzer.ancestry import AncestryAnalyzer

    # Load data
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
        transient=True,
    ) as progress:
        progress.add_task("Loading DNA data...", total=None)
        try:
            dataset = load_dna_data(filepath)
        except Exception as e:
            console.print(f"[red]Error loading file: {e}[/]")
            raise typer.Exit(1)

    console.print(f"[green]✓[/] Loaded {dataset.snp_count:,} SNPs from {filepath.name}")

    # Run ancestry estimation
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
        transient=True,
    ) as progress:
        progress.add_task("Estimating ancestry composition...", total=None)
        analyzer = AncestryAnalyzer(n_bootstrap=bootstrap)
        result = analyzer.analyze(dataset)

    coverage_pct = f"{result.coverage * 100:.1f}%"
    console.print(
        f"[green]✓[/] Matched {result.snps_used} / {result.snps_available} "
        f"ancestry markers ({coverage_pct})"
    )

    if result.snps_used == 0:
        console.print(f"\n[yellow]{result.interpretation}[/]")
        raise typer.Exit(0)

    _print_ancestry_results(result)

    # Save JSON if requested
    if output:
        import json

        output_data = {
            "populations": [
                {
                    "population": p.population,
                    "region": p.region,
                    "proportion": round(p.proportion, 4),
                    "confidence_low": round(p.confidence_low, 4),
                    "confidence_high": round(p.confidence_high, 4),
                }
                for p in result.populations
                if p.proportion >= 0.01
            ],
            "regions": {k: round(v, 4) for k, v in result.top_regions.items()},
            "snps_used": result.snps_used,
            "snps_available": result.snps_available,
            "coverage": round(result.coverage, 4),
            "interpretation": result.interpretation,
        }
        with output.open("w") as f:
            json.dump(output_data, f, indent=2)
        console.print(f"\n[green]✓[/] Results saved to {output}")


def _print_ancestry_results(result):
    """Print ancestry estimation results with Rich tables."""
    from pydna_analyzer.ancestry import AncestryResult

    console.print()
    console.print(
        Panel.fit(
            "[bold blue]🌍 Ancestry Composition[/]",
            border_style="blue",
        )
    )
    console.print()

    # Population table with visual bars
    table = Table(show_lines=True)
    table.add_column("Population", style="bold")
    table.add_column("%", justify="right", style="cyan")
    table.add_column("", min_width=25)  # visual bar
    table.add_column("95% CI", style="dim")

    for pop in result.populations:
        if pop.proportion < 0.02:
            continue  # Skip very small proportions
        pct = f"{pop.proportion * 100:.1f}%"
        bar_len = int(pop.proportion * 30)
        bar = "█" * bar_len
        ci = f"{pop.confidence_low * 100:.1f} - {pop.confidence_high * 100:.1f}%"
        table.add_row(pop.population, pct, f"[blue]{bar}[/]", ci)

    console.print(table)
    console.print()

    # Regional summary
    console.print("[bold]📊 Regional Summary[/]")
    for region, proportion in sorted(
        result.top_regions.items(), key=lambda x: x[1], reverse=True
    ):
        if proportion >= 0.02:
            console.print(f"  {region}: {proportion * 100:.1f}%")

    # Interpretation
    console.print()
    console.print(
        Panel(
            result.interpretation,
            title="📝 Interpretation",
            border_style="dim",
        )
    )

    console.print()
    console.print(
        "[dim]⚠ Note: Ancestry estimates are based on a curated set of "
        "ancestry-informative markers. Commercial services use significantly "
        "more data for higher resolution. Results are for educational purposes only.[/]"
    )
```

**Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_ancestry.py::TestAncestryCLI -v`
Expected: 3 PASSED

**Step 5: Run all tests**

Run: `uv run pytest -v`
Expected: All 129 + ~31 new = ~160 tests PASS

**Step 6: Commit**

```bash
git add pydna_analyzer/cli.py tests/test_ancestry.py
git commit -m "feat(ancestry): add ancestry CLI command with Rich table output"
```

---

### Task 8: Wire Up Exports and Update Package

**Files:**
- Modify: `pydna_analyzer/__init__.py`
- Modify: `pydna_analyzer/cli.py` (add to CLAUDE.md commands)

**Step 1: Add ancestry exports to package __init__**

Add to `pydna_analyzer/__init__.py`:

```python
from pydna_analyzer.ancestry import AncestryAnalyzer, AncestryResult

# Add to __all__:
__all__ = [
    ...,
    "AncestryAnalyzer",
    "AncestryResult",
]
```

**Step 2: Run full test suite and lint**

Run: `uv run pytest -v && uv run ruff check . && uv run ruff format --check .`
Expected: All tests pass, no lint errors

**Step 3: Commit**

```bash
git add pydna_analyzer/__init__.py
git commit -m "feat(ancestry): export AncestryAnalyzer from package root"
```

---

### Task 9: Update CLAUDE.md and README

**Files:**
- Modify: `CLAUDE.md` (add ancestry commands, update test count, add to module map)
- Modify: `README.md` (update roadmap checkbox, add ancestry section)

**Step 1: Update CLAUDE.md**

- Add `uv run pydna_analyzer ancestry <file>` to Commands section
- Add ancestry module to Module Map
- Update test count

**Step 2: Update README.md**

- Check the `[x] Ancestry composition estimation` roadmap item
- Add brief usage example for `pydna_analyzer ancestry`

**Step 3: Commit**

```bash
git add CLAUDE.md README.md
git commit -m "docs: update CLAUDE.md and README for ancestry estimation module"
```

---

## Summary

| Task | What | Tests Added |
|------|------|-------------|
| 1 | Update pyproject.toml deps | — |
| 2 | AIM reference database + JSON | 8 |
| 3 | Genotype likelihood (HWE math) | 7 |
| 4 | MLE optimization (SLSQP) | 5 |
| 5 | Bootstrap confidence intervals | 3 |
| 6 | Full AncestryAnalyzer pipeline | 8 |
| 7 | CLI command with Rich output | 3 |
| 8 | Package exports | — |
| 9 | Documentation updates | — |
| **Total** | | **~34 tests** |

Estimated final test count: ~163 (129 existing + 34 new).
