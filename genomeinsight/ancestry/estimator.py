"""Maximum likelihood ancestry estimator using AIM frequencies."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import minimize

from genomeinsight.ancestry.reference_data import AIMDatabase


@dataclass
class PopulationResult:
    """Ancestry proportion for a single population."""

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
    top_regions: dict[str, float] = field(default_factory=dict)
    interpretation: str = ""


class AncestryEstimator:
    """Estimates ancestry proportions via maximum likelihood on AIM genotypes."""

    EPSILON: float = 1e-10

    def _count_effect_alleles(self, genotype: str, effect_allele: str) -> int:
        """Count 0, 1, or 2 copies of the effect allele in a genotype."""
        return sum(1 for allele in genotype if allele == effect_allele)

    def _genotype_likelihood(
        self, genotype: str, effect_allele: str, freq: float
    ) -> float:
        """Compute P(genotype | population) under Hardy-Weinberg equilibrium.

        For a population with effect-allele frequency p:
          - homozygous effect:  p^2
          - heterozygous:       2pq
          - homozygous other:   q^2
        Frequency is clamped to [epsilon, 1-epsilon] to avoid log(0).
        """
        eps = self.EPSILON
        p = max(eps, min(freq, 1.0 - eps))
        q = 1.0 - p

        n_effect = self._count_effect_alleles(genotype, effect_allele)
        if n_effect == 2:
            return p * p
        elif n_effect == 1:
            return 2.0 * p * q
        else:
            return q * q

    def _build_likelihood_matrix(
        self,
        genotypes: dict[str, str],
        aim_database: AIMDatabase,
    ) -> tuple[NDArray[np.float64], list[str]]:
        """Build an (n_snps x n_pops) matrix of genotype likelihoods.

        Returns the matrix and the list of population names (column order).
        """
        pop_names = aim_database.get_population_names()
        aim_rsids = aim_database.get_aim_rsids()
        shared_rsids = sorted(set(genotypes.keys()) & aim_rsids)

        n_snps = len(shared_rsids)
        n_pops = len(pop_names)
        matrix = np.zeros((n_snps, n_pops), dtype=np.float64)

        for i, rsid in enumerate(shared_rsids):
            gt = genotypes[rsid]
            effect_allele = aim_database.get_effect_allele(rsid)
            if effect_allele is None:
                continue
            for j, pop in enumerate(pop_names):
                freq = aim_database.get_frequency(rsid, pop)
                if freq is None:
                    freq = 0.5
                matrix[i, j] = self._genotype_likelihood(gt, effect_allele, freq)

        return matrix, pop_names

    def _optimize_proportions(
        self, likelihood_matrix: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Find ancestry proportions Q that maximise the log-likelihood.

        Uses SLSQP with constraints: sum(Q) = 1, Q >= 0.
        """
        n_pops = likelihood_matrix.shape[1]
        q0 = np.full(n_pops, 1.0 / n_pops)

        def neg_log_likelihood(q: NDArray[np.float64]) -> float:
            mixture = likelihood_matrix @ q
            mixture = np.maximum(mixture, self.EPSILON)
            return -np.sum(np.log(mixture))

        constraints = {"type": "eq", "fun": lambda q: np.sum(q) - 1.0}
        bounds = [(0.0, 1.0)] * n_pops

        result = minimize(
            neg_log_likelihood,
            q0,
            method="SLSQP",
            bounds=bounds,
            constraints=constraints,
            options={"maxiter": 1000, "ftol": 1e-12},
        )

        proportions = np.maximum(result.x, 0.0)
        proportions /= proportions.sum()
        return proportions

    def _bootstrap_confidence(
        self,
        likelihood_matrix: NDArray[np.float64],
        n_bootstrap: int = 100,
        ci: float = 0.95,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Compute bootstrap confidence intervals by resampling SNP rows.

        Returns (lower_bounds, upper_bounds) arrays of shape (n_pops,).
        """
        rng = np.random.default_rng(seed=42)
        n_snps = likelihood_matrix.shape[0]
        n_pops = likelihood_matrix.shape[1]
        samples = np.zeros((n_bootstrap, n_pops), dtype=np.float64)

        for b in range(n_bootstrap):
            indices = rng.choice(n_snps, size=n_snps, replace=True)
            resampled = likelihood_matrix[indices, :]
            samples[b, :] = self._optimize_proportions(resampled)

        alpha = (1.0 - ci) / 2.0
        lower = np.percentile(samples, 100.0 * alpha, axis=0)
        upper = np.percentile(samples, 100.0 * (1.0 - alpha), axis=0)
        return lower, upper
