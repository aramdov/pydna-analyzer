# Ancestry Estimation Module Design

**Date:** 2026-02-14
**Phase:** 5
**Status:** Approved

## Overview

Add sub-continental ancestry composition estimation to PyDNA Analyzer using likelihood-based admixture analysis with bundled reference allele frequencies from 1000 Genomes Phase 3.

## Goals

- Estimate ancestry proportions across ~15-20 sub-continental populations
- Provide confidence intervals via bootstrap resampling
- Bundle all reference data (~100-200KB JSON) for fully offline analysis
- Follow existing module patterns (pydantic models, Rich CLI, pytest)

## Non-Goals

- Chromosome painting (requires phased data)
- Fine-grained within-country resolution (needs 1000s of reference samples)
- Ancient DNA ancestry
- Competing with commercial services on accuracy

## Algorithm

**Maximum Likelihood Estimation (MLE) admixture:**

1. Match user's SNPs against ~300-500 curated ancestry-informative markers (AIMs)
2. For each SNP, compute genotype likelihood per population using Hardy-Weinberg:
   - P(AA|k) = p_k^2, P(AB|k) = 2*p_k*(1-p_k), P(BB|k) = (1-p_k)^2
3. Find population mixture proportions Q that maximize:
   - sum(log(sum(q_k * P(genotype_i | k)))) across all SNPs
4. Constraints: q_k >= 0, sum(q_k) = 1
5. Solver: scipy.optimize.minimize with SLSQP
6. Bootstrap resample SNPs (n=100) for 95% confidence intervals

## Architecture

```
pydna_analyzer/ancestry/
├── __init__.py          # Public API: AncestryAnalyzer, AncestryResult, PopulationResult
├── reference_data.py    # AIMDatabase: loads/validates bundled JSON
├── estimator.py         # AncestryEstimator: MLE optimization + bootstrap CI
└── data/
    └── aim_frequencies.json  # Bundled reference allele frequencies
```

### Data Models

```python
@dataclass
class PopulationResult:
    population: str          # e.g. "Northern European"
    region: str              # e.g. "Europe"
    proportion: float        # 0.0 - 1.0
    confidence_low: float    # 95% CI lower
    confidence_high: float   # 95% CI upper

@dataclass
class AncestryResult:
    populations: list[PopulationResult]  # Sorted by proportion desc
    snps_used: int
    snps_available: int
    coverage: float
    log_likelihood: float
    convergence: bool
    top_regions: dict[str, float]  # Aggregated by continent
    interpretation: str
```

### Reference Data Format

```json
{
  "metadata": {
    "source": "1000 Genomes Phase 3",
    "populations": ["Northern European", ...],
    "aim_count": 350
  },
  "populations": {
    "Northern European": {"code": "NEU", "region": "Europe", "sample_size": 503}
  },
  "aims": {
    "rs1426654": {
      "chromosome": "15",
      "position": 48426484,
      "gene": "SLC24A5",
      "effect_allele": "A",
      "frequencies": {
        "Northern European": 0.98,
        "West African": 0.03
      }
    }
  }
}
```

### CLI Command

`pydna_analyzer ancestry <file>` with `--json` and `-o` options.

Output: Rich table with population percentages, visual bars, 95% CI, regional summary, and interpretation text.

## Dependencies

- numpy (array operations)
- scipy (optimization: scipy.optimize.minimize with SLSQP)
- No new external dependencies beyond what's likely already available

## Test Plan (~20-25 tests)

- Unit: genotype likelihood math (HWE)
- Unit: optimization with synthetic known-mixture data
- Integration: 100% single-population synthetic individual
- Integration: 50/50 synthetic mix
- Edge cases: low coverage, no matching AIMs, convergence failure
- Reference data integrity (all populations present, valid frequencies)

## Populations (~15-20)

| Region | Populations |
|--------|------------|
| Europe | Northern European, Southern European, Eastern European |
| Africa | West African, East African, North African |
| East Asia | Han Chinese, Japanese, Korean |
| South Asia | South Asian |
| Americas | Native American |
| Middle East | Middle Eastern |
| Oceania | Oceanian |
| Central Asia | Central Asian |
