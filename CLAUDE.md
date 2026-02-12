# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GenomeInsight is a privacy-first personal genomics toolkit that analyzes DNA data from consumer testing services (AncestryDNA, 23andMe, MyHeritage, VCF). All processing is local — no data is uploaded.

## Commands

```bash
# Install dependencies
uv sync                    # Base dependencies
uv sync --all-extras       # All optional features (ai, prs, pgx, ancestry, dev)

# Run CLI
uv run genomeinsight analyze <file>
uv run genomeinsight prs <file> --weights <weights.csv>
uv run genomeinsight pgx <file>                          # All 6 genes
uv run genomeinsight pgx <file> --gene CYP2C19            # Single gene
uv run genomeinsight pgx <file> -o results.json           # JSON export
uv run genomeinsight info <file>
uv run genomeinsight variants

# Tests
uv run pytest                                          # All tests
uv run pytest tests/test_prs.py                        # Single file
uv run pytest tests/test_prs.py::TestPRSCalculator     # Single class
uv run pytest tests/test_prs.py::TestPRSCalculator::test_calculate_returns_result  # Single test
uv run pytest --cov=genomeinsight                      # With coverage

# Linting & formatting
uv run ruff check .          # Lint
uv run ruff format .         # Format
uv run mypy genomeinsight    # Type check
```

## Architecture

### Data Flow

```
DNA File → DataLoader (auto-detect format) → DNADataset (pydantic) → Analyzer → Results → Reports
                                                  ├─→ ClinicalAnalyzer → AnalysisResult
                                                  ├─→ PGxAnalyzer      → PGxResult (diplotypes, phenotypes, drug recs)
                                                  └─→ PRSCalculator    → PRSResult
```

### Module Map

- **`cli.py`** — Typer CLI. Orchestrates all modules. Entry point: `genomeinsight.cli:app`
- **`core/data_loader.py`** — Multi-format parsers using strategy pattern (BaseLoader → AncestryDNALoader, TwentyThreeMeLoader, MyHeritageLoader, VCFLoader). `load_dna_data()` auto-detects format.
- **`clinical/`** — Variant analysis:
  - `variants.py` — Curated database of 35+ clinical SNPs as pydantic models
  - `analyzer.py` — ClinicalAnalyzer matches user genotypes against the database
  - `apoe.py` — Special-case APOE genotyping (rs429358 + rs7412 determine ε2/ε3/ε4 isoforms)
- **`polygenic/`** — PRS calculator. Loads external weight files (CSV with rsid, effect_allele, weight columns).
- **`pharmacogenomics/`** — PGx star allele calling and drug recommendations:
  - `__init__.py` — Full module: enums (`MetabolizerPhenotype`, `CPICLevel`), gene database (6 genes), drug recommendation database (~20 entries), `PGxAnalyzer` class
  - Pipeline: extract per-gene SNPs → call star alleles → form diplotype → compute activity score → predict phenotype → lookup drug recs
  - Two analysis paths: multi-SNP genes (star allele calling with CPIC activity scores) and single-SNP genes (VKORC1, COMT direct genotype map)
  - Self-contained — no external PGx dependencies; works with consumer DTC data
- **`ai/`** — LLM report generation:
  - `client.py` — Abstract LLM client with OpenAI and Anthropic implementations
  - `report_generator.py` — Generates natural language reports from analysis results
- **`reports/`** — HTML (Jinja2 + Plotly) and JSON export

### Key Patterns

- **Pydantic models** for all data structures (DNADataset, ClinicalVariant, etc.)
- **Dataclasses** for results (VariantResult, AnalysisResult, PRSResult)
- **Enums** for risk levels, evidence levels, and variant categories
- **Rich** for terminal output (tables, panels, progress spinners)
- Type hints required on all function signatures (enforced by mypy strict mode)

### Test Structure

Tests use class-based organization with pytest fixtures in `conftest.py` for sample DNA data. The `tmp_path` fixture handles temporary file creation. Markers: `slow`, `integration`. Current count: 129 tests (33 PGx, 22 AI, 8 APOE, 13 CLI, 7 clinical analyzer, 13 data loader, 9 PRS, 21 reports, 3 PGx database integrity).

## Code Style

- Ruff: line-length 100, Python 3.10+ target, double quotes, space indentation
- Ruff lint rules: E, F, W, I (isort), UP (pyupgrade), B (bugbear), SIM, C4, PTH (pathlib)
- When adding clinical SNPs: provide rsID, gene, risk allele, category, and cite peer-reviewed sources (PubMed IDs)
