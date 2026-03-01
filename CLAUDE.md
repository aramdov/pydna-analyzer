# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PyDNA Analyzer is a privacy-first personal genomics toolkit that analyzes DNA data from consumer testing services (AncestryDNA, 23andMe, MyHeritage, VCF). All processing is local — no data is uploaded.

## Commands

```bash
# Install dependencies
uv sync                    # Base dependencies
uv sync --all-extras       # All optional features (ai, prs, pgx, ancestry, api, dev)

# Run CLI
uv run pydna-analyzer analyze <file>
uv run pydna-analyzer prs <file> --weights <weights.csv>
uv run pydna-analyzer pgx <file>                          # All 6 genes
uv run pydna-analyzer pgx <file> --gene CYP2C19            # Single gene
uv run pydna-analyzer pgx <file> -o results.json           # JSON export
uv run pydna-analyzer ancestry <file>                     # Ancestry estimation
uv run pydna-analyzer ancestry <file> -o ancestry.json    # JSON export
uv run pydna-analyzer ancestry <file> --bootstrap 200     # More CI iterations
uv run pydna-analyzer info <file>
uv run pydna-analyzer variants
uv run pydna-analyzer serve                             # Start API on localhost:8000
uv run pydna-analyzer serve --port 9000                 # Custom port
uv run pydna-analyzer serve --reload                    # Dev mode with auto-reload

# Tests
uv run pytest                                          # All tests
uv run pytest tests/test_prs.py                        # Single file
uv run pytest tests/test_prs.py::TestPRSCalculator     # Single class
uv run pytest tests/test_prs.py::TestPRSCalculator::test_calculate_returns_result  # Single test
uv run pytest --cov=pydna_analyzer                      # With coverage

# Linting & formatting
uv run ruff check .          # Lint
uv run ruff format .         # Format
uv run mypy pydna_analyzer    # Type check
```

## Architecture

### Data Flow

```
DNA File → DataLoader (auto-detect format) → DNADataset (pydantic) → Analyzer → Results → Reports
                                                  ├─→ ClinicalAnalyzer → AnalysisResult
                                                  ├─→ PGxAnalyzer      → PGxResult (diplotypes, phenotypes, drug recs)
                                                  ├─→ PRSCalculator    → PRSResult
                                                  └─→ AncestryAnalyzer → AncestryResult (population proportions, CI)
```

### Module Map

- **`cli.py`** — Typer CLI. Orchestrates all modules. Entry point: `pydna_analyzer.cli:app`
- **`core/data_loader.py`** — Multi-format parsers using strategy pattern (BaseLoader → AncestryDNALoader, TwentyThreeMeLoader, MyHeritageLoader, VCFLoader). `load_dna_data()` auto-detects format. VCFLoader supports multi-sample selection, quality filtering (min_qual), and multi-allelic site handling.
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
- **`ancestry/`** — Ancestry composition estimation:
  - `reference_data.py` — AIMDatabase: loads bundled AIM allele frequencies (123 markers, 14 populations)
  - `estimator.py` — AncestryEstimator: MLE admixture via scipy SLSQP, bootstrap CI, HWE likelihood
  - `data/aim_frequencies.json` — Curated allele frequencies from 1000 Genomes Phase 3
  - Pipeline: match AIMs → build likelihood matrix → optimize proportions → bootstrap CI → interpret
- **`api/`** — REST API (FastAPI):
  - `__init__.py` — App factory (`create_app()`) with CORS configuration
  - `routes.py` — Endpoint definitions: /health, /info, /analyze, /pgx, /ancestry, /prs, /variants
  - `schemas.py` — File input resolution (path or upload) via `ResolvedFile` dataclass
- **`reports/`** — HTML (Jinja2 + Plotly) and JSON export

### Key Patterns

- **Pydantic models** for all data structures (DNADataset, ClinicalVariant, etc.)
- **Dataclasses** for results (VariantResult, AnalysisResult, PRSResult)
- **Enums** for risk levels, evidence levels, and variant categories
- **Rich** for terminal output (tables, panels, progress spinners)
- Type hints required on all function signatures (enforced by mypy strict mode)

### Test Structure

Tests use class-based organization with pytest fixtures in `conftest.py` for sample DNA data. The `tmp_path` fixture handles temporary file creation. Markers: `slow`, `integration`. Current count: 222 tests (25 API, 34 ancestry, 33 PGx, 22 AI, 8 APOE, 13 CLI, 7 clinical analyzer, 15 data loader, 34 VCF, 9 PRS, 21 reports, 3 PGx database integrity).

## Code Style

### Code Intelligence

Prefer LSP over Grep/Read for code navigation — it's faster, precise, and avoids reading entire files:
- `workspaceSymbol` to find where something is defined
- `findReferences` to see all usages across the codebase
- `goToDefinition` / `goToImplementation` to jump to source
- `hover` for type info without reading the file

Use Grep only when LSP isn't available or for text/pattern searches (comments, strings, config).

After writing or editing code, check LSP diagnostics and fix errors before proceeding.

- Ruff: line-length 100, Python 3.10+ target, double quotes, space indentation
- Ruff lint rules: E, F, W, I (isort), UP (pyupgrade), B (bugbear), SIM, C4, PTH (pathlib)
- When adding clinical SNPs: provide rsID, gene, risk allele, category, and cite peer-reviewed sources (PubMed IDs)
