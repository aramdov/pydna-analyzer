# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Initial project structure with modular architecture
- Multi-format DNA file parsing (AncestryDNA, 23andMe, MyHeritage)
- VCF file format support
- Clinical variant analysis with 35+ curated SNPs
- APOE genotyping with ε2/ε3/ε4 determination
- Gene interaction detection (MTHFR compound heterozygosity)
- Polygenic Risk Score calculator
- AI-powered natural language reports (OpenAI and Anthropic providers, consumer and technical report styles)
- HTML report generation with interactive charts
- JSON export functionality
- Command-line interface with Typer
- **Pharmacogenomics module (Phase 2)**:
  - Star allele calling for 6 pharmacogenes (CYP2C9, CYP2C19, CYP2D6, TPMT, VKORC1, COMT)
  - Diplotype calling from unphased consumer DNA data (*1 by absence)
  - CPIC-based activity score computation (0, 0.25, 0.5, 1.0, 1.5)
  - Metabolizer phenotype prediction (PM/IM/NM/RM/UM + VKORC1 sensitivity + COMT activity)
  - Drug recommendation database (~20 entries) covering clopidogrel, warfarin, codeine, tamoxifen, azathioprine, omeprazole, celecoxib, opioids
  - `genomeinsight pgx` CLI command with Rich-formatted gene summary, drug recommendations, and JSON export
  - Self-contained module — removed unused `pypgx` dependency
- Comprehensive test suite expanded to 129 tests (33 new PGx tests)

### Coming Soon
- Ancestry composition estimation
- REST API
- Web dashboard
