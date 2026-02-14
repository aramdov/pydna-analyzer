<h1 align="center">
  🧬 GenomeInsight
</h1>

<p align="center">
  <strong>A privacy-first personal genomics toolkit for analyzing DNA data from consumer testing services</strong>
</p>

<p align="center">
  <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/python-3.10+-blue.svg" alt="Python 3.10+"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-green.svg" alt="License: MIT"></a>
  <a href="https://github.com/genomeinsight/genomeinsight/issues"><img src="https://img.shields.io/github/issues/genomeinsight/genomeinsight" alt="Issues"></a>
</p>

---

GenomeInsight analyzes your raw DNA data from services like **AncestryDNA**, **23andMe**, and **MyHeritage** to provide insights into clinical variants, pharmacogenomics, and health-related genetic markers.

> **🔒 Privacy First**: All analysis runs 100% locally on your machine. Your genetic data never leaves your computer.

## ✨ Features

| Feature | Description |
|---------|-------------|
| 📁 **Multi-Format Support** | AncestryDNA, 23andMe, MyHeritage, and VCF files |
| 🏥 **Clinical Variant Analysis** | 35+ curated SNPs with evidence-based interpretations |
| 🧠 **APOE Genotyping** | Complete ε2/ε3/ε4 determination with risk assessment |
| 🔗 **Gene Interactions** | Detection of compound effects (e.g., MTHFR compound heterozygosity) |
| 💊 **Pharmacogenomics** | Drug metabolism variants (CYP2D6, CYP2C19, VKORC1, etc.) |
| 📈 **Polygenic Risk Scores** | Calculate PRS from GWAS weights (PGS Catalog compatible) |
| 🤖 **AI-Powered Reports** | Natural language reports via OpenAI or Anthropic |
| 📊 **Interactive Reports** | Beautiful HTML reports with charts and tables |
| 🌍 **Ancestry Estimation** | Sub-continental ancestry composition with confidence intervals |
| ✅ **Comprehensive Tests** | 163 tests covering all modules |

## 🚀 Quick Start

### Installation

**Using [uv](https://docs.astral.sh/uv/) (Recommended)**

```bash
# Clone the repository
git clone https://github.com/genomeinsight/genomeinsight.git
cd genomeinsight

# Create environment and install
uv sync

# Or with all optional dependencies
uv sync --all-extras
```

**Using pip**

```bash
# Clone the repository
git clone https://github.com/genomeinsight/genomeinsight.git
cd genomeinsight

# Create virtual environment
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install in development mode
pip install -e ".[dev]"
```

### Basic Usage

```bash
# Analyze your DNA data
genomeinsight analyze your_dna_file.txt

# Generate an HTML report
genomeinsight analyze your_dna_file.txt --html -o report.html

# Export to JSON
genomeinsight analyze your_dna_file.txt --json -o results.json

# Estimate ancestry composition
genomeinsight ancestry your_dna_file.txt

# Calculate Polygenic Risk Scores
genomeinsight prs your_dna_file.txt --weights cardiovascular_prs.csv

# Generate AI-powered report (requires API key)
genomeinsight analyze your_dna_file.txt --ai --ai-style consumer

# View file information
genomeinsight info your_dna_file.txt

# List all clinical variants in database
genomeinsight variants
```

### Python API

```python
from genomeinsight import load_dna_data, ClinicalAnalyzer

# Load your data
data = load_dna_data("your_dna_file.txt")
print(f"Loaded {data.snp_count:,} SNPs")

# Analyze clinical variants
analyzer = ClinicalAnalyzer()
result = analyzer.analyze(data)

# Check APOE status
if result.apoe_status:
    print(f"APOE: {result.apoe_status.genotype}")
    print(f"Risk: {result.apoe_status.risk_category}")

# View high-priority findings
for variant in result.high_priority:
    print(f"{variant.gene}: {variant.genotype} - {variant.risk_level.value}")
```

## 📊 Clinical Variant Categories

| Category | Description | Examples |
|----------|-------------|----------|
| ❤️ Cardiovascular | Heart and circulatory system | MTHFR, APOE, Factor V Leiden |
| 🔥 Metabolic | Diabetes and metabolic syndrome | TCF7L2, PPARG |
| 🎗️ Cancer | Cancer risk variants | TP53, NQO1, GSTP1 |
| 💊 Pharmacogenomics | Drug metabolism | CYP2D6, CYP2C19, VKORC1 |
| 🥗 Nutrient | Nutrient metabolism | LCT (lactose), HFE (iron), VDR (vitamin D) |
| 🧠 Neurological | Mental health and cognition | BDNF, COMT, DRD2 |

## 🤖 AI-Powered Reports

Generate natural language reports using OpenAI or Anthropic:

```bash
# Set your API key
export OPENAI_API_KEY=sk-...  # or ANTHROPIC_API_KEY

# Generate consumer-friendly report
genomeinsight analyze your_dna.txt --ai --ai-style consumer

# Generate technical report for researchers
genomeinsight analyze your_dna.txt --ai --ai-style technical

# Save report to file
genomeinsight analyze your_dna.txt --ai --ai-output report.md
```

| Style | Audience | Features |
|-------|----------|----------|
| `consumer` | General public | Plain language, lifestyle tips, encouraging tone |
| `technical` | Researchers/clinicians | Detailed analysis, citations, clinical considerations |
| `both` | All | Generates both report types |

## 🔬 Gene Interactions

GenomeInsight detects known gene-gene interactions:

- **MTHFR Compound Heterozygosity**: C677T + A1298C significantly reduces enzyme function
- **APOE Isoforms**: Combines rs429358 + rs7412 for accurate ε2/ε3/ε4 calling
- **CYP Compound Effects**: Multiple CYP variants affecting drug metabolism

## 📁 Project Structure

```
genomeinsight/
├── core/              # Data loaders and utilities
├── clinical/          # Clinical variant analysis
├── pharmacogenomics/  # Drug metabolism analysis
├── polygenic/         # Polygenic risk score calculator ✅
├── ai/                # AI-powered report generation ✅
├── ancestry/          # Ancestry composition estimation ✅
├── reports/           # HTML and JSON report generation
└── cli.py             # Command-line interface

tests/                 # Comprehensive test suite (163 tests) ✅
examples/              # Sample PRS weight files
```

## 🛠️ Development

```bash
# Install with development dependencies
uv sync --all-extras

# Run tests
uv run pytest

# Run linting
uv run ruff check .

# Format code
uv run ruff format .

# Type checking
uv run mypy genomeinsight
```

## 🗺️ Roadmap

- [x] Polygenic Risk Score calculator
- [x] Comprehensive test suite (61+ tests)
- [x] AI-powered natural language reports (OpenAI & Anthropic)
- [x] Ancestry composition estimation
- [ ] REST API for integration
- [ ] Web dashboard interface
- [ ] VCF file support improvements

## ⚠️ Disclaimer

> **This tool is for educational and research purposes only.**

- Not intended to diagnose, treat, or prevent any disease
- Genetic variants have incomplete penetrance
- Environmental factors significantly influence health outcomes
- Always consult healthcare providers for medical decisions
- Consider professional genetic counseling for significant findings

## 📄 License

MIT License - see [LICENSE](LICENSE) for details.

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📚 References

- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) - Clinical variant database
- [PharmGKB](https://www.pharmgkb.org/) - Pharmacogenomics knowledge base
- [CPIC](https://cpicpgx.org/) - Clinical pharmacogenetics implementation
- [SNPedia](https://www.snpedia.com/) - SNP wiki
- [GWAS Catalog](https://www.ebi.ac.uk/gwas/) - Genome-wide association studies

---

<p align="center">
  Made with ❤️ for the personal genomics community
</p>
