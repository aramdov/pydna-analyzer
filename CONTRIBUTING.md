# Contributing to GenomeInsight

Thank you for your interest in contributing to GenomeInsight! 🧬

## Getting Started

1. Fork the repository
2. Clone your fork: `git clone https://github.com/YOUR_USERNAME/genomeinsight.git`
3. Install with development dependencies: `uv sync --all-extras`
4. Create a feature branch: `git checkout -b feature/your-feature-name`

## Development Workflow

### Running Tests

```bash
uv run pytest
uv run pytest --cov=genomeinsight  # With coverage
```

### Code Quality

```bash
uv run ruff check .      # Linting
uv run ruff format .     # Formatting
uv run mypy genomeinsight  # Type checking
```

### Pre-commit Hooks (Optional)

```bash
uv run pre-commit install
```

## Guidelines

### Code Style

- Follow PEP 8, enforced by Ruff
- Use type hints for all function signatures
- Write docstrings for public functions and classes

### Commits

- Use clear, descriptive commit messages
- Reference issues when applicable: `Fix #123: Description`

### Pull Requests

- Keep PRs focused on a single feature or fix
- Update documentation if needed
- Add tests for new functionality
- Ensure all tests pass before submitting

## Adding Clinical Variants

When adding new SNPs to the clinical database:

1. Provide peer-reviewed sources (PubMed IDs preferred)
2. Include rsID, gene name, and risk allele
3. Add appropriate category classification
4. Document the clinical significance

## Questions?

Open an issue for any questions or discussions!
