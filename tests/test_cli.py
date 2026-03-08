"""Integration tests for the PyDNA Analyzer CLI."""

import json
from types import SimpleNamespace

import pytest
from typer.testing import CliRunner

from pydna_analyzer.cli import app

runner = CliRunner()


class TestAnalyzeCommand:
    """Tests for the 'analyze' CLI command."""

    def test_analyze_ancestrydna(self, sample_ancestrydna_file):
        """Analyze with AncestryDNA file should succeed and mention clinical variants."""
        result = runner.invoke(app, ["analyze", str(sample_ancestrydna_file)])
        assert result.exit_code == 0
        assert "clinical variants" in result.output.lower()

    def test_analyze_23andme(self, sample_23andme_file):
        """Analyze with 23andMe file should succeed and mention clinical variants."""
        result = runner.invoke(app, ["analyze", str(sample_23andme_file)])
        assert result.exit_code == 0
        assert "clinical variants" in result.output.lower()

    def test_analyze_json_output(self, sample_ancestrydna_file, tmp_path):
        """Analyze with --json flag should create a valid JSON output file."""
        output_file = tmp_path / "results.json"
        result = runner.invoke(
            app,
            ["analyze", str(sample_ancestrydna_file), "--json", "-o", str(output_file)],
        )
        assert result.exit_code == 0
        assert output_file.exists()
        # Verify it's valid JSON
        data = json.loads(output_file.read_text())
        assert "summary" in data
        assert "variants" in data

    def test_analyze_html_output(self, sample_ancestrydna_file, tmp_path):
        """Analyze with --html flag should create an HTML report file."""
        output_file = tmp_path / "report.html"
        result = runner.invoke(
            app,
            ["analyze", str(sample_ancestrydna_file), "--html", "-o", str(output_file)],
        )
        assert result.exit_code == 0
        assert output_file.exists()
        html_content = output_file.read_text()
        assert "<html" in html_content
        assert "PyDNA Analyzer" in html_content

    def test_analyze_nonexistent_file(self):
        """Analyze with a nonexistent file should exit with non-zero code."""
        result = runner.invoke(app, ["analyze", "/nonexistent/path/fake_dna.txt"])
        assert result.exit_code != 0

    def test_analyze_quiet_mode(self, sample_ancestrydna_file):
        """Analyze with --quiet should produce minimal output."""
        result_quiet = runner.invoke(
            app, ["analyze", str(sample_ancestrydna_file), "--quiet"]
        )
        result_normal = runner.invoke(
            app, ["analyze", str(sample_ancestrydna_file)]
        )
        assert result_quiet.exit_code == 0
        # Quiet mode output should be shorter than normal mode
        assert len(result_quiet.output) < len(result_normal.output)

    def test_analyze_ai_model_is_forwarded(
        self, sample_ancestrydna_file, monkeypatch, tmp_path
    ):
        """Analyze forwards --ai-model to the OpenAI client path."""
        captured: dict[str, str | None] = {}

        def fake_get_client(provider=None, api_key=None, model=None):
            captured["model"] = model
            captured["provider"] = None if provider is None else provider.value
            return SimpleNamespace(provider=SimpleNamespace(value="openai"), model=model)

        class FakeGenerator:
            def __init__(self, client):
                self.client = client

            def generate(self, result, style):
                return f"report for {style.value}"

        monkeypatch.setattr("pydna_analyzer.ai.get_client", fake_get_client)
        monkeypatch.setattr("pydna_analyzer.ai.AIReportGenerator", FakeGenerator)

        output_file = tmp_path / "ai-report.md"
        result = runner.invoke(
            app,
            [
                "analyze",
                str(sample_ancestrydna_file),
                "--ai",
                "--ai-provider",
                "openai",
                "--ai-model",
                "gpt-5-mini",
                "--ai-output",
                str(output_file),
            ],
        )

        assert result.exit_code == 0
        assert output_file.exists()
        assert captured["provider"] == "openai"
        assert captured["model"] == "gpt-5-mini"


class TestInfoCommand:
    """Tests for the 'info' CLI command."""

    def test_info_shows_file_details(self, sample_ancestrydna_file):
        """Info should display SNP count and detected format."""
        result = runner.invoke(app, ["info", str(sample_ancestrydna_file)])
        assert result.exit_code == 0
        output_lower = result.output.lower()
        assert "snp" in output_lower
        # The format should appear (ancestrydna)
        assert "ancestrydna" in output_lower or "format" in output_lower

    def test_info_nonexistent_file(self):
        """Info with a nonexistent file should exit with non-zero code."""
        result = runner.invoke(app, ["info", "/nonexistent/path/fake_dna.txt"])
        assert result.exit_code != 0


class TestVariantsCommand:
    """Tests for the 'variants' CLI command."""

    def test_variants_lists_database(self):
        """Variants command should output a table containing rsIDs."""
        result = runner.invoke(app, ["variants"])
        assert result.exit_code == 0
        # Should contain at least one known rsID from the database
        assert "rs" in result.output

    def test_variants_shows_count(self):
        """Variants command should display a total count."""
        result = runner.invoke(app, ["variants"])
        assert result.exit_code == 0
        assert "Total:" in result.output


class TestPrsCommand:
    """Tests for the 'prs' CLI command."""

    @pytest.fixture
    def weights_csv(self, tmp_path):
        """Create a minimal weights CSV for PRS testing."""
        content = (
            "rsid,effect_allele,weight\n"
            "rs429358,C,0.15\n"
            "rs7412,T,0.10\n"
            "rs1801133,T,0.05\n"
        )
        weights_file = tmp_path / "test_weights.csv"
        weights_file.write_text(content)
        return weights_file

    def test_prs_with_valid_weights(self, sample_ancestrydna_file, weights_csv):
        """PRS command with valid DNA file and weights should succeed."""
        result = runner.invoke(
            app,
            ["prs", str(sample_ancestrydna_file), "--weights", str(weights_csv)],
        )
        assert result.exit_code == 0
        output_lower = result.output.lower()
        # Should display score-related information
        assert "score" in output_lower or "raw score" in output_lower
        assert "snps used" in output_lower or "coverage" in output_lower

    def test_prs_missing_weights(self, sample_ancestrydna_file):
        """PRS command without --weights should fail."""
        result = runner.invoke(app, ["prs", str(sample_ancestrydna_file)])
        assert result.exit_code != 0


class TestVersionFlag:
    """Tests for the --version flag."""

    def test_version_flag(self):
        """--version should display the version string."""
        result = runner.invoke(app, ["--version"])
        assert result.exit_code == 0
        assert "version" in result.output.lower()
