"""Tests for the GenomeInsight REST API."""

from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from genomeinsight.api import create_app


@pytest.fixture
def client():
    """Create a test client for the API."""
    app = create_app()
    return TestClient(app)


class TestHealthEndpoint:
    def test_health_returns_ok(self, client):
        response = client.get("/health")
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "ok"

    def test_health_includes_version(self, client):
        response = client.get("/health")
        data = response.json()
        assert "version" in data
        assert data["version"] == "0.1.0"


class TestFileInput:
    def test_info_with_file_path(self, client, sample_ancestrydna_file):
        """Path-based input works."""
        response = client.post("/info", data={"file_path": str(sample_ancestrydna_file)})
        assert response.status_code == 200
        data = response.json()
        assert data["snp_count"] > 0

    def test_info_with_upload(self, client, sample_ancestrydna_content):
        """Upload-based input works."""
        response = client.post(
            "/info",
            files={"file": ("test.txt", sample_ancestrydna_content, "text/plain")},
        )
        assert response.status_code == 200
        data = response.json()
        assert data["snp_count"] > 0

    def test_info_no_input_returns_422(self, client):
        """Neither path nor upload -> 422."""
        response = client.post("/info")
        assert response.status_code == 422

    def test_info_bad_path_returns_404(self, client):
        """Nonexistent path -> 404."""
        response = client.post("/info", data={"file_path": "/nonexistent/file.txt"})
        assert response.status_code == 404


class TestVariantsEndpoint:
    def test_variants_returns_list(self, client):
        response = client.get("/variants")
        assert response.status_code == 200
        data = response.json()
        assert isinstance(data["variants"], list)
        assert len(data["variants"]) > 30

    def test_variant_has_required_fields(self, client):
        response = client.get("/variants")
        variant = response.json()["variants"][0]
        assert "rsid" in variant
        assert "gene" in variant
        assert "name" in variant
        assert "category" in variant
        assert "evidence" in variant


class TestAnalyzeEndpoint:
    def test_analyze_returns_results(self, client, sample_ancestrydna_file):
        response = client.post("/analyze", data={"file_path": str(sample_ancestrydna_file)})
        assert response.status_code == 200
        data = response.json()
        assert "snp_count" in data
        assert "variants_found" in data
        assert "variant_results" in data

    def test_analyze_with_upload(self, client, sample_ancestrydna_content):
        response = client.post(
            "/analyze",
            files={"file": ("test.txt", sample_ancestrydna_content, "text/plain")},
        )
        assert response.status_code == 200
        assert response.json()["snp_count"] > 0

    def test_analyze_html_format(self, client, sample_ancestrydna_file):
        response = client.post(
            "/analyze",
            data={"file_path": str(sample_ancestrydna_file)},
            params={"format": "html"},
        )
        assert response.status_code == 200
        assert "text/html" in response.headers["content-type"]

    def test_analyze_includes_apoe(self, client, sample_ancestrydna_file):
        response = client.post("/analyze", data={"file_path": str(sample_ancestrydna_file)})
        data = response.json()
        assert "apoe_status" in data


class TestPgxEndpoint:
    def test_pgx_returns_results(self, client, tmp_path, pgx_all_wildtype_content):
        filepath = tmp_path / "pgx.txt"
        filepath.write_text(pgx_all_wildtype_content)
        response = client.post("/pgx", data={"file_path": str(filepath)})
        assert response.status_code == 200
        data = response.json()
        assert "gene_results" in data
        assert len(data["gene_results"]) > 0

    def test_pgx_single_gene(self, client, tmp_path, pgx_all_wildtype_content):
        filepath = tmp_path / "pgx.txt"
        filepath.write_text(pgx_all_wildtype_content)
        response = client.post(
            "/pgx", data={"file_path": str(filepath)}, params={"gene": "CYP2C19"}
        )
        assert response.status_code == 200
        data = response.json()
        assert len(data["gene_results"]) == 1
        assert data["gene_results"][0]["gene"] == "CYP2C19"

    def test_pgx_unknown_gene(self, client, tmp_path, pgx_all_wildtype_content):
        filepath = tmp_path / "pgx.txt"
        filepath.write_text(pgx_all_wildtype_content)
        response = client.post(
            "/pgx", data={"file_path": str(filepath)}, params={"gene": "FAKEGENE"}
        )
        assert response.status_code == 422


@pytest.fixture
def sample_weights_file(tmp_path):
    """Create a sample PRS weights CSV file."""
    content = "rsid,effect_allele,weight\nrs429358,C,0.15\nrs7412,T,0.10\nrs1801133,T,0.05\n"
    filepath = tmp_path / "weights.csv"
    filepath.write_text(content)
    return filepath


class TestPrsEndpoint:
    def test_prs_with_paths(self, client, sample_ancestrydna_file, sample_weights_file):
        response = client.post(
            "/prs",
            data={
                "file_path": str(sample_ancestrydna_file),
                "weights_path": str(sample_weights_file),
            },
        )
        assert response.status_code == 200
        data = response.json()
        assert "raw_score" in data
        assert "snps_used" in data

    def test_prs_missing_weights_returns_422(self, client, sample_ancestrydna_file):
        response = client.post(
            "/prs",
            data={"file_path": str(sample_ancestrydna_file)},
        )
        assert response.status_code == 422

    def test_prs_with_name(self, client, sample_ancestrydna_file, sample_weights_file):
        response = client.post(
            "/prs",
            data={
                "file_path": str(sample_ancestrydna_file),
                "weights_path": str(sample_weights_file),
            },
            params={"name": "Test Score"},
        )
        data = response.json()
        assert data["score_name"] == "Test Score"


class TestAncestryEndpoint:
    def test_ancestry_returns_results(self, client, tmp_path, sample_ancestrydna_content):
        filepath = tmp_path / "ancestry.txt"
        filepath.write_text(sample_ancestrydna_content)
        response = client.post("/ancestry", data={"file_path": str(filepath)})
        assert response.status_code == 200
        data = response.json()
        assert "populations" in data
        assert "snps_used" in data
        assert "interpretation" in data

    def test_ancestry_bootstrap_param(self, client, tmp_path, sample_ancestrydna_content):
        filepath = tmp_path / "ancestry.txt"
        filepath.write_text(sample_ancestrydna_content)
        response = client.post(
            "/ancestry", data={"file_path": str(filepath)}, params={"bootstrap": 50}
        )
        assert response.status_code == 200

    def test_ancestry_with_upload(self, client, sample_ancestrydna_content):
        response = client.post(
            "/ancestry",
            files={"file": ("test.txt", sample_ancestrydna_content, "text/plain")},
        )
        assert response.status_code == 200


class TestServeCLI:
    def test_serve_command_exists(self):
        from typer.testing import CliRunner

        from genomeinsight.cli import app as cli_app

        runner = CliRunner()
        result = runner.invoke(cli_app, ["serve", "--help"])
        assert result.exit_code == 0
        assert "host" in result.output.lower()
        assert "port" in result.output.lower()


class TestErrorHandling:
    def test_analyze_corrupt_file(self, client, tmp_path):
        """Corrupt file returns 422, not 500."""
        filepath = tmp_path / "corrupt.txt"
        filepath.write_text("this is not a valid DNA file at all")
        response = client.post("/analyze", data={"file_path": str(filepath)})
        assert response.status_code == 422
        assert "detail" in response.json()

    def test_pgx_corrupt_file(self, client, tmp_path):
        filepath = tmp_path / "corrupt.txt"
        filepath.write_text("garbage data")
        response = client.post("/pgx", data={"file_path": str(filepath)})
        assert response.status_code == 422

    def test_ancestry_corrupt_file(self, client, tmp_path):
        filepath = tmp_path / "corrupt.txt"
        filepath.write_text("not dna data")
        response = client.post("/ancestry", data={"file_path": str(filepath)})
        assert response.status_code == 422
