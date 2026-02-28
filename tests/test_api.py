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
