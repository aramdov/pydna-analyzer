# REST API Implementation Plan

**Status:** Complete (merged to main 2026-03-01) — 11 tasks, 25 tests, all passing.

**Goal:** Add a local-first FastAPI REST API with full CLI parity (all 6 analysis commands as HTTP endpoints).

**Architecture:** Thin HTTP layer in `pydna_analyzer/api/` that delegates to existing analyzers. File input via local path or multipart upload. All endpoints return JSON; `/analyze` supports `?format=html`. Sync endpoints (CPU-bound work runs in FastAPI's threadpool).

**Tech Stack:** FastAPI, uvicorn, python-multipart, httpx (test client via FastAPI's TestClient)

---

### Task 1: Add API Dependencies

**Files:**
- Modify: `pyproject.toml:60-73`

**Step 1: Add the `api` optional dependency group**

In `pyproject.toml`, after the `ancestry` group (~line 64), add:

```toml
api = [
    "fastapi>=0.110",
    "uvicorn[standard]>=0.27",
    "python-multipart>=0.0.9",
]
```

Update the `all` extra to include `api`:

```toml
all = [
    "pydna_analyzer[dev,prs,ancestry,ai,api]",
]
```

**Step 2: Install**

Run: `uv sync --all-extras`
Expected: FastAPI, uvicorn, python-multipart installed successfully.

**Step 3: Verify import**

Run: `uv run python -c "import fastapi; print(fastapi.__version__)"`
Expected: Prints version >= 0.110

**Step 4: Commit**

```bash
git add pyproject.toml
git commit -m "feat(api): add FastAPI optional dependencies"
```

---

### Task 2: Create App Factory and Health Endpoint

**Files:**
- Create: `pydna_analyzer/api/__init__.py`
- Create: `pydna_analyzer/api/routes.py`
- Create: `tests/test_api.py`

**Step 1: Write failing tests for health endpoint**

Create `tests/test_api.py`:

```python
"""Tests for the PyDNA Analyzer REST API."""

from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from pydna_analyzer.api import create_app


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
```

**Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_api.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'pydna_analyzer.api'`

**Step 3: Implement app factory and health route**

Create `pydna_analyzer/api/__init__.py`:

```python
"""PyDNA Analyzer REST API."""

from __future__ import annotations

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from pydna_analyzer import __version__


def create_app(cors_origins: list[str] | None = None) -> FastAPI:
    """Create and configure the FastAPI application."""
    app = FastAPI(
        title="PyDNA Analyzer API",
        description="Privacy-first personal genomics analysis API",
        version=__version__,
    )

    # CORS for local dashboard
    origins = cors_origins or ["http://localhost:3000", "http://localhost:5173"]
    app.add_middleware(
        CORSMiddleware,
        allow_origins=origins,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    from pydna_analyzer.api.routes import router

    app.include_router(router)

    return app
```

Create `pydna_analyzer/api/routes.py`:

```python
"""API route definitions."""

from __future__ import annotations

from fastapi import APIRouter

from pydna_analyzer import __version__

router = APIRouter()


@router.get("/health")
def health():
    """Health check endpoint."""
    return {"status": "ok", "version": __version__}
```

**Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_api.py -v`
Expected: 2 passed

**Step 5: Commit**

```bash
git add pydna_analyzer/api/__init__.py pydna_analyzer/api/routes.py tests/test_api.py
git commit -m "feat(api): add app factory and health endpoint"
```

---

### Task 3: File Input Resolution Helper

**Files:**
- Create: `pydna_analyzer/api/schemas.py`
- Modify: `pydna_analyzer/api/routes.py`
- Modify: `tests/test_api.py`

This task builds the shared `resolve_file_input()` function that all POST endpoints use.

**Step 1: Write failing tests for file input**

Add to `tests/test_api.py`:

```python
from pathlib import Path
import io


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
        """Neither path nor upload → 422."""
        response = client.post("/info")
        assert response.status_code == 422

    def test_info_bad_path_returns_404(self, client):
        """Nonexistent path → 404."""
        response = client.post("/info", data={"file_path": "/nonexistent/file.txt"})
        assert response.status_code == 404
```

**Step 2: Run to verify they fail**

Run: `uv run pytest tests/test_api.py::TestFileInput -v`
Expected: FAIL — no `/info` endpoint

**Step 3: Implement schemas and file resolution**

Create `pydna_analyzer/api/schemas.py`:

```python
"""API request/response schemas."""

from __future__ import annotations

import tempfile
from pathlib import Path

from fastapi import Form, HTTPException, UploadFile


async def resolve_file_input(
    file: UploadFile | None = None,
    file_path: str | None = Form(None),
) -> Path:
    """Resolve file input from upload or local path.

    Returns the path to the file. For uploads, writes to a temp file.
    Caller is responsible for cleanup of temp files.
    """
    if file is not None:
        suffix = Path(file.filename).suffix if file.filename else ".txt"
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
        content = await file.read()
        tmp.write(content)
        tmp.close()
        return Path(tmp.name)

    if file_path is not None:
        path = Path(file_path)
        if not path.exists():
            raise HTTPException(status_code=404, detail=f"File not found: {file_path}")
        return path

    raise HTTPException(status_code=422, detail="Provide either 'file' upload or 'file_path'")
```

**Step 4: Add the `/info` endpoint to routes.py**

Add to `pydna_analyzer/api/routes.py`:

```python
import os
from pathlib import Path

from fastapi import APIRouter, Depends, HTTPException, UploadFile

from pydna_analyzer import __version__
from pydna_analyzer.api.schemas import resolve_file_input

router = APIRouter()


@router.get("/health")
def health():
    return {"status": "ok", "version": __version__}


@router.post("/info")
def info(resolved_path: Path = Depends(resolve_file_input)):
    """Get information about a DNA data file."""
    tmp_file = not str(resolved_path).startswith("/") or "tmp" in str(resolved_path)
    try:
        from pydna_analyzer.core.data_loader import load_dna_data

        dataset = load_dna_data(resolved_path)
        return {
            "format": dataset.format.value,
            "source_file": dataset.source_file,
            "snp_count": dataset.snp_count,
        }
    except Exception as e:
        raise HTTPException(status_code=422, detail=str(e))
    finally:
        # Clean up temp files from uploads
        if resolved_path.name.startswith("tmp") and resolved_path.exists():
            os.unlink(resolved_path)
```

Note: The temp file cleanup logic checks if the path looks like a tempfile. A cleaner approach is to track whether the file was uploaded. Implement this by returning a tuple from `resolve_file_input` or setting an attribute. The implementer should use the approach that keeps the code simplest — e.g., a small dataclass:

```python
@dataclass
class ResolvedFile:
    path: Path
    is_temp: bool
```

**Step 5: Run tests**

Run: `uv run pytest tests/test_api.py -v`
Expected: All 6 tests pass (2 health + 4 file input)

**Step 6: Commit**

```bash
git add pydna_analyzer/api/schemas.py pydna_analyzer/api/routes.py tests/test_api.py
git commit -m "feat(api): add file input resolution and /info endpoint"
```

---

### Task 4: Variants Endpoint (GET)

**Files:**
- Modify: `pydna_analyzer/api/routes.py`
- Modify: `tests/test_api.py`

**Step 1: Write failing test**

Add to `tests/test_api.py`:

```python
class TestVariantsEndpoint:
    def test_variants_returns_list(self, client):
        response = client.get("/variants")
        assert response.status_code == 200
        data = response.json()
        assert isinstance(data["variants"], list)
        assert len(data["variants"]) > 30  # We have 35+ curated SNPs

    def test_variant_has_required_fields(self, client):
        response = client.get("/variants")
        variant = response.json()["variants"][0]
        assert "rsid" in variant
        assert "gene" in variant
        assert "name" in variant
        assert "category" in variant
        assert "evidence" in variant
```

**Step 2: Run to verify fail**

Run: `uv run pytest tests/test_api.py::TestVariantsEndpoint -v`
Expected: FAIL — 404 or no route

**Step 3: Implement**

Add to `pydna_analyzer/api/routes.py`:

```python
@router.get("/variants")
def variants():
    """List all clinical variants in the database."""
    from pydna_analyzer.clinical.variants import CLINICAL_VARIANTS

    return {
        "variants": [
            {
                "rsid": rsid,
                "gene": v.gene,
                "name": v.name,
                "category": v.category.value,
                "evidence": v.evidence.value,
                "risk_allele": v.risk_allele,
            }
            for rsid, v in CLINICAL_VARIANTS.items()
        ],
        "total": len(CLINICAL_VARIANTS),
    }
```

**Step 4: Run tests**

Run: `uv run pytest tests/test_api.py -v`
Expected: All 8 tests pass

**Step 5: Commit**

```bash
git add pydna_analyzer/api/routes.py tests/test_api.py
git commit -m "feat(api): add GET /variants endpoint"
```

---

### Task 5: Analyze Endpoint

**Files:**
- Modify: `pydna_analyzer/api/routes.py`
- Modify: `tests/test_api.py`

**Step 1: Write failing tests**

Add to `tests/test_api.py`:

```python
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
```

**Step 2: Run to verify fail**

Run: `uv run pytest tests/test_api.py::TestAnalyzeEndpoint -v`
Expected: FAIL

**Step 3: Implement**

Add to `pydna_analyzer/api/routes.py`:

```python
from fastapi.responses import HTMLResponse


@router.post("/analyze")
def analyze(
    resolved_path: Path = Depends(resolve_file_input),
    format: str | None = None,
):
    """Analyze DNA data for clinical variants."""
    try:
        from pydna_analyzer.core.data_loader import load_dna_data
        from pydna_analyzer.clinical.analyzer import ClinicalAnalyzer

        dataset = load_dna_data(resolved_path)
        analyzer = ClinicalAnalyzer()
        result = analyzer.analyze(dataset)

        if format == "html":
            from pydna_analyzer.reports.html_report import generate_html_report
            import tempfile

            tmp_html = Path(tempfile.mktemp(suffix=".html"))
            generate_html_report(result, tmp_html)
            html_content = tmp_html.read_text()
            os.unlink(tmp_html)
            return HTMLResponse(content=html_content)

        return _serialize_analysis_result(result)
    except Exception as e:
        raise HTTPException(status_code=422, detail=str(e))
    finally:
        _cleanup_temp(resolved_path)
```

The `_serialize_analysis_result` helper converts the dataclass to a dict. The implementer should handle nested dataclasses (VariantResult, APOEResult) by converting enums to `.value` strings. Use `dataclasses.asdict()` as a starting point but handle enum serialization.

**Step 4: Run tests**

Run: `uv run pytest tests/test_api.py -v`
Expected: All 12 tests pass

**Step 5: Commit**

```bash
git add pydna_analyzer/api/routes.py tests/test_api.py
git commit -m "feat(api): add POST /analyze endpoint with HTML support"
```

---

### Task 6: PGx Endpoint

**Files:**
- Modify: `pydna_analyzer/api/routes.py`
- Modify: `tests/test_api.py`

**Step 1: Write failing tests**

Add to `tests/test_api.py`:

```python
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
            "/pgx",
            data={"file_path": str(filepath)},
            params={"gene": "CYP2C19"},
        )
        assert response.status_code == 200
        data = response.json()
        assert len(data["gene_results"]) == 1
        assert data["gene_results"][0]["gene"] == "CYP2C19"

    def test_pgx_unknown_gene(self, client, tmp_path, pgx_all_wildtype_content):
        filepath = tmp_path / "pgx.txt"
        filepath.write_text(pgx_all_wildtype_content)
        response = client.post(
            "/pgx",
            data={"file_path": str(filepath)},
            params={"gene": "FAKEGENE"},
        )
        assert response.status_code == 422
```

**Step 2: Run to verify fail**

Run: `uv run pytest tests/test_api.py::TestPgxEndpoint -v`
Expected: FAIL

**Step 3: Implement**

Add to `pydna_analyzer/api/routes.py`:

```python
@router.post("/pgx")
def pgx(
    resolved_path: Path = Depends(resolve_file_input),
    gene: str | None = None,
):
    """Pharmacogenomics analysis."""
    try:
        from pydna_analyzer.core.data_loader import load_dna_data
        from pydna_analyzer.pharmacogenomics import PGxAnalyzer, PGxResult

        dataset = load_dna_data(resolved_path)
        analyzer = PGxAnalyzer()

        if gene:
            single = analyzer.analyze_gene(dataset, gene.upper())
            if single is None:
                raise HTTPException(
                    status_code=422,
                    detail=f"Unknown gene: {gene}. Available: {', '.join(sorted(analyzer.gene_definitions))}",
                )
            result = PGxResult(
                gene_results=[single],
                total_genes_tested=1,
                total_snps_tested=single.snps_tested,
                total_snps_missing=single.snps_missing,
            )
        else:
            result = analyzer.analyze(dataset)

        return _serialize_pgx_result(result)
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=422, detail=str(e))
    finally:
        _cleanup_temp(resolved_path)
```

The `_serialize_pgx_result` helper should convert GeneResult dataclasses to dicts, including `drug_recommendations` with enum values as strings. Mirror the JSON structure from `cli.py:632-661`.

**Step 4: Run tests**

Run: `uv run pytest tests/test_api.py -v`
Expected: All 15 tests pass

**Step 5: Commit**

```bash
git add pydna_analyzer/api/routes.py tests/test_api.py
git commit -m "feat(api): add POST /pgx endpoint with gene filter"
```

---

### Task 7: Ancestry Endpoint

**Files:**
- Modify: `pydna_analyzer/api/routes.py`
- Modify: `tests/test_api.py`

**Step 1: Write failing tests**

Add to `tests/test_api.py`:

```python
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
            "/ancestry",
            data={"file_path": str(filepath)},
            params={"bootstrap": 50},
        )
        assert response.status_code == 200

    def test_ancestry_with_upload(self, client, sample_ancestrydna_content):
        response = client.post(
            "/ancestry",
            files={"file": ("test.txt", sample_ancestrydna_content, "text/plain")},
        )
        assert response.status_code == 200
```

**Step 2: Run to verify fail**

Run: `uv run pytest tests/test_api.py::TestAncestryEndpoint -v`
Expected: FAIL

**Step 3: Implement**

Add to `pydna_analyzer/api/routes.py`:

```python
@router.post("/ancestry")
def ancestry(
    resolved_path: Path = Depends(resolve_file_input),
    bootstrap: int = 100,
):
    """Estimate ancestry composition."""
    try:
        from pydna_analyzer.core.data_loader import load_dna_data
        from pydna_analyzer.ancestry import AncestryAnalyzer

        dataset = load_dna_data(resolved_path)
        analyzer = AncestryAnalyzer(n_bootstrap=bootstrap)
        result = analyzer.analyze(dataset)

        return {
            "populations": [
                {
                    "population": p.population,
                    "region": p.region,
                    "proportion": round(p.proportion, 4),
                    "confidence_low": round(p.confidence_low, 4),
                    "confidence_high": round(p.confidence_high, 4),
                }
                for p in result.populations
                if p.proportion >= 0.01
            ],
            "regions": {k: round(v, 4) for k, v in result.top_regions.items()},
            "snps_used": result.snps_used,
            "snps_available": result.snps_available,
            "coverage": round(result.coverage, 4),
            "convergence": result.convergence,
            "interpretation": result.interpretation,
        }
    except Exception as e:
        raise HTTPException(status_code=422, detail=str(e))
    finally:
        _cleanup_temp(resolved_path)
```

**Step 4: Run tests**

Run: `uv run pytest tests/test_api.py -v`
Expected: All 18 tests pass

**Step 5: Commit**

```bash
git add pydna_analyzer/api/routes.py tests/test_api.py
git commit -m "feat(api): add POST /ancestry endpoint"
```

---

### Task 8: PRS Endpoint (Dual File Input)

**Files:**
- Modify: `pydna_analyzer/api/schemas.py`
- Modify: `pydna_analyzer/api/routes.py`
- Modify: `tests/test_api.py`

This is the trickiest endpoint — it needs two files (DNA data + weights).

**Step 1: Write failing tests**

First, add a PRS weights fixture to `tests/test_api.py`:

```python
@pytest.fixture
def sample_weights_file(tmp_path):
    """Create a sample PRS weights CSV file."""
    content = "rsid,effect_allele,weight\nrs429358,C,0.15\nrs7412,T,0.10\nrs1801133,T,0.05\n"
    filepath = tmp_path / "weights.csv"
    filepath.write_text(content)
    return filepath
```

Then add tests:

```python
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
```

**Step 2: Run to verify fail**

Run: `uv run pytest tests/test_api.py::TestPrsEndpoint -v`
Expected: FAIL

**Step 3: Implement dual file resolution**

Add to `pydna_analyzer/api/schemas.py` a second resolver for the weights file:

```python
async def resolve_weights_input(
    weights_file: UploadFile | None = None,
    weights_path: str | None = Form(None),
) -> Path:
    """Resolve weights file input."""
    if weights_file is not None:
        suffix = Path(weights_file.filename).suffix if weights_file.filename else ".csv"
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
        content = await weights_file.read()
        tmp.write(content)
        tmp.close()
        return Path(tmp.name)

    if weights_path is not None:
        path = Path(weights_path)
        if not path.exists():
            raise HTTPException(status_code=404, detail=f"Weights file not found: {weights_path}")
        return path

    raise HTTPException(status_code=422, detail="Provide either 'weights_file' upload or 'weights_path'")
```

Add PRS endpoint to `pydna_analyzer/api/routes.py`:

```python
from pydna_analyzer.api.schemas import resolve_file_input, resolve_weights_input


@router.post("/prs")
def prs(
    resolved_path: Path = Depends(resolve_file_input),
    weights_resolved: Path = Depends(resolve_weights_input),
    name: str | None = None,
):
    """Calculate polygenic risk score."""
    try:
        from pydna_analyzer.core.data_loader import load_dna_data
        from pydna_analyzer.polygenic import PRSCalculator

        dataset = load_dna_data(resolved_path)
        calculator = PRSCalculator()
        result = calculator.calculate_from_file(dataset, weights_resolved, score_name=name)

        return {
            "score_name": result.score_name,
            "raw_score": result.raw_score,
            "normalized_score": result.normalized_score,
            "percentile": result.percentile,
            "risk_category": result.risk_category,
            "snps_used": result.snps_used,
            "snps_available": result.snps_available,
            "coverage": result.coverage,
            "interpretation": result.interpretation,
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=422, detail=str(e))
    finally:
        _cleanup_temp(resolved_path)
        _cleanup_temp(weights_resolved)
```

**Step 4: Run tests**

Run: `uv run pytest tests/test_api.py -v`
Expected: All 21 tests pass

**Step 5: Commit**

```bash
git add pydna_analyzer/api/schemas.py pydna_analyzer/api/routes.py tests/test_api.py
git commit -m "feat(api): add POST /prs endpoint with dual file input"
```

---

### Task 9: CLI `serve` Command

**Files:**
- Modify: `pydna_analyzer/cli.py`
- Modify: `tests/test_api.py`

**Step 1: Write failing test**

Add to `tests/test_api.py`:

```python
from typer.testing import CliRunner
from pydna_analyzer.cli import app as cli_app


class TestServeCLI:
    def test_serve_command_exists(self):
        runner = CliRunner()
        result = runner.invoke(cli_app, ["serve", "--help"])
        assert result.exit_code == 0
        assert "host" in result.output.lower()
        assert "port" in result.output.lower()
```

**Step 2: Run to verify fail**

Run: `uv run pytest tests/test_api.py::TestServeCLI -v`
Expected: FAIL — no `serve` command

**Step 3: Implement**

Add to `pydna_analyzer/cli.py` (after the ancestry command):

```python
@app.command()
def serve(
    host: str = typer.Option("127.0.0.1", "--host", "-h", help="Host to bind to"),
    port: int = typer.Option(8000, "--port", "-p", help="Port to bind to"),
    reload: bool = typer.Option(False, "--reload", help="Enable auto-reload for development"),
):
    """
    Start the PyDNA Analyzer REST API server.

    Examples:
        pydna_analyzer serve
        pydna_analyzer serve --port 9000
        pydna_analyzer serve --reload
    """
    try:
        import uvicorn
    except ImportError:
        console.print("[red]Missing dependency: uvicorn[/]")
        console.print("[dim]Install with: uv sync --extra api[/]")
        raise typer.Exit(1)

    try:
        from pydna_analyzer.api import create_app  # noqa: F401
    except ImportError:
        console.print("[red]Missing dependency: fastapi[/]")
        console.print("[dim]Install with: uv sync --extra api[/]")
        raise typer.Exit(1)

    console.print(f"[bold blue]🧬 PyDNA Analyzer API[/] starting on http://{host}:{port}")
    console.print("[dim]API docs: http://{}:{}/docs[/]".format(host, port))

    uvicorn.run(
        "pydna_analyzer.api:create_app",
        host=host,
        port=port,
        reload=reload,
        factory=True,
    )
```

**Step 4: Run tests**

Run: `uv run pytest tests/test_api.py -v`
Expected: All 22 tests pass

**Step 5: Commit**

```bash
git add pydna_analyzer/cli.py tests/test_api.py
git commit -m "feat(api): add 'serve' CLI command"
```

---

### Task 10: Error Handling and Edge Cases

**Files:**
- Modify: `pydna_analyzer/api/routes.py`
- Modify: `tests/test_api.py`

**Step 1: Write failing tests for error scenarios**

Add to `tests/test_api.py`:

```python
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
```

**Step 2: Run tests — some may already pass if error handling is solid**

Run: `uv run pytest tests/test_api.py::TestErrorHandling -v`

If any return 500 instead of 422, add a global exception handler in `pydna_analyzer/api/__init__.py`:

```python
from fastapi import Request
from fastapi.responses import JSONResponse

@app.exception_handler(Exception)
async def general_exception_handler(request: Request, exc: Exception):
    return JSONResponse(status_code=500, content={"detail": "Internal server error"})
```

**Step 3: Run all tests**

Run: `uv run pytest tests/test_api.py -v`
Expected: All ~25 tests pass

**Step 4: Run the full test suite to verify no regressions**

Run: `uv run pytest -v`
Expected: All ~188 tests pass (163 existing + ~25 new)

**Step 5: Commit**

```bash
git add pydna_analyzer/api/ tests/test_api.py
git commit -m "feat(api): add error handling and edge case tests"
```

---

### Task 11: Update Docs and Finalize

**Files:**
- Modify: `CLAUDE.md`
- Modify: `README.md`
- Modify: `pydna_analyzer/__init__.py`

**Step 1: Update CLAUDE.md**

Add to the Commands section:

```bash
# API Server
uv run pydna_analyzer serve                          # Start on localhost:8000
uv run pydna_analyzer serve --port 9000              # Custom port
uv run pydna_analyzer serve --reload                 # Dev mode with auto-reload
```

Add to Architecture → Module Map:

```
- **`api/`** — REST API (FastAPI):
  - `__init__.py` — App factory with CORS configuration
  - `routes.py` — Endpoint definitions: /health, /info, /analyze, /pgx, /ancestry, /prs, /variants
  - `schemas.py` — File input resolution (path or upload)
```

Update test count.

**Step 2: Update README**

Check the REST API roadmap item:

```
- [x] REST API for integration
```

Add API usage example to README.

**Step 3: Commit**

```bash
git add CLAUDE.md README.md
git commit -m "docs: update CLAUDE.md and README for REST API"
```

**Step 4: Run full suite one final time**

Run: `uv run pytest -v && uv run ruff check . && uv run ruff format --check .`
Expected: All pass, clean lint, clean format.
