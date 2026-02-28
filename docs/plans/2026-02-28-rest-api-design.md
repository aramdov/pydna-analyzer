# REST API Design — GenomeInsight

**Date**: 2026-02-28
**Status**: Approved

## Overview

Add a local-first REST API to GenomeInsight using FastAPI. The API wraps existing analysis
pipelines behind HTTP endpoints, enabling integration with local AI agents, MCP servers,
scripts, and a future web dashboard. All processing remains local — consistent with the
project's privacy-first philosophy.

## Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Use case | Local-first server | Privacy-first; no auth/multi-user complexity |
| Framework | FastAPI | Pydantic-native, auto OpenAPI docs, async threadpool |
| File input | Path + upload | Path for agents/scripts, upload for dashboard/curl |
| Endpoints | All 6 CLI commands | Full CLI parity |
| Response format | JSON only | HTML via explicit `?format=html` on `/analyze` |

## Module Structure

```
genomeinsight/api/
├── __init__.py      # create_app() factory
├── routes.py        # All endpoint definitions
└── schemas.py       # Request/response Pydantic models
```

The API is a thin HTTP layer. No business logic lives here — every endpoint delegates to the
same analyzers the CLI uses.

## Endpoints

| Method | Path | Input | Output |
|--------|------|-------|--------|
| GET | `/health` | — | `{"status": "ok", "version": "..."}` |
| GET | `/variants` | — | List of all clinical variants |
| POST | `/info` | file (path or upload) | SNP count, format metadata |
| POST | `/analyze` | file, `?format=html` | Clinical analysis (JSON or HTML) |
| POST | `/pgx` | file, optional `gene` | Diplotypes, phenotypes, drug recs |
| POST | `/ancestry` | file, optional `bootstrap` | Population proportions with CIs |
| POST | `/prs` | file + weights file, optional `name` | PRS score, percentile |

### File Input

POST endpoints accept either:
- `file_path` form field (string) — server reads file directly (zero copy, local use)
- `file` upload via multipart/form-data — server saves to temp, processes, cleans up

When both provided, upload takes precedence. Neither provided → 422 error.

PRS endpoint additionally accepts `weights_path` or `weights_file` using the same pattern.

## Request/Response Flow

```
Client POST /analyze (file upload or path)
  → routes.py: resolve_file_input()
    → Save UploadFile to tempfile, OR validate path exists
    → Returns Path object
  → load_dna_data(path) → DNADataset
  → ClinicalAnalyzer().analyze(dataset) → AnalysisResult
  → Serialize → JSON response
  → Cleanup temp file if upload
```

## Error Handling

- `FileNotFoundError` → 404
- `ValueError` (bad format, invalid params) → 422
- Unexpected errors → 500 (generic message, no stack traces)

## Sync Endpoints

All endpoints use regular `def` (not `async def`). Analysis functions are CPU-bound
(numpy/scipy), not I/O-bound. FastAPI runs sync endpoints in a threadpool, preventing
event loop blocking during long computations like ancestry bootstrap.

## CLI Integration

New command:
```
genomeinsight serve [--host localhost] [--port 8000] [--reload]
```

Lazy imports uvicorn and create_app. Friendly error if FastAPI not installed.

## Dependencies

```toml
[project.optional-dependencies]
api = [
    "fastapi>=0.110",
    "uvicorn[standard]>=0.27",
    "python-multipart>=0.0.9",
]
```

`all` extra updated to include `api`.

## CORS

Allow `localhost` origins by default (for future web dashboard). Configurable via
`create_app(cors_origins=[...])`.

## Test Strategy

New `tests/test_api.py` using FastAPI's `TestClient` (~20-25 tests):

- **Health/meta**: GET /health, GET /variants
- **File input**: path works, upload works, missing → 422, bad path → 404
- **Analyze**: valid results, HTML format flag, schema structure
- **PGx**: full analysis, single gene, unknown gene → 422
- **Ancestry**: populations with CIs, bootstrap param
- **PRS**: dual file input, missing weights → 422
- **Errors**: bad format, corrupt file, graceful responses

Reuses existing conftest.py fixtures plus new `api_client` fixture.
