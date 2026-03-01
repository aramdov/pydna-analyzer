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
