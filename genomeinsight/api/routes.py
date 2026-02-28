"""API route definitions."""

from __future__ import annotations

from fastapi import APIRouter

from genomeinsight import __version__

router = APIRouter()


@router.get("/health")
def health():
    """Health check endpoint."""
    return {"status": "ok", "version": __version__}
