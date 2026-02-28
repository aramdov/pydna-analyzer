"""Shared schemas and dependencies for the API."""

from __future__ import annotations

import tempfile
from dataclasses import dataclass
from pathlib import Path

from fastapi import Form, HTTPException, UploadFile


@dataclass
class ResolvedFile:
    """Tracks a resolved file path and whether it needs cleanup."""

    path: Path
    is_temp: bool


async def resolve_file_input(
    file: UploadFile | None = None,
    file_path: str | None = Form(None),
) -> ResolvedFile:
    """Resolve file input from either an upload or a local path.

    Used as a FastAPI dependency by all POST endpoints.
    """
    if file is not None:
        # Write upload to a temp file
        suffix = Path(file.filename).suffix if file.filename else ".txt"
        content = await file.read()
        with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
            tmp.write(content)
        return ResolvedFile(path=Path(tmp.name), is_temp=True)

    if file_path is not None:
        p = Path(file_path)
        if not p.exists():
            raise HTTPException(status_code=404, detail=f"File not found: {file_path}")
        return ResolvedFile(path=p, is_temp=False)

    raise HTTPException(
        status_code=422,
        detail="Provide either a file upload or a file_path form field.",
    )


def cleanup_temp(resolved: ResolvedFile) -> None:
    """Remove the temp file if it was created by resolve_file_input."""
    if resolved.is_temp and resolved.path.exists():
        resolved.path.unlink()
