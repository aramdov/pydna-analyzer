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


async def _resolve_input(
    upload: UploadFile | None,
    path_str: str | None,
    default_suffix: str,
    label: str,
) -> ResolvedFile:
    """Shared logic for resolving file input from upload or local path."""
    if upload is not None:
        suffix = Path(upload.filename).suffix if upload.filename else default_suffix
        content = await upload.read()
        with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
            tmp.write(content)
        return ResolvedFile(path=Path(tmp.name), is_temp=True)

    if path_str is not None:
        p = Path(path_str).resolve()
        if not p.exists():
            raise HTTPException(status_code=404, detail=f"File not found: {path_str}")
        if not p.is_file():
            raise HTTPException(status_code=422, detail=f"Not a regular file: {path_str}")
        return ResolvedFile(path=p, is_temp=False)

    raise HTTPException(
        status_code=422,
        detail=f"Provide either a {label} upload or a {label}_path form field.",
    )


async def resolve_file_input(
    file: UploadFile | None = None,
    file_path: str | None = Form(None),
) -> ResolvedFile:
    """Resolve DNA file input from either an upload or a local path."""
    return await _resolve_input(file, file_path, ".txt", "file")


async def resolve_weights_input(
    weights_file: UploadFile | None = None,
    weights_path: str | None = Form(None),
) -> ResolvedFile:
    """Resolve weights file input from either an upload or a local path."""
    return await _resolve_input(weights_file, weights_path, ".csv", "weights")


def cleanup_temp(resolved: ResolvedFile) -> None:
    """Remove the temp file if it was created during upload resolution."""
    if resolved.is_temp and resolved.path.exists():
        resolved.path.unlink()
