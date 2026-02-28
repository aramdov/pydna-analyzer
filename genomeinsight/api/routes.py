"""API route definitions."""

from __future__ import annotations

from fastapi import APIRouter, Depends, HTTPException

from genomeinsight import __version__
from genomeinsight.api.schemas import ResolvedFile, cleanup_temp, resolve_file_input

router = APIRouter()


@router.get("/health")
def health():
    """Health check endpoint."""
    return {"status": "ok", "version": __version__}


@router.post("/info")
async def info(resolved: ResolvedFile = Depends(resolve_file_input)):  # noqa: B008
    """Return basic information about a DNA data file."""
    from genomeinsight.core.data_loader import load_dna_data

    try:
        dataset = load_dna_data(resolved.path)
        return {
            "format": dataset.format.value,
            "source_file": dataset.source_file,
            "snp_count": dataset.snp_count,
        }
    except Exception as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    finally:
        cleanup_temp(resolved)


@router.get("/variants")
def variants():
    """Return all clinical variants from the database."""
    from genomeinsight.clinical.variants import CLINICAL_VARIANTS

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
