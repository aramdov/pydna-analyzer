"""API route definitions."""

from __future__ import annotations

import dataclasses
from enum import Enum

from fastapi import APIRouter, Depends, HTTPException, Query
from fastapi.responses import HTMLResponse

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


def _serialize_enum(value: object) -> object:
    """Convert an Enum to its .value; pass through anything else."""
    return value.value if isinstance(value, Enum) else value


def _serialize_analysis_result(result) -> dict:
    """Serialize an AnalysisResult dataclass to a JSON-safe dict."""
    from genomeinsight.clinical.apoe import APOEResult

    def _serialize_variant(vr) -> dict:
        d = dataclasses.asdict(vr)
        for key, val in d.items():
            d[key] = _serialize_enum(val)
        return d

    variant_results = [_serialize_variant(vr) for vr in result.variant_results]
    by_category = {
        _serialize_enum(cat): [_serialize_variant(vr) for vr in vrs]
        for cat, vrs in result.by_category.items()
    }

    return {
        "source_file": result.source_file,
        "snp_count": result.snp_count,
        "variants_analyzed": result.variants_analyzed,
        "variants_found": result.variants_found,
        "apoe_status": (
            result.apoe_status.model_dump(mode="json")
            if isinstance(result.apoe_status, APOEResult)
            else result.apoe_status
        ),
        "variant_results": variant_results,
        "high_priority": [_serialize_variant(vr) for vr in result.high_priority],
        "moderate_priority": [_serialize_variant(vr) for vr in result.moderate_priority],
        "low_priority": [_serialize_variant(vr) for vr in result.low_priority],
        "by_category": by_category,
        "gene_interactions": result.gene_interactions,
    }


@router.post("/analyze")
async def analyze(
    resolved: ResolvedFile = Depends(resolve_file_input),  # noqa: B008
    format: str = Query("json", alias="format"),
):
    """Analyze a DNA file for clinical variants."""
    from genomeinsight.clinical.analyzer import ClinicalAnalyzer
    from genomeinsight.core.data_loader import load_dna_data

    try:
        dataset = load_dna_data(resolved.path)
        analyzer = ClinicalAnalyzer()
        result = analyzer.analyze(dataset)

        if format == "html":
            import tempfile
            from pathlib import Path

            from genomeinsight.reports.html_report import generate_html_report

            with tempfile.NamedTemporaryFile(
                delete=False, suffix=".html"
            ) as tmp:
                tmp_path = Path(tmp.name)
            generate_html_report(result, tmp_path)
            html_content = tmp_path.read_text()
            tmp_path.unlink()
            return HTMLResponse(content=html_content)

        return _serialize_analysis_result(result)
    except Exception as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    finally:
        cleanup_temp(resolved)
