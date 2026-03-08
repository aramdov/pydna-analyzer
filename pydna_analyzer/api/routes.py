"""API route definitions."""

from __future__ import annotations

import dataclasses
from enum import Enum

from fastapi import APIRouter, Depends, HTTPException, Query
from fastapi.responses import HTMLResponse

from pydna_analyzer import __version__
from pydna_analyzer.api.schemas import (
    ResolvedFile,
    cleanup_temp,
    resolve_file_input,
    resolve_weights_input,
)

router = APIRouter()


@router.get("/health")
def health():
    """Health check endpoint."""
    return {"status": "ok", "version": __version__}


@router.post("/info")
async def info(resolved: ResolvedFile = Depends(resolve_file_input)):  # noqa: B008
    """Return basic information about a DNA data file."""
    from pydna_analyzer.core.data_loader import load_dna_data

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


def _deep_serialize(obj: object) -> object:
    """Recursively convert Enums to .value in nested dicts/lists."""
    if isinstance(obj, Enum):
        return obj.value
    if isinstance(obj, dict):
        return {_deep_serialize(k): _deep_serialize(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_deep_serialize(v) for v in obj]
    return obj


def _serialize_analysis_result(result) -> dict:
    """Serialize an AnalysisResult dataclass to a JSON-safe dict."""
    from pydna_analyzer.clinical.apoe import APOEResult

    def _serialize_variant(vr) -> dict:
        return _deep_serialize(dataclasses.asdict(vr))

    variant_results = [_serialize_variant(vr) for vr in result.variant_results]
    by_category = {
        (cat.value if isinstance(cat, Enum) else cat): [_serialize_variant(vr) for vr in vrs]
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
    from pydna_analyzer.clinical.analyzer import ClinicalAnalyzer
    from pydna_analyzer.core.data_loader import load_dna_data

    try:
        dataset = load_dna_data(resolved.path)
        analyzer = ClinicalAnalyzer()
        result = analyzer.analyze(dataset)

        if format == "html":
            import tempfile
            from pathlib import Path

            from pydna_analyzer.reports.html_report import generate_html_report

            tmp_path = None
            try:
                with tempfile.NamedTemporaryFile(
                    delete=False, suffix=".html"
                ) as tmp:
                    tmp_path = Path(tmp.name)
                generate_html_report(result, tmp_path)
                html_content = tmp_path.read_text()
                return HTMLResponse(content=html_content)
            finally:
                if tmp_path and tmp_path.exists():
                    tmp_path.unlink()

        return _serialize_analysis_result(result)
    except Exception as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    finally:
        cleanup_temp(resolved)


def _serialize_pgx_result(result) -> dict:
    """Serialize a PGxResult dataclass to a JSON-safe dict."""
    return {
        "gene_results": [
            {
                "gene": gr.gene,
                "diplotype": gr.diplotype,
                "phenotype": gr.phenotype.value,
                "activity_score": gr.activity_score,
                "confidence": gr.confidence,
                "snps_tested": gr.snps_tested,
                "snps_missing": gr.snps_missing,
                "interpretation": gr.interpretation,
                "drug_recommendations": [
                    {
                        "drug": dr.drug_name,
                        "drug_class": dr.drug_class,
                        "recommendation": dr.recommendation,
                        "cpic_level": dr.cpic_level.value,
                        "alternatives": dr.alternatives,
                    }
                    for dr in gr.drug_recommendations
                ],
            }
            for gr in result.gene_results
        ],
        "summary": result.summary,
        "actionable_count": result.actionable_count,
    }


@router.post("/pgx")
async def pgx(
    resolved: ResolvedFile = Depends(resolve_file_input),  # noqa: B008
    gene: str | None = Query(None),
):
    """Run pharmacogenomics analysis on a DNA file."""
    from pydna_analyzer.core.data_loader import load_dna_data
    from pydna_analyzer.pharmacogenomics import PGxAnalyzer, PGxResult

    try:
        dataset = load_dna_data(resolved.path)
        analyzer = PGxAnalyzer()

        if gene is not None:
            single_result = analyzer.analyze_gene(dataset, gene.upper())
            if single_result is None:
                available = ", ".join(sorted(analyzer.gene_definitions))
                raise HTTPException(
                    status_code=422,
                    detail=f"Unknown gene: {gene}. Available: {available}",
                )
            result = PGxResult(
                gene_results=[single_result],
                total_genes_tested=1,
                total_snps_tested=single_result.snps_tested,
                total_snps_missing=single_result.snps_missing,
            )
        else:
            result = analyzer.analyze(dataset)

        return _serialize_pgx_result(result)
    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    finally:
        cleanup_temp(resolved)


@router.post("/prs")
async def prs(
    resolved: ResolvedFile = Depends(resolve_file_input),  # noqa: B008
    weights_resolved: ResolvedFile = Depends(resolve_weights_input),  # noqa: B008
    name: str | None = None,
):
    """Calculate a Polygenic Risk Score from a DNA file and weights CSV."""
    from pydna_analyzer.core.data_loader import load_dna_data
    from pydna_analyzer.polygenic import PRSCalculator

    try:
        dataset = load_dna_data(resolved.path)
        calculator = PRSCalculator()
        result = calculator.calculate_from_file(dataset, weights_resolved.path, score_name=name)
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
    except Exception as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    finally:
        cleanup_temp(resolved)
        cleanup_temp(weights_resolved)


@router.post("/ancestry")
async def ancestry(
    resolved: ResolvedFile = Depends(resolve_file_input),  # noqa: B008
    bootstrap: int = Query(100, ge=10, le=10000),
):
    """Return an experimental ancestry approximation from a DNA file."""
    from pydna_analyzer.ancestry import AncestryAnalyzer
    from pydna_analyzer.core.data_loader import load_dna_data

    try:
        dataset = load_dna_data(resolved.path)
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
            "experimental": True,
            "disclaimer": (
                "This ancestry output is experimental and based on a simplified "
                "reference panel. It may be inaccurate and is not comparable to "
                "commercial ancestry services."
            ),
            "interpretation": result.interpretation,
        }
    except Exception as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    finally:
        cleanup_temp(resolved)
