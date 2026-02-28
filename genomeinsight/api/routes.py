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
    from genomeinsight.core.data_loader import load_dna_data
    from genomeinsight.pharmacogenomics import PGxAnalyzer, PGxResult

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


@router.post("/ancestry")
async def ancestry(
    resolved: ResolvedFile = Depends(resolve_file_input),  # noqa: B008
    bootstrap: int = Query(100),
):
    """Estimate ancestry composition from a DNA file."""
    from genomeinsight.ancestry import AncestryAnalyzer
    from genomeinsight.core.data_loader import load_dna_data

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
            "interpretation": result.interpretation,
        }
    except Exception as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    finally:
        cleanup_temp(resolved)
