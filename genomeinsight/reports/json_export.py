"""
JSON export for analysis results.
"""

from __future__ import annotations

import json
from pathlib import Path
from datetime import datetime
from dataclasses import asdict

from genomeinsight.clinical.analyzer import AnalysisResult


def export_to_json(result: AnalysisResult, output_path: Path):
    """
    Export analysis results to JSON.
    
    Args:
        result: AnalysisResult to export
        output_path: Path for output JSON file
    """
    from genomeinsight import __version__
    
    # Convert dataclasses to dicts
    def serialize_variant(v):
        return {
            "rsid": v.rsid,
            "gene": v.gene,
            "name": v.name,
            "category": v.category.value,
            "genotype": v.genotype,
            "risk_level": v.risk_level.value,
            "risk_allele": v.risk_allele,
            "description": v.description,
            "effect": v.effect,
            "evidence": v.evidence.value,
            "interventions": v.interventions,
            "chromosome": v.chromosome,
            "position": v.position,
            "has_interactions": v.has_interactions,
            "interaction_notes": v.interaction_notes,
        }
    
    output = {
        "metadata": {
            "version": __version__,
            "generated_at": datetime.now().isoformat(),
            "source_file": result.source_file,
        },
        "summary": {
            "snp_count": result.snp_count,
            "variants_analyzed": result.variants_analyzed,
            "variants_found": result.variants_found,
            "high_priority_count": len(result.high_priority),
            "moderate_priority_count": len(result.moderate_priority),
            "low_priority_count": len(result.low_priority),
            "gene_interactions_count": len(result.gene_interactions),
        },
        "apoe_status": None,
        "gene_interactions": result.gene_interactions,
        "variants": {
            "high_priority": [serialize_variant(v) for v in result.high_priority],
            "moderate_priority": [serialize_variant(v) for v in result.moderate_priority],
            "low_priority": [serialize_variant(v) for v in result.low_priority],
        },
        "by_category": {
            cat.value: [serialize_variant(v) for v in variants]
            for cat, variants in result.by_category.items()
        },
    }
    
    # Add APOE status
    if result.apoe_status:
        apoe = result.apoe_status
        output["apoe_status"] = {
            "genotype": apoe.genotype,
            "allele1": apoe.allele1.value,
            "allele2": apoe.allele2.value,
            "rs429358_genotype": apoe.rs429358_genotype,
            "rs7412_genotype": apoe.rs7412_genotype,
            "risk_category": apoe.risk_category,
            "interpretation": apoe.interpretation,
            "recommendations": apoe.recommendations,
            "has_e4": apoe.has_e4,
            "is_e4_homozygous": apoe.is_e4_homozygous,
            "has_e2": apoe.has_e2,
        }
    
    output_path = Path(output_path)
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
