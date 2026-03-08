"""
Clinical variant analyzer.

Analyzes DNA data for clinically significant variants and generates
annotated results with risk assessments and recommendations.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import pandas as pd

from pydna_analyzer.clinical.variants import (
    CLINICAL_VARIANTS,
    Category,
    ClinicalVariant,
    EvidenceLevel,
    RiskLevel,
    get_all_rsids,
    get_interacting_variants,
)
from pydna_analyzer.clinical.apoe import APOEResult, determine_apoe_from_dataset
from pydna_analyzer.core.data_loader import DNADataset


@dataclass
class VariantResult:
    """Result of analyzing a single variant."""
    
    rsid: str
    gene: str
    name: str
    category: Category
    genotype: str
    matched_genotype: Optional[str]
    risk_level: RiskLevel
    risk_allele: str
    description: str
    effect: str
    evidence: EvidenceLevel
    interventions: list[str]
    chromosome: str
    position: int
    
    # Gene interaction info
    has_interactions: bool = False
    interaction_notes: list[str] = field(default_factory=list)


@dataclass
class AnalysisResult:
    """Complete analysis result."""
    
    source_file: str
    snp_count: int
    variants_analyzed: int
    variants_found: int
    
    apoe_status: Optional[APOEResult]
    variant_results: list[VariantResult]
    
    # Categorized results
    high_priority: list[VariantResult] = field(default_factory=list)
    moderate_priority: list[VariantResult] = field(default_factory=list)
    low_priority: list[VariantResult] = field(default_factory=list)
    
    # By category
    by_category: dict[Category, list[VariantResult]] = field(default_factory=dict)
    
    # Gene interactions detected
    gene_interactions: list[dict] = field(default_factory=list)


class ClinicalAnalyzer:
    """
    Analyze DNA data for clinically significant variants.
    
    Example:
        >>> from pydna_analyzer import load_dna_data, ClinicalAnalyzer
        >>> data = load_dna_data("AncestryDNA.txt")
        >>> analyzer = ClinicalAnalyzer()
        >>> result = analyzer.analyze(data)
        >>> print(f"Found {result.variants_found} clinical variants")
    """
    
    def __init__(self, include_low_evidence: bool = True):
        """
        Initialize the analyzer.
        
        Args:
            include_low_evidence: Include variants with limited evidence
        """
        self.include_low_evidence = include_low_evidence
    
    def analyze(self, dataset: DNADataset) -> AnalysisResult:
        """
        Analyze a DNA dataset for clinical variants.
        
        Args:
            dataset: DNADataset to analyze
        
        Returns:
            AnalysisResult with all findings
        """
        df = dataset.dataframe
        variant_results = []
        rsids_in_data = set(df['rsid'].values)
        
        # Analyze each clinical variant
        for rsid, variant in CLINICAL_VARIANTS.items():
            if rsid not in rsids_in_data:
                continue
            
            # Skip limited evidence if requested
            if not self.include_low_evidence and variant.evidence == EvidenceLevel.LIMITED:
                continue
            
            row = df[df['rsid'] == rsid].iloc[0]
            genotype = row['genotype']
            
            # Find matching genotype interpretation
            interpretation = self._match_genotype(genotype, variant)
            
            result = VariantResult(
                rsid=rsid,
                gene=variant.gene,
                name=variant.name,
                category=variant.category,
                genotype=genotype,
                matched_genotype=interpretation[0] if interpretation else None,
                risk_level=interpretation[1].risk if interpretation else RiskLevel.UNKNOWN,
                risk_allele=variant.risk_allele,
                description=interpretation[1].description if interpretation else "Genotype not in database",
                effect=variant.effect,
                evidence=variant.evidence,
                interventions=variant.interventions,
                chromosome=str(row['chromosome']),
                position=int(row['position']),
                has_interactions=variant.interacts_with is not None,
            )
            variant_results.append(result)
        
        # Determine APOE status
        apoe_status = determine_apoe_from_dataset(dataset)
        
        # Categorize results
        high_priority = []
        moderate_priority = []
        low_priority = []
        by_category: dict[Category, list[VariantResult]] = {}
        
        for r in variant_results:
            # By risk level
            if r.risk_level in [RiskLevel.HIGH, RiskLevel.ELEVATED]:
                high_priority.append(r)
            elif r.risk_level == RiskLevel.MODERATE:
                moderate_priority.append(r)
            else:
                low_priority.append(r)
            
            # By category
            if r.category not in by_category:
                by_category[r.category] = []
            by_category[r.category].append(r)
        
        # Sort by evidence strength
        high_priority.sort(key=lambda x: (x.evidence != EvidenceLevel.STRONG, x.category.value))
        moderate_priority.sort(key=lambda x: (x.evidence != EvidenceLevel.STRONG, x.category.value))
        
        # Detect gene interactions
        gene_interactions = self._detect_interactions(variant_results)
        
        # Add interaction notes to affected variants
        for interaction in gene_interactions:
            for r in variant_results:
                if r.rsid in interaction['rsids']:
                    r.interaction_notes.append(interaction['note'])
        
        return AnalysisResult(
            source_file=dataset.source_file,
            snp_count=dataset.snp_count,
            variants_analyzed=len(CLINICAL_VARIANTS),
            variants_found=len(variant_results),
            apoe_status=apoe_status,
            variant_results=variant_results,
            high_priority=high_priority,
            moderate_priority=moderate_priority,
            low_priority=low_priority,
            by_category=by_category,
            gene_interactions=gene_interactions,
        )
    
    def _match_genotype(self, genotype: str, variant: ClinicalVariant):
        """Match a genotype to the variant's known genotypes."""
        for candidate in self._genotype_candidates(genotype):
            if candidate in variant.genotypes:
                return (candidate, variant.genotypes[candidate])

        return None

    def _genotype_candidates(self, genotype: str) -> list[str]:
        """Generate normalized genotype candidates, including reverse complements."""
        candidates: list[str] = []
        seen: set[str] = set()

        def add(candidate: str) -> None:
            if candidate and candidate not in seen:
                seen.add(candidate)
                candidates.append(candidate)

        add(genotype)
        add(genotype[::-1])
        add("".join(sorted(genotype)))

        complemented = self._reverse_complement_genotype(genotype)
        if complemented is not None:
            add(complemented)
            add(complemented[::-1])
            add("".join(sorted(complemented)))

        return candidates

    def _reverse_complement_genotype(self, genotype: str) -> Optional[str]:
        """Return the reverse-complemented genotype for standard nucleotide calls."""
        complement_map = {"A": "T", "T": "A", "C": "G", "G": "C"}

        try:
            return "".join(complement_map[allele] for allele in genotype)
        except KeyError:
            return None
    
    def _detect_interactions(self, results: list[VariantResult]) -> list[dict]:
        """Detect known gene-gene interactions in results."""
        interactions = []
        analyzed_pairs = set()
        
        for r in results:
            if not r.has_interactions:
                continue
            
            variant = CLINICAL_VARIANTS.get(r.rsid)
            if not variant or not variant.interacts_with:
                continue
            
            for interacting_rsid in variant.interacts_with:
                # Skip if already analyzed this pair
                pair = tuple(sorted([r.rsid, interacting_rsid]))
                if pair in analyzed_pairs:
                    continue
                analyzed_pairs.add(pair)
                
                # Check if interacting variant is in results
                interacting_result = next(
                    (x for x in results if x.rsid == interacting_rsid), None
                )
                if not interacting_result:
                    continue
                
                # Check for compound effects
                interaction = self._analyze_interaction(r, interacting_result, variant)
                if interaction:
                    interactions.append(interaction)
        
        return interactions
    
    def _analyze_interaction(
        self,
        result1: VariantResult,
        result2: VariantResult,
        variant1: ClinicalVariant,
    ) -> Optional[dict]:
        """Analyze a specific gene-gene interaction."""
        # MTHFR compound heterozygosity
        if result1.gene == "MTHFR" and result2.gene == "MTHFR":
            if (result1.risk_level in [RiskLevel.MODERATE, RiskLevel.ELEVATED] and
                result2.risk_level in [RiskLevel.MODERATE, RiskLevel.ELEVATED, RiskLevel.LOW]):
                return {
                    "type": "compound_heterozygosity",
                    "genes": ["MTHFR"],
                    "rsids": [result1.rsid, result2.rsid],
                    "severity": "elevated",
                    "note": (
                        f"Compound heterozygosity detected: {result1.name} ({result1.genotype}) + "
                        f"{result2.name} ({result2.genotype}). This combination significantly reduces "
                        "MTHFR enzyme activity. Methylfolate supplementation strongly recommended."
                    ),
                    "recommendations": [
                        "Methylfolate (L-methylfolate) supplementation",
                        "Avoid folic acid supplements (use methylfolate instead)",
                        "Monitor homocysteine levels",
                        "B12 (methylcobalamin) co-supplementation",
                    ],
                }
        
        # APOE interaction already handled separately
        if result1.gene == "APOE" and result2.gene == "APOE":
            return None
        
        # Generic interaction note
        if variant1.interaction_effect:
            return {
                "type": "gene_interaction",
                "genes": [result1.gene, result2.gene],
                "rsids": [result1.rsid, result2.rsid],
                "severity": "moderate",
                "note": variant1.interaction_effect,
                "recommendations": list(set(result1.interventions + result2.interventions)),
            }
        
        return None
