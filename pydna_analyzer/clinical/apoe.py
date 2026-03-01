"""
APOE genotype determination from rs429358 and rs7412.

APOE alleles:
- e2: rs429358=T, rs7412=T (protective for Alzheimer's)
- e3: rs429358=T, rs7412=C (most common, neutral)
- e4: rs429358=C, rs7412=C (risk for Alzheimer's and cardiovascular disease)
"""

from __future__ import annotations

from enum import Enum
from typing import Optional

from pydantic import BaseModel


class APOEAllele(str, Enum):
    """APOE allele types."""
    E2 = "e2"
    E3 = "e3"
    E4 = "e4"
    UNKNOWN = "unknown"


class APOEResult(BaseModel):
    """APOE genotype determination result."""
    
    genotype: str  # e.g., "e3/e4"
    allele1: APOEAllele
    allele2: APOEAllele
    rs429358_genotype: Optional[str]
    rs7412_genotype: Optional[str]
    risk_category: str
    interpretation: str
    recommendations: list[str]
    
    @property
    def has_e4(self) -> bool:
        """Check if genotype contains e4 allele."""
        return self.allele1 == APOEAllele.E4 or self.allele2 == APOEAllele.E4
    
    @property
    def is_e4_homozygous(self) -> bool:
        """Check if homozygous for e4."""
        return self.allele1 == APOEAllele.E4 and self.allele2 == APOEAllele.E4
    
    @property
    def has_e2(self) -> bool:
        """Check if genotype contains e2 allele."""
        return self.allele1 == APOEAllele.E2 or self.allele2 == APOEAllele.E2


# APOE risk interpretations
APOE_INTERPRETATIONS = {
    "e2/e2": {
        "risk": "reduced",
        "interpretation": "Protective genotype for Alzheimer's disease. Slightly increased risk for type III hyperlipoproteinemia.",
        "recommendations": [
            "Monitor lipid levels annually",
            "Standard cardiovascular care",
            "Generally protective cognitive profile",
        ],
    },
    "e2/e3": {
        "risk": "reduced",
        "interpretation": "Below average risk for Alzheimer's disease due to e2 allele.",
        "recommendations": [
            "Standard preventive care",
            "Maintain healthy lifestyle",
        ],
    },
    "e2/e4": {
        "risk": "average",
        "interpretation": "Mixed genotype - e2 may partially offset e4 risk. Overall near-average risk.",
        "recommendations": [
            "Regular cardiovascular monitoring",
            "Brain-healthy lifestyle recommended",
            "Discuss with healthcare provider",
        ],
    },
    "e3/e3": {
        "risk": "average",
        "interpretation": "Most common genotype (60-70% of population). Average baseline risk for Alzheimer's and cardiovascular disease.",
        "recommendations": [
            "Standard preventive care",
            "Healthy diet and exercise",
        ],
    },
    "e3/e4": {
        "risk": "elevated",
        "interpretation": "One e4 allele increases Alzheimer's risk ~2-3x compared to e3/e3. Also elevated cardiovascular risk.",
        "recommendations": [
            "Regular cognitive screening after age 50",
            "Mediterranean/MIND diet recommended",
            "Regular cardiovascular monitoring",
            "Consider omega-3 supplementation",
            "Prioritize quality sleep and exercise",
        ],
    },
    "e4/e4": {
        "risk": "high",
        "interpretation": "Highest genetic risk for late-onset Alzheimer's (~8-12x increased). Significant cardiovascular implications.",
        "recommendations": [
            "Consider genetic counseling",
            "Proactive cognitive health measures",
            "Mediterranean/MIND diet strongly recommended",
            "Regular cognitive and cardiovascular screening",
            "Discuss risk reduction strategies with specialist",
            "Exercise and sleep optimization critical",
        ],
    },
}


def _determine_allele(rs429358_base: str, rs7412_base: str) -> APOEAllele:
    """Determine APOE allele from individual bases at each position."""
    if rs429358_base == 'T' and rs7412_base == 'T':
        return APOEAllele.E2
    elif rs429358_base == 'T' and rs7412_base == 'C':
        return APOEAllele.E3
    elif rs429358_base == 'C' and rs7412_base == 'C':
        return APOEAllele.E4
    else:
        return APOEAllele.UNKNOWN


def determine_apoe_status(
    rs429358_genotype: Optional[str],
    rs7412_genotype: Optional[str],
) -> Optional[APOEResult]:
    """
    Determine APOE genotype from rs429358 and rs7412 genotypes.
    
    Args:
        rs429358_genotype: Genotype at rs429358 (e.g., "TC", "TT", "CC")
        rs7412_genotype: Genotype at rs7412 (e.g., "CT", "CC", "TT")
    
    Returns:
        APOEResult with full interpretation, or None if insufficient data
    """
    if not rs429358_genotype or not rs7412_genotype:
        return None
    
    if len(rs429358_genotype) < 2 or len(rs7412_genotype) < 2:
        return None
    
    # Determine each allele
    allele1 = _determine_allele(rs429358_genotype[0], rs7412_genotype[0])
    allele2 = _determine_allele(rs429358_genotype[1], rs7412_genotype[1])
    
    # Sort alleles for consistent genotype string (e2 < e3 < e4)
    allele_order = {APOEAllele.E2: 0, APOEAllele.E3: 1, APOEAllele.E4: 2, APOEAllele.UNKNOWN: 3}
    if allele_order[allele1] > allele_order[allele2]:
        allele1, allele2 = allele2, allele1
    
    genotype = f"{allele1.value}/{allele2.value}"
    
    # Get interpretation
    interp = APOE_INTERPRETATIONS.get(genotype, {
        "risk": "unknown",
        "interpretation": "Unable to determine risk profile.",
        "recommendations": ["Consult genetic counselor"],
    })
    
    return APOEResult(
        genotype=genotype,
        allele1=allele1,
        allele2=allele2,
        rs429358_genotype=rs429358_genotype,
        rs7412_genotype=rs7412_genotype,
        risk_category=interp["risk"],
        interpretation=interp["interpretation"],
        recommendations=interp["recommendations"],
    )


def determine_apoe_from_dataset(dataset) -> Optional[APOEResult]:
    """
    Determine APOE status from a DNADataset.
    
    Args:
        dataset: DNADataset containing genotype data
    
    Returns:
        APOEResult or None if rs429358/rs7412 not found
    """
    rs429358 = dataset.get_genotype("rs429358")
    rs7412 = dataset.get_genotype("rs7412")
    
    return determine_apoe_status(rs429358, rs7412)
