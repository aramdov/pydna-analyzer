"""
Curated database of clinically significant SNPs.

Categories:
- cardiovascular: Heart and circulatory system
- metabolic: Diabetes and metabolic syndrome
- cancer: Cancer risk variants
- pharmacogenomics: Drug metabolism
- nutrient: Nutrient metabolism
- neuro: Neurological and mental health
"""

from __future__ import annotations

from enum import Enum
from typing import Optional

from pydantic import BaseModel


class RiskLevel(str, Enum):
    """Risk level classifications."""
    NORMAL = "normal"
    LOW = "low"
    MODERATE = "moderate"
    ELEVATED = "elevated"
    HIGH = "high"
    UNKNOWN = "unknown"


class EvidenceLevel(str, Enum):
    """Evidence strength for clinical associations."""
    STRONG = "strong"
    MODERATE = "moderate"
    LIMITED = "limited"


class Category(str, Enum):
    """Clinical variant categories."""
    CARDIOVASCULAR = "cardiovascular"
    METABOLIC = "metabolic"
    CANCER = "cancer"
    PHARMACOGENOMICS = "pharmacogenomics"
    NUTRIENT = "nutrient"
    NEURO = "neuro"


class GenotypeInfo(BaseModel):
    """Information about a specific genotype."""
    risk: RiskLevel
    description: str


class ClinicalVariant(BaseModel):
    """A clinically significant genetic variant."""
    rsid: str
    gene: str
    name: str
    category: Category
    risk_allele: str
    effect: str
    genotypes: dict[str, GenotypeInfo]
    evidence: EvidenceLevel
    interventions: list[str]
    
    # Optional gene interaction info
    interacts_with: Optional[list[str]] = None
    interaction_effect: Optional[str] = None


# ============================================================================
# CLINICAL VARIANTS DATABASE
# ============================================================================

CLINICAL_VARIANTS: dict[str, ClinicalVariant] = {
    # =========================================================================
    # CARDIOVASCULAR
    # =========================================================================
    "rs1801133": ClinicalVariant(
        rsid="rs1801133",
        gene="MTHFR",
        name="MTHFR C677T",
        category=Category.CARDIOVASCULAR,
        risk_allele="T",
        effect="Reduced enzyme activity affecting folate metabolism and homocysteine levels",
        genotypes={
            "CC": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal MTHFR activity"),
            "CT": GenotypeInfo(risk=RiskLevel.MODERATE, description="~35% reduced enzyme activity"),
            "TT": GenotypeInfo(risk=RiskLevel.ELEVATED, description="~70% reduced enzyme activity, elevated homocysteine risk"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["Methylfolate supplementation", "B12 monitoring", "Homocysteine testing"],
        interacts_with=["rs1801131"],
        interaction_effect="Compound heterozygosity with A1298C significantly reduces MTHFR function",
    ),
    "rs1801131": ClinicalVariant(
        rsid="rs1801131",
        gene="MTHFR",
        name="MTHFR A1298C",
        category=Category.CARDIOVASCULAR,
        risk_allele="C",
        effect="Reduced enzyme activity, less severe than C677T",
        genotypes={
            "AA": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal activity"),
            "AC": GenotypeInfo(risk=RiskLevel.LOW, description="Slightly reduced activity"),
            "CC": GenotypeInfo(risk=RiskLevel.MODERATE, description="Reduced activity, compound effect with C677T"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["Consider methylfolate if compound heterozygous with C677T"],
        interacts_with=["rs1801133"],
        interaction_effect="Compound heterozygosity with C677T significantly reduces MTHFR function",
    ),
    "rs1799983": ClinicalVariant(
        rsid="rs1799983",
        gene="NOS3",
        name="eNOS G894T",
        category=Category.CARDIOVASCULAR,
        risk_allele="T",
        effect="Reduced nitric oxide production, affecting blood vessel dilation",
        genotypes={
            "GG": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal NO production"),
            "GT": GenotypeInfo(risk=RiskLevel.MODERATE, description="Reduced NO production"),
            "TT": GenotypeInfo(risk=RiskLevel.ELEVATED, description="Significantly reduced NO, hypertension risk"),
        },
        evidence=EvidenceLevel.MODERATE,
        interventions=["L-arginine supplementation", "Regular BP monitoring", "Nitrate-rich foods"],
    ),
    "rs429358": ClinicalVariant(
        rsid="rs429358",
        gene="APOE",
        name="APOE e4 determinant 1",
        category=Category.CARDIOVASCULAR,
        risk_allele="C",
        effect="Part of APOE e4 haplotype - cardiovascular and Alzheimer's risk",
        genotypes={
            "TT": GenotypeInfo(risk=RiskLevel.NORMAL, description="Not e4 carrier at this position"),
            "TC": GenotypeInfo(risk=RiskLevel.ELEVATED, description="One e4 allele possible"),
            "CT": GenotypeInfo(risk=RiskLevel.ELEVATED, description="One e4 allele possible"),
            "CC": GenotypeInfo(risk=RiskLevel.HIGH, description="e4/e4 possible - check rs7412"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["Lipid monitoring", "Mediterranean diet", "Cognitive screening after 50"],
        interacts_with=["rs7412"],
        interaction_effect="Combined with rs7412 determines APOE isoform (e2/e3/e4)",
    ),
    "rs7412": ClinicalVariant(
        rsid="rs7412",
        gene="APOE",
        name="APOE e4 determinant 2",
        category=Category.CARDIOVASCULAR,
        risk_allele="C",
        effect="Second determinant of APOE isoform",
        genotypes={
            "CC": GenotypeInfo(risk=RiskLevel.NORMAL, description="e3 or e4 (check rs429358)"),
            "CT": GenotypeInfo(risk=RiskLevel.LOW, description="One e2 allele - may be protective"),
            "TC": GenotypeInfo(risk=RiskLevel.LOW, description="One e2 allele - may be protective"),
            "TT": GenotypeInfo(risk=RiskLevel.LOW, description="e2/e2 - protective for Alzheimer's, hyperlipoproteinemia risk"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["Lipid panel monitoring", "Varies by full APOE genotype"],
        interacts_with=["rs429358"],
        interaction_effect="Combined with rs429358 determines APOE isoform (e2/e3/e4)",
    ),
    "rs6025": ClinicalVariant(
        rsid="rs6025",
        gene="F5",
        name="Factor V Leiden",
        category=Category.CARDIOVASCULAR,
        risk_allele="A",
        effect="Increased blood clotting risk (thrombophilia)",
        genotypes={
            "GG": GenotypeInfo(risk=RiskLevel.NORMAL, description="No Factor V Leiden"),
            "GA": GenotypeInfo(risk=RiskLevel.ELEVATED, description="Heterozygous - 5-10x clot risk"),
            "AG": GenotypeInfo(risk=RiskLevel.ELEVATED, description="Heterozygous - 5-10x clot risk"),
            "AA": GenotypeInfo(risk=RiskLevel.HIGH, description="Homozygous - 50-100x clot risk"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["Avoid prolonged immobility", "Discuss anticoagulation before surgery", "Avoid estrogen-containing contraceptives"],
    ),
    "rs1799963": ClinicalVariant(
        rsid="rs1799963",
        gene="F2",
        name="Prothrombin G20210A",
        category=Category.CARDIOVASCULAR,
        risk_allele="A",
        effect="Increased prothrombin levels and clotting risk",
        genotypes={
            "GG": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal prothrombin"),
            "GA": GenotypeInfo(risk=RiskLevel.ELEVATED, description="2-3x increased clot risk"),
            "AG": GenotypeInfo(risk=RiskLevel.ELEVATED, description="2-3x increased clot risk"),
            "AA": GenotypeInfo(risk=RiskLevel.HIGH, description="Significantly elevated clot risk"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["Same as Factor V Leiden", "Genetic counseling recommended"],
    ),

    # =========================================================================
    # METABOLIC / DIABETES
    # =========================================================================
    "rs7903146": ClinicalVariant(
        rsid="rs7903146",
        gene="TCF7L2",
        name="TCF7L2 diabetes risk",
        category=Category.METABOLIC,
        risk_allele="T",
        effect="Strongest common genetic risk factor for type 2 diabetes",
        genotypes={
            "CC": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal diabetes risk"),
            "CT": GenotypeInfo(risk=RiskLevel.MODERATE, description="~1.4x increased T2D risk"),
            "TC": GenotypeInfo(risk=RiskLevel.MODERATE, description="~1.4x increased T2D risk"),
            "TT": GenotypeInfo(risk=RiskLevel.ELEVATED, description="~2x increased T2D risk"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["Regular glucose monitoring", "Low glycemic diet", "Exercise", "Maintain healthy weight"],
    ),
    "rs1801282": ClinicalVariant(
        rsid="rs1801282",
        gene="PPARG",
        name="PPARG Pro12Ala",
        category=Category.METABOLIC,
        risk_allele="C",
        effect="Affects insulin sensitivity and fat metabolism",
        genotypes={
            "CC": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal insulin sensitivity"),
            "CG": GenotypeInfo(risk=RiskLevel.LOW, description="Slightly improved insulin sensitivity"),
            "GC": GenotypeInfo(risk=RiskLevel.LOW, description="Slightly improved insulin sensitivity"),
            "GG": GenotypeInfo(risk=RiskLevel.LOW, description="Better insulin sensitivity, lower T2D risk"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["Lifestyle optimization still beneficial"],
    ),
    "rs5219": ClinicalVariant(
        rsid="rs5219",
        gene="KCNJ11",
        name="KCNJ11 E23K",
        category=Category.METABOLIC,
        risk_allele="T",
        effect="Affects insulin secretion from pancreatic beta cells",
        genotypes={
            "CC": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal insulin secretion"),
            "CT": GenotypeInfo(risk=RiskLevel.LOW, description="Slightly reduced insulin secretion"),
            "TC": GenotypeInfo(risk=RiskLevel.LOW, description="Slightly reduced insulin secretion"),
            "TT": GenotypeInfo(risk=RiskLevel.MODERATE, description="Reduced insulin secretion, T2D risk"),
        },
        evidence=EvidenceLevel.MODERATE,
        interventions=["Monitor fasting glucose", "Sulfonylureas may be effective"],
    ),
    "rs1801260": ClinicalVariant(
        rsid="rs1801260",
        gene="CLOCK",
        name="CLOCK circadian gene",
        category=Category.METABOLIC,
        risk_allele="C",
        effect="Affects circadian rhythm and metabolic regulation",
        genotypes={
            "TT": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal circadian function"),
            "TC": GenotypeInfo(risk=RiskLevel.LOW, description="Mild circadian disruption tendency"),
            "CT": GenotypeInfo(risk=RiskLevel.LOW, description="Mild circadian disruption tendency"),
            "CC": GenotypeInfo(risk=RiskLevel.MODERATE, description="Evening chronotype, metabolic syndrome risk"),
        },
        evidence=EvidenceLevel.MODERATE,
        interventions=["Consistent sleep schedule", "Time-restricted eating", "Morning light exposure"],
    ),

    # =========================================================================
    # CANCER RISK
    # =========================================================================
    "rs1042522": ClinicalVariant(
        rsid="rs1042522",
        gene="TP53",
        name="p53 Arg72Pro",
        category=Category.CANCER,
        risk_allele="C",
        effect="Affects p53 tumor suppressor function",
        genotypes={
            "GG": GenotypeInfo(risk=RiskLevel.NORMAL, description="Arginine variant - better apoptosis"),
            "GC": GenotypeInfo(risk=RiskLevel.LOW, description="Heterozygous"),
            "CG": GenotypeInfo(risk=RiskLevel.LOW, description="Heterozygous"),
            "CC": GenotypeInfo(risk=RiskLevel.MODERATE, description="Proline variant - altered cancer risk profile"),
        },
        evidence=EvidenceLevel.MODERATE,
        interventions=["Standard cancer screening", "Avoid carcinogens", "Antioxidant-rich diet"],
    ),
    "rs1800566": ClinicalVariant(
        rsid="rs1800566",
        gene="NQO1",
        name="NQO1 Pro187Ser",
        category=Category.CANCER,
        risk_allele="T",
        effect="Affects detoxification of quinones and carcinogens",
        genotypes={
            "CC": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal NQO1 activity"),
            "CT": GenotypeInfo(risk=RiskLevel.MODERATE, description="Reduced detoxification"),
            "TC": GenotypeInfo(risk=RiskLevel.MODERATE, description="Reduced detoxification"),
            "TT": GenotypeInfo(risk=RiskLevel.ELEVATED, description="No NQO1 activity - benzene sensitivity, leukemia risk"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["Avoid benzene exposure", "Avoid smoking", "Cruciferous vegetables"],
    ),
    "rs1048943": ClinicalVariant(
        rsid="rs1048943",
        gene="CYP1A1",
        name="CYP1A1 Ile462Val",
        category=Category.CANCER,
        risk_allele="G",
        effect="Affects metabolism of polycyclic aromatic hydrocarbons",
        genotypes={
            "AA": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal PAH metabolism"),
            "AG": GenotypeInfo(risk=RiskLevel.MODERATE, description="Increased carcinogen activation"),
            "GA": GenotypeInfo(risk=RiskLevel.MODERATE, description="Increased carcinogen activation"),
            "GG": GenotypeInfo(risk=RiskLevel.ELEVATED, description="Higher lung cancer risk if smoking"),
        },
        evidence=EvidenceLevel.MODERATE,
        interventions=["Never smoke", "Avoid air pollution", "Cruciferous vegetables"],
    ),
    "rs1695": ClinicalVariant(
        rsid="rs1695",
        gene="GSTP1",
        name="GSTP1 Ile105Val",
        category=Category.CANCER,
        risk_allele="G",
        effect="Affects glutathione-mediated detoxification",
        genotypes={
            "AA": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal GSTP1 activity"),
            "AG": GenotypeInfo(risk=RiskLevel.LOW, description="Slightly reduced detox capacity"),
            "GA": GenotypeInfo(risk=RiskLevel.LOW, description="Slightly reduced detox capacity"),
            "GG": GenotypeInfo(risk=RiskLevel.MODERATE, description="Reduced detoxification, cancer risk"),
        },
        evidence=EvidenceLevel.MODERATE,
        interventions=["Sulforaphane (broccoli sprouts)", "NAC supplementation", "Reduce toxin exposure"],
    ),

    # =========================================================================
    # PHARMACOGENOMICS
    # =========================================================================
    "rs1799853": ClinicalVariant(
        rsid="rs1799853",
        gene="CYP2C9",
        name="CYP2C9*2",
        category=Category.PHARMACOGENOMICS,
        risk_allele="T",
        effect="Reduced metabolism of warfarin, NSAIDs, sulfonylureas",
        genotypes={
            "CC": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal CYP2C9 activity (*1/*1)"),
            "CT": GenotypeInfo(risk=RiskLevel.MODERATE, description="Intermediate metabolizer (*1/*2)"),
            "TC": GenotypeInfo(risk=RiskLevel.MODERATE, description="Intermediate metabolizer (*1/*2)"),
            "TT": GenotypeInfo(risk=RiskLevel.ELEVATED, description="Poor metabolizer (*2/*2) - dose reduction needed"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["Lower warfarin doses", "NSAID caution", "Pharmacist consultation"],
    ),
    "rs1057910": ClinicalVariant(
        rsid="rs1057910",
        gene="CYP2C9",
        name="CYP2C9*3",
        category=Category.PHARMACOGENOMICS,
        risk_allele="C",
        effect="More severely reduced CYP2C9 activity than *2",
        genotypes={
            "AA": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal activity (*1/*1)"),
            "AC": GenotypeInfo(risk=RiskLevel.ELEVATED, description="Intermediate metabolizer (*1/*3)"),
            "CA": GenotypeInfo(risk=RiskLevel.ELEVATED, description="Intermediate metabolizer (*1/*3)"),
            "CC": GenotypeInfo(risk=RiskLevel.HIGH, description="Poor metabolizer (*3/*3) - significant dose reduction"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["50-75% warfarin dose reduction", "Genetic dosing algorithms", "INR monitoring"],
    ),
    "rs4244285": ClinicalVariant(
        rsid="rs4244285",
        gene="CYP2C19",
        name="CYP2C19*2",
        category=Category.PHARMACOGENOMICS,
        risk_allele="A",
        effect="Loss of function - affects clopidogrel (Plavix), PPIs, antidepressants",
        genotypes={
            "GG": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal metabolizer"),
            "GA": GenotypeInfo(risk=RiskLevel.MODERATE, description="Intermediate metabolizer - reduced clopidogrel effect"),
            "AG": GenotypeInfo(risk=RiskLevel.MODERATE, description="Intermediate metabolizer - reduced clopidogrel effect"),
            "AA": GenotypeInfo(risk=RiskLevel.ELEVATED, description="Poor metabolizer - clopidogrel ineffective"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["Alternative to clopidogrel (prasugrel, ticagrelor)", "PPI dose adjustment"],
    ),
    "rs4986893": ClinicalVariant(
        rsid="rs4986893",
        gene="CYP2C19",
        name="CYP2C19*3",
        category=Category.PHARMACOGENOMICS,
        risk_allele="A",
        effect="Loss of function allele (more common in Asian populations)",
        genotypes={
            "GG": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal metabolizer"),
            "GA": GenotypeInfo(risk=RiskLevel.MODERATE, description="Intermediate metabolizer"),
            "AG": GenotypeInfo(risk=RiskLevel.MODERATE, description="Intermediate metabolizer"),
            "AA": GenotypeInfo(risk=RiskLevel.ELEVATED, description="Poor metabolizer"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["Same as CYP2C19*2"],
    ),
    "rs1065852": ClinicalVariant(
        rsid="rs1065852",
        gene="CYP2D6",
        name="CYP2D6*10",
        category=Category.PHARMACOGENOMICS,
        risk_allele="A",
        effect="Reduced activity - affects codeine, tamoxifen, many antidepressants",
        genotypes={
            "GG": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal CYP2D6"),
            "GA": GenotypeInfo(risk=RiskLevel.MODERATE, description="Intermediate metabolizer"),
            "AG": GenotypeInfo(risk=RiskLevel.MODERATE, description="Intermediate metabolizer"),
            "AA": GenotypeInfo(risk=RiskLevel.ELEVATED, description="Poor metabolizer - codeine ineffective, tamoxifen reduced"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["Alternative pain meds to codeine", "Consider aromatase inhibitors over tamoxifen"],
    ),
    "rs3892097": ClinicalVariant(
        rsid="rs3892097",
        gene="CYP2D6",
        name="CYP2D6*4",
        category=Category.PHARMACOGENOMICS,
        risk_allele="A",
        effect="Non-functional allele - most common cause of poor metabolizer status",
        genotypes={
            "GG": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal metabolizer"),
            "GA": GenotypeInfo(risk=RiskLevel.MODERATE, description="Intermediate metabolizer"),
            "AG": GenotypeInfo(risk=RiskLevel.MODERATE, description="Intermediate metabolizer"),
            "AA": GenotypeInfo(risk=RiskLevel.ELEVATED, description="Poor metabolizer - many drugs affected"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["Avoid codeine/tramadol", "Antidepressant dose adjustments", "Pharmacogenomic card"],
    ),
    "rs9923231": ClinicalVariant(
        rsid="rs9923231",
        gene="VKORC1",
        name="VKORC1 -1639G>A",
        category=Category.PHARMACOGENOMICS,
        risk_allele="A",
        effect="Affects warfarin sensitivity - key for dosing",
        genotypes={
            "GG": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal warfarin dose needed"),
            "GA": GenotypeInfo(risk=RiskLevel.MODERATE, description="Intermediate sensitivity - lower dose"),
            "AG": GenotypeInfo(risk=RiskLevel.MODERATE, description="Intermediate sensitivity - lower dose"),
            "AA": GenotypeInfo(risk=RiskLevel.ELEVATED, description="High sensitivity - much lower dose needed"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["Genetic-guided warfarin dosing", "More frequent INR monitoring initially"],
    ),
    "rs1800460": ClinicalVariant(
        rsid="rs1800460",
        gene="TPMT",
        name="TPMT*3B",
        category=Category.PHARMACOGENOMICS,
        risk_allele="A",
        effect="Reduced thiopurine metabolism - azathioprine, 6-MP toxicity risk",
        genotypes={
            "GG": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal TPMT activity"),
            "GA": GenotypeInfo(risk=RiskLevel.ELEVATED, description="Intermediate - 50% dose reduction"),
            "AG": GenotypeInfo(risk=RiskLevel.ELEVATED, description="Intermediate - 50% dose reduction"),
            "AA": GenotypeInfo(risk=RiskLevel.HIGH, description="Poor metabolizer - severe toxicity risk, avoid thiopurines"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["TPMT testing before thiopurine therapy", "Dose reduction or alternative drugs"],
    ),
    "rs1142345": ClinicalVariant(
        rsid="rs1142345",
        gene="TPMT",
        name="TPMT*3C",
        category=Category.PHARMACOGENOMICS,
        risk_allele="C",
        effect="Reduced TPMT activity",
        genotypes={
            "AA": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal TPMT activity"),
            "AC": GenotypeInfo(risk=RiskLevel.ELEVATED, description="Intermediate metabolizer"),
            "CA": GenotypeInfo(risk=RiskLevel.ELEVATED, description="Intermediate metabolizer"),
            "CC": GenotypeInfo(risk=RiskLevel.HIGH, description="Poor metabolizer - thiopurine toxicity risk"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["Same as TPMT*3B"],
    ),
    "rs4680": ClinicalVariant(
        rsid="rs4680",
        gene="COMT",
        name="COMT Val158Met",
        category=Category.PHARMACOGENOMICS,
        risk_allele="A",
        effect="Affects dopamine/catecholamine metabolism - pain sensitivity, stress response",
        genotypes={
            "GG": GenotypeInfo(risk=RiskLevel.NORMAL, description="Val/Val - faster dopamine clearance, stress resilient"),
            "GA": GenotypeInfo(risk=RiskLevel.LOW, description="Val/Met - intermediate"),
            "AG": GenotypeInfo(risk=RiskLevel.LOW, description="Val/Met - intermediate"),
            "AA": GenotypeInfo(risk=RiskLevel.MODERATE, description="Met/Met - slower clearance, higher pain sensitivity, anxiety risk"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["Stress management", "May need lower opioid doses", "Magnesium supplementation"],
    ),

    # =========================================================================
    # NUTRIENT METABOLISM
    # =========================================================================
    "rs4988235": ClinicalVariant(
        rsid="rs4988235",
        gene="MCM6/LCT",
        name="Lactase persistence",
        category=Category.NUTRIENT,
        risk_allele="G",
        effect="Determines ability to digest lactose in adulthood",
        genotypes={
            "AA": GenotypeInfo(risk=RiskLevel.ELEVATED, description="Lactose intolerant - lactase non-persistent"),
            "AG": GenotypeInfo(risk=RiskLevel.NORMAL, description="Lactose tolerant - one persistence allele"),
            "GA": GenotypeInfo(risk=RiskLevel.NORMAL, description="Lactose tolerant - one persistence allele"),
            "GG": GenotypeInfo(risk=RiskLevel.NORMAL, description="Lactose tolerant - lactase persistent"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["Lactase supplements if intolerant", "Calcium from non-dairy sources"],
    ),
    "rs602662": ClinicalVariant(
        rsid="rs602662",
        gene="FUT2",
        name="FUT2 secretor status",
        category=Category.NUTRIENT,
        risk_allele="A",
        effect="Affects B12 absorption and gut microbiome",
        genotypes={
            "GG": GenotypeInfo(risk=RiskLevel.NORMAL, description="Secretor - normal B12 absorption"),
            "GA": GenotypeInfo(risk=RiskLevel.NORMAL, description="Secretor"),
            "AG": GenotypeInfo(risk=RiskLevel.NORMAL, description="Secretor"),
            "AA": GenotypeInfo(risk=RiskLevel.MODERATE, description="Non-secretor - reduced B12 absorption, different microbiome"),
        },
        evidence=EvidenceLevel.MODERATE,
        interventions=["Monitor B12 levels", "Consider methylcobalamin supplementation"],
    ),
    "rs1799945": ClinicalVariant(
        rsid="rs1799945",
        gene="HFE",
        name="HFE H63D",
        category=Category.NUTRIENT,
        risk_allele="G",
        effect="Affects iron absorption - hereditary hemochromatosis risk",
        genotypes={
            "CC": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal iron absorption"),
            "CG": GenotypeInfo(risk=RiskLevel.LOW, description="Carrier - mild increased absorption"),
            "GC": GenotypeInfo(risk=RiskLevel.LOW, description="Carrier - mild increased absorption"),
            "GG": GenotypeInfo(risk=RiskLevel.MODERATE, description="Homozygous - monitor iron levels"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["Regular ferritin monitoring", "Avoid iron supplements unless deficient", "Blood donation"],
    ),
    "rs1800562": ClinicalVariant(
        rsid="rs1800562",
        gene="HFE",
        name="HFE C282Y",
        category=Category.NUTRIENT,
        risk_allele="A",
        effect="Major hereditary hemochromatosis mutation",
        genotypes={
            "GG": GenotypeInfo(risk=RiskLevel.NORMAL, description="No C282Y mutation"),
            "GA": GenotypeInfo(risk=RiskLevel.MODERATE, description="Carrier - monitor iron, especially with H63D"),
            "AG": GenotypeInfo(risk=RiskLevel.MODERATE, description="Carrier - monitor iron, especially with H63D"),
            "AA": GenotypeInfo(risk=RiskLevel.HIGH, description="Homozygous - hemochromatosis risk, regular phlebotomy may be needed"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["Ferritin and transferrin saturation monitoring", "Therapeutic phlebotomy if elevated", "Avoid vitamin C with meals"],
    ),
    "rs12934922": ClinicalVariant(
        rsid="rs12934922",
        gene="BCMO1",
        name="BCMO1 beta-carotene conversion",
        category=Category.NUTRIENT,
        risk_allele="T",
        effect="Reduced conversion of beta-carotene to vitamin A",
        genotypes={
            "AA": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal conversion"),
            "AT": GenotypeInfo(risk=RiskLevel.MODERATE, description="~30% reduced conversion"),
            "TA": GenotypeInfo(risk=RiskLevel.MODERATE, description="~30% reduced conversion"),
            "TT": GenotypeInfo(risk=RiskLevel.ELEVATED, description="~60% reduced conversion - need preformed vitamin A"),
        },
        evidence=EvidenceLevel.MODERATE,
        interventions=["Include preformed vitamin A (retinol) from animal sources", "Cod liver oil", "Egg yolks"],
    ),
    "rs7946": ClinicalVariant(
        rsid="rs7946",
        gene="PEMT",
        name="PEMT choline synthesis",
        category=Category.NUTRIENT,
        risk_allele="T",
        effect="Reduced endogenous choline production",
        genotypes={
            "CC": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal choline synthesis"),
            "CT": GenotypeInfo(risk=RiskLevel.LOW, description="Slightly reduced synthesis"),
            "TC": GenotypeInfo(risk=RiskLevel.LOW, description="Slightly reduced synthesis"),
            "TT": GenotypeInfo(risk=RiskLevel.MODERATE, description="Reduced synthesis - higher dietary choline needs"),
        },
        evidence=EvidenceLevel.MODERATE,
        interventions=["Eggs (high choline)", "Liver", "Choline supplementation"],
    ),
    "rs2282679": ClinicalVariant(
        rsid="rs2282679",
        gene="GC",
        name="Vitamin D binding protein",
        category=Category.NUTRIENT,
        risk_allele="C",
        effect="Affects vitamin D transport and availability",
        genotypes={
            "AA": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal vitamin D binding"),
            "AC": GenotypeInfo(risk=RiskLevel.MODERATE, description="Lower circulating 25(OH)D"),
            "CA": GenotypeInfo(risk=RiskLevel.MODERATE, description="Lower circulating 25(OH)D"),
            "CC": GenotypeInfo(risk=RiskLevel.ELEVATED, description="Lower vitamin D levels - higher supplementation needs"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["Regular vitamin D testing", "Higher supplementation doses", "Sun exposure"],
    ),
    "rs12785878": ClinicalVariant(
        rsid="rs12785878",
        gene="DHCR7/NADSYN1",
        name="Vitamin D synthesis",
        category=Category.NUTRIENT,
        risk_allele="G",
        effect="Affects vitamin D production in skin",
        genotypes={
            "TT": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal vitamin D synthesis"),
            "GT": GenotypeInfo(risk=RiskLevel.LOW, description="Slightly reduced synthesis"),
            "TG": GenotypeInfo(risk=RiskLevel.LOW, description="Slightly reduced synthesis"),
            "GG": GenotypeInfo(risk=RiskLevel.MODERATE, description="Reduced synthesis - supplementation important"),
        },
        evidence=EvidenceLevel.MODERATE,
        interventions=["Vitamin D supplementation", "More sun exposure needed"],
    ),

    # =========================================================================
    # NEUROLOGICAL
    # =========================================================================
    "rs6265": ClinicalVariant(
        rsid="rs6265",
        gene="BDNF",
        name="BDNF Val66Met",
        category=Category.NEURO,
        risk_allele="T",
        effect="Affects brain-derived neurotrophic factor secretion - memory, mood",
        genotypes={
            "CC": GenotypeInfo(risk=RiskLevel.NORMAL, description="Val/Val - normal BDNF secretion"),
            "CT": GenotypeInfo(risk=RiskLevel.LOW, description="Val/Met - reduced activity-dependent secretion"),
            "TC": GenotypeInfo(risk=RiskLevel.LOW, description="Val/Met - reduced activity-dependent secretion"),
            "TT": GenotypeInfo(risk=RiskLevel.MODERATE, description="Met/Met - reduced BDNF, depression/anxiety risk, memory effects"),
        },
        evidence=EvidenceLevel.STRONG,
        interventions=["Regular exercise (increases BDNF)", "Omega-3 fatty acids", "Stress reduction", "Sleep optimization"],
    ),
    "rs1800497": ClinicalVariant(
        rsid="rs1800497",
        gene="ANKK1/DRD2",
        name="DRD2 Taq1A",
        category=Category.NEURO,
        risk_allele="T",
        effect="Affects dopamine receptor density - reward sensitivity, addiction risk",
        genotypes={
            "CC": GenotypeInfo(risk=RiskLevel.NORMAL, description="A2/A2 - normal dopamine receptors"),
            "CT": GenotypeInfo(risk=RiskLevel.MODERATE, description="A1/A2 - reduced receptor density"),
            "TC": GenotypeInfo(risk=RiskLevel.MODERATE, description="A1/A2 - reduced receptor density"),
            "TT": GenotypeInfo(risk=RiskLevel.ELEVATED, description="A1/A1 - ~30% fewer D2 receptors, addiction vulnerability"),
        },
        evidence=EvidenceLevel.MODERATE,
        interventions=["Avoid addictive substances", "Behavioral addiction awareness", "Dopamine-supporting supplements"],
    ),
    "rs6311": ClinicalVariant(
        rsid="rs6311",
        gene="HTR2A",
        name="Serotonin receptor 2A",
        category=Category.NEURO,
        risk_allele="T",
        effect="Affects serotonin signaling - antidepressant response, mood",
        genotypes={
            "CC": GenotypeInfo(risk=RiskLevel.NORMAL, description="Normal serotonin receptor expression"),
            "CT": GenotypeInfo(risk=RiskLevel.LOW, description="Altered expression"),
            "TC": GenotypeInfo(risk=RiskLevel.LOW, description="Altered expression"),
            "TT": GenotypeInfo(risk=RiskLevel.MODERATE, description="Different antidepressant response profile"),
        },
        evidence=EvidenceLevel.MODERATE,
        interventions=["May respond differently to SSRIs", "Tryptophan-rich foods"],
    ),
    "rs25531": ClinicalVariant(
        rsid="rs25531",
        gene="SLC6A4",
        name="Serotonin transporter (5-HTTLPR related)",
        category=Category.NEURO,
        risk_allele="G",
        effect="Affects serotonin reuptake - stress response, depression risk",
        genotypes={
            "AA": GenotypeInfo(risk=RiskLevel.NORMAL, description="Long/Long - resilient to stress"),
            "AG": GenotypeInfo(risk=RiskLevel.MODERATE, description="Long/Short - intermediate stress sensitivity"),
            "GA": GenotypeInfo(risk=RiskLevel.MODERATE, description="Long/Short - intermediate stress sensitivity"),
            "GG": GenotypeInfo(risk=RiskLevel.ELEVATED, description="Short/Short - higher stress sensitivity, depression risk"),
        },
        evidence=EvidenceLevel.MODERATE,
        interventions=["Stress management essential", "Strong social support", "Mindfulness practices"],
    ),
}


def get_variant(rsid: str) -> Optional[ClinicalVariant]:
    """Get a clinical variant by rsID."""
    return CLINICAL_VARIANTS.get(rsid)


def get_variants_by_category(category: Category) -> list[ClinicalVariant]:
    """Get all variants in a category."""
    return [v for v in CLINICAL_VARIANTS.values() if v.category == category]


def get_all_rsids() -> set[str]:
    """Get set of all rsIDs in the database."""
    return set(CLINICAL_VARIANTS.keys())


def get_interacting_variants(rsid: str) -> list[ClinicalVariant]:
    """Get variants that interact with the given rsID."""
    variant = CLINICAL_VARIANTS.get(rsid)
    if variant and variant.interacts_with:
        return [CLINICAL_VARIANTS[r] for r in variant.interacts_with if r in CLINICAL_VARIANTS]
    return []
