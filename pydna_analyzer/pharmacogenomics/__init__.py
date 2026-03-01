"""
Pharmacogenomics module for star allele calling and drug response prediction.

Pipeline:
    DNADataset -> extract per-gene SNPs -> call star alleles -> form diplotype
               -> compute activity score -> predict phenotype -> lookup drug recs -> PGxResult

Key design decisions:
    - No pypgx dependency: self-contained, works with consumer DNA data
    - *1 by absence: wild-type allele assigned when no variant alleles match
    - CPIC numeric activity scores (0, 0.25, 0.5, 1.0)
    - Unphased data: assumes variant alleles on different chromosomes (most likely)
    - Single-SNP genes (VKORC1, COMT): direct genotype map, no star allele system
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from pydna_analyzer.core.data_loader import DNADataset

# =============================================================================
# Enums
# =============================================================================


class MetabolizerPhenotype(str, Enum):
    """Predicted metabolizer phenotype based on CPIC guidelines."""

    POOR = "Poor Metabolizer"
    INTERMEDIATE = "Intermediate Metabolizer"
    NORMAL = "Normal Metabolizer"
    RAPID = "Rapid Metabolizer"
    ULTRARAPID = "Ultrarapid Metabolizer"
    # VKORC1 sensitivity phenotypes
    HIGH_SENSITIVITY = "High Sensitivity"
    INTERMEDIATE_SENSITIVITY = "Intermediate Sensitivity"
    NORMAL_SENSITIVITY = "Normal Sensitivity"
    # COMT activity phenotypes
    HIGH_ACTIVITY = "High Activity"
    INTERMEDIATE_ACTIVITY = "Intermediate Activity"
    LOW_ACTIVITY = "Low Activity"


class CPICLevel(str, Enum):
    """CPIC recommendation strength level."""

    STRONG = "strong"
    MODERATE = "moderate"
    OPTIONAL = "optional"


# =============================================================================
# Data Definitions (frozen dataclasses)
# =============================================================================


@dataclass(frozen=True)
class StarAlleleDefinition:
    """Definition of a star allele by its defining SNPs."""

    name: str
    defining_snps: dict[str, str]  # rsid -> variant allele
    activity_score: float
    function: str  # e.g. "no function", "decreased function", "increased function"


@dataclass(frozen=True)
class GeneDefinition:
    """Definition of a pharmacogene with its star alleles and phenotype thresholds."""

    gene: str
    relevant_rsids: list[str]
    star_alleles: list[StarAlleleDefinition]
    is_single_snp_gene: bool
    phenotype_thresholds: dict[str, tuple[float, float]]  # phenotype -> (min_score, max_score)


@dataclass(frozen=True)
class DrugRecommendation:
    """CPIC-based drug recommendation for a gene-phenotype pair."""

    drug_name: str
    drug_class: str
    gene: str
    phenotype: MetabolizerPhenotype
    recommendation: str
    cpic_level: CPICLevel
    alternatives: list[str]


# =============================================================================
# Gene Database (6 genes)
# =============================================================================

GENE_DEFINITIONS: dict[str, GeneDefinition] = {
    "CYP2C9": GeneDefinition(
        gene="CYP2C9",
        relevant_rsids=["rs1799853", "rs1057910"],
        star_alleles=[
            StarAlleleDefinition(
                name="*2",
                defining_snps={"rs1799853": "T"},
                activity_score=0.5,
                function="decreased function",
            ),
            StarAlleleDefinition(
                name="*3",
                defining_snps={"rs1057910": "C"},
                activity_score=0.0,
                function="no function",
            ),
        ],
        is_single_snp_gene=False,
        phenotype_thresholds={
            MetabolizerPhenotype.POOR.value: (0.0, 0.0),
            MetabolizerPhenotype.INTERMEDIATE.value: (0.25, 1.0),
            MetabolizerPhenotype.NORMAL.value: (1.5, 2.0),
        },
    ),
    "CYP2C19": GeneDefinition(
        gene="CYP2C19",
        relevant_rsids=["rs4244285", "rs4986893", "rs12248560"],
        star_alleles=[
            StarAlleleDefinition(
                name="*2",
                defining_snps={"rs4244285": "A"},
                activity_score=0.0,
                function="no function",
            ),
            StarAlleleDefinition(
                name="*3",
                defining_snps={"rs4986893": "A"},
                activity_score=0.0,
                function="no function",
            ),
            StarAlleleDefinition(
                name="*17",
                defining_snps={"rs12248560": "T"},
                activity_score=1.5,
                function="increased function",
            ),
        ],
        is_single_snp_gene=False,
        phenotype_thresholds={
            MetabolizerPhenotype.POOR.value: (0.0, 0.0),
            MetabolizerPhenotype.INTERMEDIATE.value: (0.5, 1.0),
            MetabolizerPhenotype.NORMAL.value: (1.5, 2.0),
            MetabolizerPhenotype.RAPID.value: (2.5, 2.5),
            MetabolizerPhenotype.ULTRARAPID.value: (3.0, 3.0),
        },
    ),
    "CYP2D6": GeneDefinition(
        gene="CYP2D6",
        relevant_rsids=["rs3892097", "rs1065852"],
        star_alleles=[
            StarAlleleDefinition(
                name="*4",
                defining_snps={"rs3892097": "A"},
                activity_score=0.0,
                function="no function",
            ),
            StarAlleleDefinition(
                name="*10",
                defining_snps={"rs1065852": "A"},
                activity_score=0.25,
                function="decreased function",
            ),
        ],
        is_single_snp_gene=False,
        phenotype_thresholds={
            MetabolizerPhenotype.POOR.value: (0.0, 0.0),
            MetabolizerPhenotype.INTERMEDIATE.value: (0.25, 1.0),
            MetabolizerPhenotype.NORMAL.value: (1.5, 2.0),
        },
    ),
    "TPMT": GeneDefinition(
        gene="TPMT",
        relevant_rsids=["rs1800460", "rs1142345"],
        star_alleles=[
            StarAlleleDefinition(
                name="*3B",
                defining_snps={"rs1800460": "A"},
                activity_score=0.0,
                function="no function",
            ),
            StarAlleleDefinition(
                name="*3C",
                defining_snps={"rs1142345": "C"},
                activity_score=0.0,
                function="no function",
            ),
        ],
        is_single_snp_gene=False,
        phenotype_thresholds={
            MetabolizerPhenotype.POOR.value: (0.0, 0.0),
            MetabolizerPhenotype.INTERMEDIATE.value: (0.5, 1.0),
            MetabolizerPhenotype.NORMAL.value: (1.5, 2.0),
        },
    ),
    "VKORC1": GeneDefinition(
        gene="VKORC1",
        relevant_rsids=["rs9923231"],
        star_alleles=[],
        is_single_snp_gene=True,
        phenotype_thresholds={},
    ),
    "COMT": GeneDefinition(
        gene="COMT",
        relevant_rsids=["rs4680"],
        star_alleles=[],
        is_single_snp_gene=True,
        phenotype_thresholds={},
    ),
}

# Single-SNP genotype-to-phenotype maps
VKORC1_PHENOTYPE_MAP: dict[str, MetabolizerPhenotype] = {
    "GG": MetabolizerPhenotype.NORMAL_SENSITIVITY,
    "GA": MetabolizerPhenotype.INTERMEDIATE_SENSITIVITY,
    "AG": MetabolizerPhenotype.INTERMEDIATE_SENSITIVITY,
    "AA": MetabolizerPhenotype.HIGH_SENSITIVITY,
}

COMT_PHENOTYPE_MAP: dict[str, MetabolizerPhenotype] = {
    "GG": MetabolizerPhenotype.HIGH_ACTIVITY,
    "GA": MetabolizerPhenotype.INTERMEDIATE_ACTIVITY,
    "AG": MetabolizerPhenotype.INTERMEDIATE_ACTIVITY,
    "AA": MetabolizerPhenotype.LOW_ACTIVITY,
}


# =============================================================================
# Drug Recommendations Database
# =============================================================================

DRUG_RECOMMENDATIONS: list[DrugRecommendation] = [
    # CYP2C19 — clopidogrel
    DrugRecommendation(
        drug_name="Clopidogrel",
        drug_class="Antiplatelet",
        gene="CYP2C19",
        phenotype=MetabolizerPhenotype.POOR,
        recommendation="Avoid clopidogrel. Use alternative antiplatelet therapy.",
        cpic_level=CPICLevel.STRONG,
        alternatives=["prasugrel", "ticagrelor"],
    ),
    DrugRecommendation(
        drug_name="Clopidogrel",
        drug_class="Antiplatelet",
        gene="CYP2C19",
        phenotype=MetabolizerPhenotype.INTERMEDIATE,
        recommendation="Reduced clopidogrel activation. Consider alternative antiplatelet.",
        cpic_level=CPICLevel.MODERATE,
        alternatives=["prasugrel", "ticagrelor"],
    ),
    DrugRecommendation(
        drug_name="Clopidogrel",
        drug_class="Antiplatelet",
        gene="CYP2C19",
        phenotype=MetabolizerPhenotype.ULTRARAPID,
        recommendation="Standard clopidogrel dosing expected to be effective.",
        cpic_level=CPICLevel.STRONG,
        alternatives=[],
    ),
    # CYP2C19 — PPIs
    DrugRecommendation(
        drug_name="Omeprazole",
        drug_class="Proton Pump Inhibitor",
        gene="CYP2C19",
        phenotype=MetabolizerPhenotype.POOR,
        recommendation="Decreased omeprazole clearance. Consider 50% dose reduction.",
        cpic_level=CPICLevel.MODERATE,
        alternatives=["rabeprazole (less CYP2C19 dependent)"],
    ),
    DrugRecommendation(
        drug_name="Omeprazole",
        drug_class="Proton Pump Inhibitor",
        gene="CYP2C19",
        phenotype=MetabolizerPhenotype.ULTRARAPID,
        recommendation="Increased omeprazole clearance. May need higher dose or alternative.",
        cpic_level=CPICLevel.OPTIONAL,
        alternatives=["rabeprazole", "pantoprazole"],
    ),
    # CYP2C9 — warfarin
    DrugRecommendation(
        drug_name="Warfarin",
        drug_class="Anticoagulant",
        gene="CYP2C9",
        phenotype=MetabolizerPhenotype.POOR,
        recommendation="Significantly reduced warfarin clearance. Use pharmacogenetic dosing algorithm. "
        "Expect 50-75% dose reduction from standard.",
        cpic_level=CPICLevel.STRONG,
        alternatives=["direct oral anticoagulants (DOACs)"],
    ),
    DrugRecommendation(
        drug_name="Warfarin",
        drug_class="Anticoagulant",
        gene="CYP2C9",
        phenotype=MetabolizerPhenotype.INTERMEDIATE,
        recommendation="Reduced warfarin clearance. Use pharmacogenetic dosing algorithm. "
        "Expect 20-40% dose reduction.",
        cpic_level=CPICLevel.STRONG,
        alternatives=["direct oral anticoagulants (DOACs)"],
    ),
    # CYP2C9 — NSAIDs
    DrugRecommendation(
        drug_name="Celecoxib",
        drug_class="NSAID",
        gene="CYP2C9",
        phenotype=MetabolizerPhenotype.POOR,
        recommendation="Reduced celecoxib metabolism. Use lowest effective dose or alternative.",
        cpic_level=CPICLevel.MODERATE,
        alternatives=["ibuprofen (partially CYP2C9)", "acetaminophen"],
    ),
    DrugRecommendation(
        drug_name="Celecoxib",
        drug_class="NSAID",
        gene="CYP2C9",
        phenotype=MetabolizerPhenotype.INTERMEDIATE,
        recommendation="Moderately reduced celecoxib metabolism. Use with caution.",
        cpic_level=CPICLevel.OPTIONAL,
        alternatives=["acetaminophen"],
    ),
    # CYP2D6 — codeine
    DrugRecommendation(
        drug_name="Codeine",
        drug_class="Opioid Analgesic",
        gene="CYP2D6",
        phenotype=MetabolizerPhenotype.POOR,
        recommendation="Codeine is ineffective (cannot convert to morphine). Avoid codeine and tramadol.",
        cpic_level=CPICLevel.STRONG,
        alternatives=["non-opioid analgesics", "morphine (with caution)"],
    ),
    DrugRecommendation(
        drug_name="Codeine",
        drug_class="Opioid Analgesic",
        gene="CYP2D6",
        phenotype=MetabolizerPhenotype.INTERMEDIATE,
        recommendation="Reduced codeine activation to morphine. May have decreased efficacy.",
        cpic_level=CPICLevel.MODERATE,
        alternatives=["non-opioid analgesics"],
    ),
    # CYP2D6 — tamoxifen
    DrugRecommendation(
        drug_name="Tamoxifen",
        drug_class="Selective Estrogen Receptor Modulator",
        gene="CYP2D6",
        phenotype=MetabolizerPhenotype.POOR,
        recommendation="Reduced tamoxifen activation to endoxifen. Consider alternative therapy.",
        cpic_level=CPICLevel.STRONG,
        alternatives=["aromatase inhibitors (postmenopausal)"],
    ),
    DrugRecommendation(
        drug_name="Tamoxifen",
        drug_class="Selective Estrogen Receptor Modulator",
        gene="CYP2D6",
        phenotype=MetabolizerPhenotype.INTERMEDIATE,
        recommendation="Moderately reduced endoxifen levels. Consider dose increase or alternative.",
        cpic_level=CPICLevel.MODERATE,
        alternatives=["aromatase inhibitors (postmenopausal)"],
    ),
    # TPMT — azathioprine / thiopurines
    DrugRecommendation(
        drug_name="Azathioprine",
        drug_class="Thiopurine Immunosuppressant",
        gene="TPMT",
        phenotype=MetabolizerPhenotype.POOR,
        recommendation="Extremely high risk of life-threatening myelosuppression. "
        "Avoid thiopurines or reduce dose by 90% with close monitoring.",
        cpic_level=CPICLevel.STRONG,
        alternatives=["mycophenolate mofetil", "alternative immunosuppressants"],
    ),
    DrugRecommendation(
        drug_name="Azathioprine",
        drug_class="Thiopurine Immunosuppressant",
        gene="TPMT",
        phenotype=MetabolizerPhenotype.INTERMEDIATE,
        recommendation="Increased risk of myelosuppression. Reduce starting dose by 30-50%.",
        cpic_level=CPICLevel.STRONG,
        alternatives=["dose-adjusted thiopurines with monitoring"],
    ),
    # VKORC1 — warfarin sensitivity
    DrugRecommendation(
        drug_name="Warfarin",
        drug_class="Anticoagulant",
        gene="VKORC1",
        phenotype=MetabolizerPhenotype.HIGH_SENSITIVITY,
        recommendation="High warfarin sensitivity. Requires significantly lower dose. "
        "Use pharmacogenetic dosing algorithm combining CYP2C9 and VKORC1.",
        cpic_level=CPICLevel.STRONG,
        alternatives=["direct oral anticoagulants (DOACs)"],
    ),
    DrugRecommendation(
        drug_name="Warfarin",
        drug_class="Anticoagulant",
        gene="VKORC1",
        phenotype=MetabolizerPhenotype.INTERMEDIATE_SENSITIVITY,
        recommendation="Intermediate warfarin sensitivity. Likely needs lower-than-average dose.",
        cpic_level=CPICLevel.STRONG,
        alternatives=["direct oral anticoagulants (DOACs)"],
    ),
    # COMT — opioid response
    DrugRecommendation(
        drug_name="Opioid Analgesics",
        drug_class="Pain Management",
        gene="COMT",
        phenotype=MetabolizerPhenotype.LOW_ACTIVITY,
        recommendation="Higher pain sensitivity and increased opioid requirements. "
        "May need dose adjustment. Monitor for efficacy.",
        cpic_level=CPICLevel.OPTIONAL,
        alternatives=["multimodal pain management", "non-opioid strategies"],
    ),
]


# =============================================================================
# Results
# =============================================================================


@dataclass
class GeneResult:
    """Result of pharmacogenomic analysis for a single gene."""

    gene: str
    diplotype: str
    phenotype: MetabolizerPhenotype
    activity_score: float | None
    snps_tested: int
    snps_missing: int
    drug_recommendations: list[DrugRecommendation]
    interpretation: str
    confidence: str  # "high", "moderate", "low"

    @property
    def is_actionable(self) -> bool:
        """Whether this result has clinically actionable drug recommendations."""
        return len(self.drug_recommendations) > 0 and self.phenotype not in (
            MetabolizerPhenotype.NORMAL,
            MetabolizerPhenotype.NORMAL_SENSITIVITY,
            MetabolizerPhenotype.HIGH_ACTIVITY,
        )


@dataclass
class PGxResult:
    """Complete pharmacogenomics analysis result."""

    gene_results: list[GeneResult]
    total_genes_tested: int
    total_snps_tested: int
    total_snps_missing: int

    @property
    def all_drug_recommendations(self) -> list[DrugRecommendation]:
        """All drug recommendations across all genes."""
        recs: list[DrugRecommendation] = []
        for result in self.gene_results:
            recs.extend(result.drug_recommendations)
        return recs

    @property
    def actionable_count(self) -> int:
        """Number of genes with actionable (non-normal) results."""
        return sum(1 for r in self.gene_results if r.is_actionable)

    @property
    def summary(self) -> str:
        """Human-readable summary of results."""
        actionable = self.actionable_count
        total = self.total_genes_tested
        if actionable == 0:
            return (
                f"Analyzed {total} pharmacogenes. "
                "No actionable drug-gene interactions detected."
            )
        return (
            f"Analyzed {total} pharmacogenes. "
            f"Found {actionable} gene(s) with actionable drug recommendations."
        )


# =============================================================================
# PGxAnalyzer
# =============================================================================


class PGxAnalyzer:
    """Pharmacogenomics analyzer for star allele calling and drug recommendations."""

    def __init__(self) -> None:
        self.gene_definitions = GENE_DEFINITIONS
        self.drug_recommendations = DRUG_RECOMMENDATIONS

    def analyze(self, dataset: DNADataset) -> PGxResult:
        """Run pharmacogenomics analysis on a DNA dataset.

        Args:
            dataset: Loaded DNA dataset with genotype data.

        Returns:
            PGxResult with per-gene diplotypes, phenotypes, and drug recommendations.
        """
        gene_results: list[GeneResult] = []
        total_snps_tested = 0
        total_snps_missing = 0

        for _gene_name, gene_def in self.gene_definitions.items():
            result = self._analyze_gene(dataset, gene_def)
            gene_results.append(result)
            total_snps_tested += result.snps_tested
            total_snps_missing += result.snps_missing

        return PGxResult(
            gene_results=gene_results,
            total_genes_tested=len(gene_results),
            total_snps_tested=total_snps_tested,
            total_snps_missing=total_snps_missing,
        )

    def analyze_gene(self, dataset: DNADataset, gene: str) -> GeneResult | None:
        """Analyze a single gene.

        Args:
            dataset: Loaded DNA dataset.
            gene: Gene name (e.g. "CYP2C19").

        Returns:
            GeneResult for the specified gene, or None if gene not found.
        """
        gene_def = self.gene_definitions.get(gene)
        if gene_def is None:
            return None
        return self._analyze_gene(dataset, gene_def)

    def _analyze_gene(self, dataset: DNADataset, gene_def: GeneDefinition) -> GeneResult:
        """Analyze a single gene definition against the dataset."""
        # Extract genotypes for this gene's SNPs
        genotypes: dict[str, str | None] = {}
        snps_tested = 0
        snps_missing = 0

        for rsid in gene_def.relevant_rsids:
            gt = dataset.get_genotype(rsid)
            genotypes[rsid] = gt
            if gt is not None:
                snps_tested += 1
            else:
                snps_missing += 1

        # Route to appropriate analysis path
        if gene_def.is_single_snp_gene:
            return self._analyze_single_snp_gene(gene_def, genotypes, snps_tested, snps_missing)
        return self._analyze_multi_snp_gene(gene_def, genotypes, snps_tested, snps_missing)

    def _analyze_single_snp_gene(
        self,
        gene_def: GeneDefinition,
        genotypes: dict[str, str | None],
        snps_tested: int,
        snps_missing: int,
    ) -> GeneResult:
        """Analyze a single-SNP gene (VKORC1, COMT) using direct genotype mapping."""
        rsid = gene_def.relevant_rsids[0]
        gt = genotypes.get(rsid)

        if gene_def.gene == "VKORC1":
            phenotype_map = VKORC1_PHENOTYPE_MAP
            default_phenotype = MetabolizerPhenotype.NORMAL_SENSITIVITY
        else:  # COMT
            phenotype_map = COMT_PHENOTYPE_MAP
            default_phenotype = MetabolizerPhenotype.HIGH_ACTIVITY

        if gt is None:
            phenotype = default_phenotype
            diplotype = "Unknown"
            confidence = "low"
            interpretation = (
                f"{gene_def.gene}: SNP {rsid} not found in dataset. "
                "Assuming normal/wild-type."
            )
        else:
            phenotype = phenotype_map.get(gt, default_phenotype)
            diplotype = gt
            confidence = "high"
            interpretation = f"{gene_def.gene} genotype {gt}: {phenotype.value}."

        # Lookup drug recommendations
        drug_recs = self._get_drug_recommendations(gene_def.gene, phenotype)

        return GeneResult(
            gene=gene_def.gene,
            diplotype=diplotype,
            phenotype=phenotype,
            activity_score=None,
            snps_tested=snps_tested,
            snps_missing=snps_missing,
            drug_recommendations=drug_recs,
            interpretation=interpretation,
            confidence=confidence,
        )

    def _analyze_multi_snp_gene(
        self,
        gene_def: GeneDefinition,
        genotypes: dict[str, str | None],
        snps_tested: int,
        snps_missing: int,
    ) -> GeneResult:
        """Analyze a multi-SNP gene using star allele calling."""
        # Call diplotype
        allele1, allele2 = self._call_diplotype(gene_def, genotypes)

        # Compute activity score
        score1 = self._allele_activity_score(gene_def, allele1)
        score2 = self._allele_activity_score(gene_def, allele2)
        activity_score = score1 + score2

        # Predict phenotype
        phenotype = self._score_to_phenotype(gene_def, activity_score)

        # Determine confidence
        if snps_missing == 0:
            confidence = "high"
        elif snps_missing < len(gene_def.relevant_rsids):
            confidence = "moderate"
        else:
            confidence = "low"

        diplotype_str = f"{allele1}/{allele2}"

        # Build interpretation
        interpretation = (
            f"{gene_def.gene} {diplotype_str}: {phenotype.value} "
            f"(activity score {activity_score:.1f})."
        )
        if snps_missing > 0:
            interpretation += (
                f" Note: {snps_missing} of {len(gene_def.relevant_rsids)} "
                "SNPs not found in dataset."
            )

        # Lookup drug recommendations
        drug_recs = self._get_drug_recommendations(gene_def.gene, phenotype)

        return GeneResult(
            gene=gene_def.gene,
            diplotype=diplotype_str,
            phenotype=phenotype,
            activity_score=activity_score,
            snps_tested=snps_tested,
            snps_missing=snps_missing,
            drug_recommendations=drug_recs,
            interpretation=interpretation,
            confidence=confidence,
        )

    def _call_diplotype(
        self,
        gene_def: GeneDefinition,
        genotypes: dict[str, str | None],
    ) -> tuple[str, str]:
        """Call diplotype from observed genotypes.

        For unphased consumer data, we count variant alleles per star allele definition
        and assign *1 (wild-type) as default. When multiple variant alleles are detected,
        we assume they are on different chromosomes (most likely for unphased data).

        Returns:
            Tuple of (allele1, allele2) star allele names.
        """
        # Count variant alleles observed for each star allele
        variant_alleles_found: list[str] = []

        for star_allele in gene_def.star_alleles:
            count = self._count_variant_alleles(star_allele, genotypes)
            for _ in range(count):
                variant_alleles_found.append(star_allele.name)

        # Assign diplotype
        if len(variant_alleles_found) == 0:
            return ("*1", "*1")
        elif len(variant_alleles_found) == 1:
            return ("*1", variant_alleles_found[0])
        else:
            # Two or more variant alleles — take first two
            return (variant_alleles_found[0], variant_alleles_found[1])

    def _count_variant_alleles(
        self,
        star_allele: StarAlleleDefinition,
        genotypes: dict[str, str | None],
    ) -> int:
        """Count how many copies of a variant allele are present.

        For star alleles defined by a single SNP, counts the number of variant
        alleles in the genotype (0, 1, or 2).
        """
        count = 0
        for rsid, variant_allele in star_allele.defining_snps.items():
            gt = genotypes.get(rsid)
            if gt is None:
                continue
            # Count occurrences of variant allele in genotype
            allele1 = gt[0] if len(gt) >= 1 else ""
            allele2 = gt[1] if len(gt) >= 2 else ""
            if allele1 == variant_allele:
                count += 1
            if allele2 == variant_allele:
                count += 1
        return count

    def _allele_activity_score(self, gene_def: GeneDefinition, allele_name: str) -> float:
        """Get the activity score for a star allele. *1 (wild-type) = 1.0."""
        if allele_name == "*1":
            return 1.0
        for star_allele in gene_def.star_alleles:
            if star_allele.name == allele_name:
                return star_allele.activity_score
        return 1.0  # Default to wild-type if unknown

    def _score_to_phenotype(
        self, gene_def: GeneDefinition, activity_score: float
    ) -> MetabolizerPhenotype:
        """Map activity score to phenotype using gene-specific thresholds."""
        for phenotype_value, (min_score, max_score) in gene_def.phenotype_thresholds.items():
            if min_score <= activity_score <= max_score:
                return MetabolizerPhenotype(phenotype_value)
        # Fallback: if score doesn't match any range, use nearest
        return MetabolizerPhenotype.NORMAL

    def _get_drug_recommendations(
        self, gene: str, phenotype: MetabolizerPhenotype
    ) -> list[DrugRecommendation]:
        """Look up drug recommendations for a gene-phenotype pair."""
        return [
            rec
            for rec in self.drug_recommendations
            if rec.gene == gene and rec.phenotype == phenotype
        ]
