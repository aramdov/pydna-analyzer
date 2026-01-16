#!/usr/bin/env python3
"""
AncestryDNA Clinical Variant Analyzer
Extracts and annotates clinically significant SNPs from raw DNA data.
"""

import pandas as pd
import requests
import json
import time
from typing import Dict, List, Optional, Tuple

# Curated database of clinically significant SNPs with strong evidence
# Categories: cardiovascular, metabolic, cancer, neuro, pharmacogenomics, nutrient
CLINICAL_SNPS = {
    # === CARDIOVASCULAR ===
    "rs1801133": {
        "gene": "MTHFR",
        "category": "cardiovascular",
        "name": "MTHFR C677T",
        "risk_allele": "T",
        "effect": "Reduced enzyme activity affecting folate metabolism and homocysteine levels",
        "genotypes": {
            "CC": {"risk": "normal", "description": "Normal MTHFR activity"},
            "CT": {"risk": "moderate", "description": "~35% reduced enzyme activity"},
            "TT": {"risk": "elevated", "description": "~70% reduced enzyme activity, elevated homocysteine risk"}
        },
        "evidence": "strong",
        "interventions": ["Methylfolate supplementation", "B12 monitoring", "Homocysteine testing"]
    },
    "rs1801131": {
        "gene": "MTHFR",
        "category": "cardiovascular",
        "name": "MTHFR A1298C",
        "risk_allele": "C",
        "effect": "Reduced enzyme activity, less severe than C677T",
        "genotypes": {
            "AA": {"risk": "normal", "description": "Normal activity"},
            "AC": {"risk": "low", "description": "Slightly reduced activity"},
            "CC": {"risk": "moderate", "description": "Reduced activity, compound effect with C677T"}
        },
        "evidence": "strong",
        "interventions": ["Consider methylfolate if compound heterozygous with C677T"]
    },
    "rs1799983": {
        "gene": "NOS3",
        "category": "cardiovascular",
        "name": "eNOS G894T",
        "risk_allele": "T",
        "effect": "Reduced nitric oxide production, affecting blood vessel dilation",
        "genotypes": {
            "GG": {"risk": "normal", "description": "Normal NO production"},
            "GT": {"risk": "moderate", "description": "Reduced NO production"},
            "TT": {"risk": "elevated", "description": "Significantly reduced NO, hypertension risk"}
        },
        "evidence": "moderate",
        "interventions": ["L-arginine supplementation", "Regular BP monitoring", "Nitrate-rich foods"]
    },
    "rs4340": {
        "gene": "ACE",
        "category": "cardiovascular",
        "name": "ACE I/D",
        "risk_allele": "D",
        "effect": "Insertion/deletion affecting ACE levels and cardiovascular risk",
        "genotypes": {
            "II": {"risk": "normal", "description": "Lower ACE activity"},
            "ID": {"risk": "moderate", "description": "Intermediate ACE activity"},
            "DD": {"risk": "elevated", "description": "Higher ACE activity, hypertension risk"}
        },
        "evidence": "strong",
        "interventions": ["ACE inhibitors may be more effective", "Monitor blood pressure"]
    },
    "rs1799752": {
        "gene": "ACE",
        "category": "cardiovascular",
        "name": "ACE I/D (alternate marker)",
        "risk_allele": "D",
        "effect": "ACE insertion/deletion polymorphism",
        "genotypes": {
            "II": {"risk": "normal", "description": "Lower ACE levels"},
            "ID": {"risk": "moderate", "description": "Intermediate ACE levels"},
            "DD": {"risk": "elevated", "description": "Higher ACE levels, cardiovascular risk"}
        },
        "evidence": "strong",
        "interventions": ["BP monitoring", "Consider ACE inhibitor response"]
    },
    "rs429358": {
        "gene": "APOE",
        "category": "cardiovascular",
        "name": "APOE e4 determinant 1",
        "risk_allele": "C",
        "effect": "Part of APOE e4 haplotype - cardiovascular and Alzheimer's risk",
        "genotypes": {
            "TT": {"risk": "normal", "description": "Not e4 carrier at this position"},
            "TC": {"risk": "elevated", "description": "One e4 allele possible"},
            "CC": {"risk": "high", "description": "e4/e4 possible - check rs7412"}
        },
        "evidence": "strong",
        "interventions": ["Lipid monitoring", "Mediterranean diet", "Cognitive screening after 50"]
    },
    "rs7412": {
        "gene": "APOE",
        "category": "cardiovascular",
        "name": "APOE e4 determinant 2",
        "risk_allele": "C",
        "effect": "Second determinant of APOE isoform",
        "genotypes": {
            "CC": {"risk": "normal", "description": "e3 or e4 (check rs429358)"},
            "CT": {"risk": "low", "description": "One e2 allele - may be protective"},
            "TT": {"risk": "low", "description": "e2/e2 - protective for Alzheimer's, hyperlipoproteinemia risk"}
        },
        "evidence": "strong",
        "interventions": ["Lipid panel monitoring", "Varies by full APOE genotype"]
    },
    "rs6025": {
        "gene": "F5",
        "category": "cardiovascular",
        "name": "Factor V Leiden",
        "risk_allele": "A",
        "effect": "Increased blood clotting risk (thrombophilia)",
        "genotypes": {
            "GG": {"risk": "normal", "description": "No Factor V Leiden"},
            "GA": {"risk": "elevated", "description": "Heterozygous - 5-10x clot risk"},
            "AA": {"risk": "high", "description": "Homozygous - 50-100x clot risk"}
        },
        "evidence": "strong",
        "interventions": ["Avoid prolonged immobility", "Discuss anticoagulation before surgery", "Avoid estrogen-containing contraceptives"]
    },
    "rs1799963": {
        "gene": "F2",
        "category": "cardiovascular",
        "name": "Prothrombin G20210A",
        "risk_allele": "A",
        "effect": "Increased prothrombin levels and clotting risk",
        "genotypes": {
            "GG": {"risk": "normal", "description": "Normal prothrombin"},
            "GA": {"risk": "elevated", "description": "2-3x increased clot risk"},
            "AA": {"risk": "high", "description": "Significantly elevated clot risk"}
        },
        "evidence": "strong",
        "interventions": ["Same as Factor V Leiden", "Genetic counseling recommended"]
    },

    # === METABOLIC / DIABETES ===
    "rs7903146": {
        "gene": "TCF7L2",
        "category": "metabolic",
        "name": "TCF7L2 diabetes risk",
        "risk_allele": "T",
        "effect": "Strongest common genetic risk factor for type 2 diabetes",
        "genotypes": {
            "CC": {"risk": "normal", "description": "Normal diabetes risk"},
            "CT": {"risk": "moderate", "description": "~1.4x increased T2D risk"},
            "TT": {"risk": "elevated", "description": "~2x increased T2D risk"}
        },
        "evidence": "strong",
        "interventions": ["Regular glucose monitoring", "Low glycemic diet", "Exercise", "Maintain healthy weight"]
    },
    "rs1801282": {
        "gene": "PPARG",
        "category": "metabolic",
        "name": "PPARG Pro12Ala",
        "risk_allele": "C",
        "effect": "Affects insulin sensitivity and fat metabolism",
        "genotypes": {
            "CC": {"risk": "normal", "description": "Normal insulin sensitivity"},
            "CG": {"risk": "low", "description": "Slightly improved insulin sensitivity"},
            "GG": {"risk": "low", "description": "Better insulin sensitivity, lower T2D risk"}
        },
        "evidence": "strong",
        "interventions": ["Lifestyle optimization still beneficial"]
    },
    "rs5219": {
        "gene": "KCNJ11",
        "category": "metabolic",
        "name": "KCNJ11 E23K",
        "risk_allele": "T",
        "effect": "Affects insulin secretion from pancreatic beta cells",
        "genotypes": {
            "CC": {"risk": "normal", "description": "Normal insulin secretion"},
            "CT": {"risk": "low", "description": "Slightly reduced insulin secretion"},
            "TT": {"risk": "moderate", "description": "Reduced insulin secretion, T2D risk"}
        },
        "evidence": "moderate",
        "interventions": ["Monitor fasting glucose", "Sulfonylureas may be effective"]
    },
    "rs1801260": {
        "gene": "CLOCK",
        "category": "metabolic",
        "name": "CLOCK circadian gene",
        "risk_allele": "C",
        "effect": "Affects circadian rhythm and metabolic regulation",
        "genotypes": {
            "TT": {"risk": "normal", "description": "Normal circadian function"},
            "TC": {"risk": "low", "description": "Mild circadian disruption tendency"},
            "CC": {"risk": "moderate", "description": "Evening chronotype, metabolic syndrome risk"}
        },
        "evidence": "moderate",
        "interventions": ["Consistent sleep schedule", "Time-restricted eating", "Morning light exposure"]
    },

    # === CANCER RISK ===
    "rs1042522": {
        "gene": "TP53",
        "category": "cancer",
        "name": "p53 Arg72Pro",
        "risk_allele": "C",
        "effect": "Affects p53 tumor suppressor function",
        "genotypes": {
            "GG": {"risk": "normal", "description": "Arginine variant - better apoptosis"},
            "GC": {"risk": "low", "description": "Heterozygous"},
            "CC": {"risk": "moderate", "description": "Proline variant - altered cancer risk profile"}
        },
        "evidence": "moderate",
        "interventions": ["Standard cancer screening", "Avoid carcinogens", "Antioxidant-rich diet"]
    },
    "rs1800566": {
        "gene": "NQO1",
        "category": "cancer",
        "name": "NQO1 Pro187Ser",
        "risk_allele": "T",
        "effect": "Affects detoxification of quinones and carcinogens",
        "genotypes": {
            "CC": {"risk": "normal", "description": "Normal NQO1 activity"},
            "CT": {"risk": "moderate", "description": "Reduced detoxification"},
            "TT": {"risk": "elevated", "description": "No NQO1 activity - benzene sensitivity, leukemia risk"}
        },
        "evidence": "strong",
        "interventions": ["Avoid benzene exposure", "Avoid smoking", "Cruciferous vegetables"]
    },
    "rs1048943": {
        "gene": "CYP1A1",
        "category": "cancer",
        "name": "CYP1A1 Ile462Val",
        "risk_allele": "G",
        "effect": "Affects metabolism of polycyclic aromatic hydrocarbons",
        "genotypes": {
            "AA": {"risk": "normal", "description": "Normal PAH metabolism"},
            "AG": {"risk": "moderate", "description": "Increased carcinogen activation"},
            "GG": {"risk": "elevated", "description": "Higher lung cancer risk if smoking"}
        },
        "evidence": "moderate",
        "interventions": ["Never smoke", "Avoid air pollution", "Cruciferous vegetables"]
    },
    "rs1695": {
        "gene": "GSTP1",
        "category": "cancer",
        "name": "GSTP1 Ile105Val",
        "risk_allele": "G",
        "effect": "Affects glutathione-mediated detoxification",
        "genotypes": {
            "AA": {"risk": "normal", "description": "Normal GSTP1 activity"},
            "AG": {"risk": "low", "description": "Slightly reduced detox capacity"},
            "GG": {"risk": "moderate", "description": "Reduced detoxification, cancer risk"}
        },
        "evidence": "moderate",
        "interventions": ["Sulforaphane (broccoli sprouts)", "NAC supplementation", "Reduce toxin exposure"]
    },
    "rs4986850": {
        "gene": "BRCA1",
        "category": "cancer",
        "name": "BRCA1 variant",
        "risk_allele": "varies",
        "effect": "DNA repair gene - mutations increase breast/ovarian cancer risk",
        "genotypes": {
            "common": {"risk": "normal", "description": "Common variant"},
            "rare": {"risk": "high", "description": "If pathogenic - significantly elevated cancer risk"}
        },
        "evidence": "strong",
        "interventions": ["Genetic counseling if concerning variant", "Enhanced screening"]
    },

    # === PHARMACOGENOMICS ===
    "rs1799853": {
        "gene": "CYP2C9",
        "category": "pharmacogenomics",
        "name": "CYP2C9*2",
        "risk_allele": "T",
        "effect": "Reduced metabolism of warfarin, NSAIDs, sulfonylureas",
        "genotypes": {
            "CC": {"risk": "normal", "description": "Normal CYP2C9 activity (*1/*1)"},
            "CT": {"risk": "moderate", "description": "Intermediate metabolizer (*1/*2)"},
            "TT": {"risk": "elevated", "description": "Poor metabolizer (*2/*2) - dose reduction needed"}
        },
        "evidence": "strong",
        "interventions": ["Lower warfarin doses", "NSAID caution", "Pharmacist consultation"]
    },
    "rs1057910": {
        "gene": "CYP2C9",
        "category": "pharmacogenomics",
        "name": "CYP2C9*3",
        "risk_allele": "C",
        "effect": "More severely reduced CYP2C9 activity than *2",
        "genotypes": {
            "AA": {"risk": "normal", "description": "Normal activity (*1/*1)"},
            "AC": {"risk": "elevated", "description": "Intermediate metabolizer (*1/*3)"},
            "CC": {"risk": "high", "description": "Poor metabolizer (*3/*3) - significant dose reduction"}
        },
        "evidence": "strong",
        "interventions": ["50-75% warfarin dose reduction", "Genetic dosing algorithms", "INR monitoring"]
    },
    "rs4244285": {
        "gene": "CYP2C19",
        "category": "pharmacogenomics",
        "name": "CYP2C19*2",
        "risk_allele": "A",
        "effect": "Loss of function - affects clopidogrel (Plavix), PPIs, antidepressants",
        "genotypes": {
            "GG": {"risk": "normal", "description": "Normal metabolizer"},
            "GA": {"risk": "moderate", "description": "Intermediate metabolizer - reduced clopidogrel effect"},
            "AA": {"risk": "elevated", "description": "Poor metabolizer - clopidogrel ineffective"}
        },
        "evidence": "strong",
        "interventions": ["Alternative to clopidogrel (prasugrel, ticagrelor)", "PPI dose adjustment"]
    },
    "rs4986893": {
        "gene": "CYP2C19",
        "category": "pharmacogenomics",
        "name": "CYP2C19*3",
        "risk_allele": "A",
        "effect": "Loss of function allele (more common in Asian populations)",
        "genotypes": {
            "GG": {"risk": "normal", "description": "Normal metabolizer"},
            "GA": {"risk": "moderate", "description": "Intermediate metabolizer"},
            "AA": {"risk": "elevated", "description": "Poor metabolizer"}
        },
        "evidence": "strong",
        "interventions": ["Same as CYP2C19*2"]
    },
    "rs1065852": {
        "gene": "CYP2D6",
        "category": "pharmacogenomics",
        "name": "CYP2D6*10",
        "risk_allele": "A",
        "effect": "Reduced activity - affects codeine, tamoxifen, many antidepressants",
        "genotypes": {
            "GG": {"risk": "normal", "description": "Normal CYP2D6"},
            "GA": {"risk": "moderate", "description": "Intermediate metabolizer"},
            "AA": {"risk": "elevated", "description": "Poor metabolizer - codeine ineffective, tamoxifen reduced"}
        },
        "evidence": "strong",
        "interventions": ["Alternative pain meds to codeine", "Consider aromatase inhibitors over tamoxifen"]
    },
    "rs3892097": {
        "gene": "CYP2D6",
        "category": "pharmacogenomics",
        "name": "CYP2D6*4",
        "risk_allele": "A",
        "effect": "Non-functional allele - most common cause of poor metabolizer status",
        "genotypes": {
            "GG": {"risk": "normal", "description": "Normal metabolizer"},
            "GA": {"risk": "moderate", "description": "Intermediate metabolizer"},
            "AA": {"risk": "elevated", "description": "Poor metabolizer - many drugs affected"}
        },
        "evidence": "strong",
        "interventions": ["Avoid codeine/tramadol", "Antidepressant dose adjustments", "Pharmacogenomic card"]
    },
    "rs9923231": {
        "gene": "VKORC1",
        "category": "pharmacogenomics",
        "name": "VKORC1 -1639G>A",
        "risk_allele": "A",
        "effect": "Affects warfarin sensitivity - key for dosing",
        "genotypes": {
            "GG": {"risk": "normal", "description": "Normal warfarin dose needed"},
            "GA": {"risk": "moderate", "description": "Intermediate sensitivity - lower dose"},
            "AA": {"risk": "elevated", "description": "High sensitivity - much lower dose needed"}
        },
        "evidence": "strong",
        "interventions": ["Genetic-guided warfarin dosing", "More frequent INR monitoring initially"]
    },
    "rs1800460": {
        "gene": "TPMT",
        "category": "pharmacogenomics",
        "name": "TPMT*3B",
        "risk_allele": "A",
        "effect": "Reduced thiopurine metabolism - azathioprine, 6-MP toxicity risk",
        "genotypes": {
            "GG": {"risk": "normal", "description": "Normal TPMT activity"},
            "GA": {"risk": "elevated", "description": "Intermediate - 50% dose reduction"},
            "AA": {"risk": "high", "description": "Poor metabolizer - severe toxicity risk, avoid thiopurines"}
        },
        "evidence": "strong",
        "interventions": ["TPMT testing before thiopurine therapy", "Dose reduction or alternative drugs"]
    },
    "rs1142345": {
        "gene": "TPMT",
        "category": "pharmacogenomics",
        "name": "TPMT*3C",
        "risk_allele": "C",
        "effect": "Reduced TPMT activity",
        "genotypes": {
            "AA": {"risk": "normal", "description": "Normal TPMT activity"},
            "AC": {"risk": "elevated", "description": "Intermediate metabolizer"},
            "CC": {"risk": "high", "description": "Poor metabolizer - thiopurine toxicity risk"}
        },
        "evidence": "strong",
        "interventions": ["Same as TPMT*3B"]
    },
    "rs4680": {
        "gene": "COMT",
        "category": "pharmacogenomics",
        "name": "COMT Val158Met",
        "risk_allele": "A",
        "effect": "Affects dopamine/catecholamine metabolism - pain sensitivity, stress response",
        "genotypes": {
            "GG": {"risk": "normal", "description": "Val/Val - faster dopamine clearance, stress resilient"},
            "GA": {"risk": "low", "description": "Val/Met - intermediate"},
            "AA": {"risk": "moderate", "description": "Met/Met - slower clearance, higher pain sensitivity, anxiety risk"}
        },
        "evidence": "strong",
        "interventions": ["Stress management", "May need lower opioid doses", "Magnesium supplementation"]
    },

    # === NUTRIENT METABOLISM ===
    "rs4988235": {
        "gene": "MCM6/LCT",
        "category": "nutrient",
        "name": "Lactase persistence",
        "risk_allele": "G",
        "effect": "Determines ability to digest lactose in adulthood",
        "genotypes": {
            "AA": {"risk": "elevated", "description": "Lactose intolerant - lactase non-persistent"},
            "AG": {"risk": "normal", "description": "Lactose tolerant - one persistence allele"},
            "GG": {"risk": "normal", "description": "Lactose tolerant - lactase persistent"}
        },
        "evidence": "strong",
        "interventions": ["Lactase supplements if intolerant", "Calcium from non-dairy sources"]
    },
    "rs602662": {
        "gene": "FUT2",
        "category": "nutrient",
        "name": "FUT2 secretor status",
        "risk_allele": "A",
        "effect": "Affects B12 absorption and gut microbiome",
        "genotypes": {
            "GG": {"risk": "normal", "description": "Secretor - normal B12 absorption"},
            "GA": {"risk": "normal", "description": "Secretor"},
            "AA": {"risk": "moderate", "description": "Non-secretor - reduced B12 absorption, different microbiome"}
        },
        "evidence": "moderate",
        "interventions": ["Monitor B12 levels", "Consider methylcobalamin supplementation"]
    },
    "rs1799945": {
        "gene": "HFE",
        "category": "nutrient",
        "name": "HFE H63D",
        "risk_allele": "G",
        "effect": "Affects iron absorption - hereditary hemochromatosis risk",
        "genotypes": {
            "CC": {"risk": "normal", "description": "Normal iron absorption"},
            "CG": {"risk": "low", "description": "Carrier - mild increased absorption"},
            "GG": {"risk": "moderate", "description": "Homozygous - monitor iron levels"}
        },
        "evidence": "strong",
        "interventions": ["Regular ferritin monitoring", "Avoid iron supplements unless deficient", "Blood donation"]
    },
    "rs1800562": {
        "gene": "HFE",
        "category": "nutrient",
        "name": "HFE C282Y",
        "risk_allele": "A",
        "effect": "Major hereditary hemochromatosis mutation",
        "genotypes": {
            "GG": {"risk": "normal", "description": "No C282Y mutation"},
            "GA": {"risk": "moderate", "description": "Carrier - monitor iron, especially with H63D"},
            "AA": {"risk": "high", "description": "Homozygous - hemochromatosis risk, regular phlebotomy may be needed"}
        },
        "evidence": "strong",
        "interventions": ["Ferritin and transferrin saturation monitoring", "Therapeutic phlebotomy if elevated", "Avoid vitamin C with meals"]
    },
    "rs12934922": {
        "gene": "BCMO1",
        "category": "nutrient",
        "name": "BCMO1 beta-carotene conversion",
        "risk_allele": "T",
        "effect": "Reduced conversion of beta-carotene to vitamin A",
        "genotypes": {
            "AA": {"risk": "normal", "description": "Normal conversion"},
            "AT": {"risk": "moderate", "description": "~30% reduced conversion"},
            "TT": {"risk": "elevated", "description": "~60% reduced conversion - need preformed vitamin A"}
        },
        "evidence": "moderate",
        "interventions": ["Include preformed vitamin A (retinol) from animal sources", "Cod liver oil", "Egg yolks"]
    },
    "rs7946": {
        "gene": "PEMT",
        "category": "nutrient",
        "name": "PEMT choline synthesis",
        "risk_allele": "T",
        "effect": "Reduced endogenous choline production",
        "genotypes": {
            "CC": {"risk": "normal", "description": "Normal choline synthesis"},
            "CT": {"risk": "low", "description": "Slightly reduced synthesis"},
            "TT": {"risk": "moderate", "description": "Reduced synthesis - higher dietary choline needs"}
        },
        "evidence": "moderate",
        "interventions": ["Eggs (high choline)", "Liver", "Choline supplementation"]
    },
    "rs2282679": {
        "gene": "GC",
        "category": "nutrient",
        "name": "Vitamin D binding protein",
        "risk_allele": "C",
        "effect": "Affects vitamin D transport and availability",
        "genotypes": {
            "AA": {"risk": "normal", "description": "Normal vitamin D binding"},
            "AC": {"risk": "moderate", "description": "Lower circulating 25(OH)D"},
            "CC": {"risk": "elevated", "description": "Lower vitamin D levels - higher supplementation needs"}
        },
        "evidence": "strong",
        "interventions": ["Regular vitamin D testing", "Higher supplementation doses", "Sun exposure"]
    },
    "rs12785878": {
        "gene": "DHCR7/NADSYN1",
        "category": "nutrient",
        "name": "Vitamin D synthesis",
        "risk_allele": "G",
        "effect": "Affects vitamin D production in skin",
        "genotypes": {
            "TT": {"risk": "normal", "description": "Normal vitamin D synthesis"},
            "GT": {"risk": "low", "description": "Slightly reduced synthesis"},
            "GG": {"risk": "moderate", "description": "Reduced synthesis - supplementation important"}
        },
        "evidence": "moderate",
        "interventions": ["Vitamin D supplementation", "More sun exposure needed"]
    },

    # === NEUROLOGICAL ===
    "rs6265": {
        "gene": "BDNF",
        "category": "neuro",
        "name": "BDNF Val66Met",
        "risk_allele": "T",
        "effect": "Affects brain-derived neurotrophic factor secretion - memory, mood",
        "genotypes": {
            "CC": {"risk": "normal", "description": "Val/Val - normal BDNF secretion"},
            "CT": {"risk": "low", "description": "Val/Met - reduced activity-dependent secretion"},
            "TT": {"risk": "moderate", "description": "Met/Met - reduced BDNF, depression/anxiety risk, memory effects"}
        },
        "evidence": "strong",
        "interventions": ["Regular exercise (increases BDNF)", "Omega-3 fatty acids", "Stress reduction", "Sleep optimization"]
    },
    "rs1800497": {
        "gene": "ANKK1/DRD2",
        "category": "neuro",
        "name": "DRD2 Taq1A",
        "risk_allele": "T",
        "effect": "Affects dopamine receptor density - reward sensitivity, addiction risk",
        "genotypes": {
            "CC": {"risk": "normal", "description": "A2/A2 - normal dopamine receptors"},
            "CT": {"risk": "moderate", "description": "A1/A2 - reduced receptor density"},
            "TT": {"risk": "elevated", "description": "A1/A1 - ~30% fewer D2 receptors, addiction vulnerability"}
        },
        "evidence": "moderate",
        "interventions": ["Avoid addictive substances", "Behavioral addiction awareness", "Dopamine-supporting supplements"]
    },
    "rs6311": {
        "gene": "HTR2A",
        "category": "neuro",
        "name": "Serotonin receptor 2A",
        "risk_allele": "T",
        "effect": "Affects serotonin signaling - antidepressant response, mood",
        "genotypes": {
            "CC": {"risk": "normal", "description": "Normal serotonin receptor expression"},
            "CT": {"risk": "low", "description": "Altered expression"},
            "TT": {"risk": "moderate", "description": "Different antidepressant response profile"}
        },
        "evidence": "moderate",
        "interventions": ["May respond differently to SSRIs", "Tryptophan-rich foods"]
    },
    "rs25531": {
        "gene": "SLC6A4",
        "category": "neuro",
        "name": "Serotonin transporter (5-HTTLPR related)",
        "risk_allele": "G",
        "effect": "Affects serotonin reuptake - stress response, depression risk",
        "genotypes": {
            "AA": {"risk": "normal", "description": "Long/Long - resilient to stress"},
            "AG": {"risk": "moderate", "description": "Long/Short - intermediate stress sensitivity"},
            "GG": {"risk": "elevated", "description": "Short/Short - higher stress sensitivity, depression risk"}
        },
        "evidence": "moderate",
        "interventions": ["Stress management essential", "Strong social support", "Mindfulness practices"]
    }
}


def load_ancestry_data(filepath: str) -> pd.DataFrame:
    """Load AncestryDNA raw data file."""
    print(f"Loading DNA data from {filepath}...")

    # Skip header comments
    df = pd.read_csv(filepath, sep='\t', comment='#',
                     names=['rsid', 'chromosome', 'position', 'allele1', 'allele2'])

    # Create genotype column
    df['genotype'] = df['allele1'] + df['allele2']

    print(f"Loaded {len(df):,} SNPs")
    return df


def analyze_clinical_variants(df: pd.DataFrame) -> List[Dict]:
    """Find and annotate clinically significant variants."""
    results = []

    rsids_in_data = set(df['rsid'].values)

    for rsid, info in CLINICAL_SNPS.items():
        if rsid in rsids_in_data:
            row = df[df['rsid'] == rsid].iloc[0]
            genotype = row['genotype']

            # Normalize genotype (sort alleles for matching)
            genotype_sorted = ''.join(sorted(genotype))
            genotype_options = [genotype, genotype_sorted, genotype[::-1]]

            # Find matching genotype interpretation
            interpretation = None
            matched_genotype = None
            for gt in genotype_options:
                if gt in info['genotypes']:
                    interpretation = info['genotypes'][gt]
                    matched_genotype = gt
                    break

            # Check for reversed notation
            if not interpretation:
                for gt_key, gt_info in info['genotypes'].items():
                    if set(genotype) == set(gt_key):
                        interpretation = gt_info
                        matched_genotype = gt_key
                        break

            result = {
                'rsid': rsid,
                'gene': info['gene'],
                'name': info['name'],
                'category': info['category'],
                'your_genotype': genotype,
                'matched_genotype': matched_genotype,
                'risk_allele': info['risk_allele'],
                'risk_level': interpretation['risk'] if interpretation else 'unknown',
                'description': interpretation['description'] if interpretation else 'Genotype not in database',
                'effect': info['effect'],
                'evidence': info['evidence'],
                'interventions': info['interventions'],
                'chromosome': row['chromosome'],
                'position': row['position']
            }
            results.append(result)

    return results


def determine_apoe_status(results: List[Dict]) -> Optional[str]:
    """Determine APOE genotype from rs429358 and rs7412."""
    rs429358 = None
    rs7412 = None

    for r in results:
        if r['rsid'] == 'rs429358':
            rs429358 = r['your_genotype']
        elif r['rsid'] == 'rs7412':
            rs7412 = r['your_genotype']

    if not rs429358 or not rs7412:
        return None

    # APOE allele determination:
    # e2: rs429358=T, rs7412=T
    # e3: rs429358=T, rs7412=C
    # e4: rs429358=C, rs7412=C

    apoe_alleles = []
    for i in range(2):
        a1 = rs429358[i] if i < len(rs429358) else rs429358[0]
        a2 = rs7412[i] if i < len(rs7412) else rs7412[0]

        if a1 == 'T' and a2 == 'T':
            apoe_alleles.append('e2')
        elif a1 == 'T' and a2 == 'C':
            apoe_alleles.append('e3')
        elif a1 == 'C' and a2 == 'C':
            apoe_alleles.append('e4')
        else:
            apoe_alleles.append('unknown')

    return f"{apoe_alleles[0]}/{apoe_alleles[1]}"


def query_clinvar(rsid: str) -> Optional[Dict]:
    """Query NCBI ClinVar for variant information."""
    try:
        # Use NCBI E-utilities
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

        # Search for the rsID
        search_url = f"{base_url}esearch.fcgi?db=clinvar&term={rsid}&retmode=json"
        response = requests.get(search_url, timeout=10)

        if response.status_code == 200:
            data = response.json()
            if data.get('esearchresult', {}).get('count', '0') != '0':
                return {'found': True, 'source': 'ClinVar'}

        return None
    except Exception as e:
        return None


def categorize_results(results: List[Dict]) -> Dict[str, List[Dict]]:
    """Categorize results by risk level and category."""
    high_priority = []
    moderate_priority = []
    low_priority = []

    for r in results:
        risk = r['risk_level']
        if risk in ['high', 'elevated']:
            high_priority.append(r)
        elif risk == 'moderate':
            moderate_priority.append(r)
        else:
            low_priority.append(r)

    return {
        'high_priority': sorted(high_priority, key=lambda x: (x['evidence'] != 'strong', x['category'])),
        'moderate_priority': sorted(moderate_priority, key=lambda x: (x['evidence'] != 'strong', x['category'])),
        'low_priority': sorted(low_priority, key=lambda x: x['category'])
    }


def generate_report(results: List[Dict], apoe_status: Optional[str]) -> str:
    """Generate a formatted report."""
    categorized = categorize_results(results)

    report = []
    report.append("=" * 80)
    report.append("ANCESTRYDNA CLINICAL VARIANT ANALYSIS REPORT")
    report.append("=" * 80)
    report.append("")

    if apoe_status:
        report.append(f"APOE STATUS: {apoe_status}")
        if 'e4' in apoe_status:
            report.append("  ⚠️  APOE e4 carrier - associated with increased cardiovascular and Alzheimer's risk")
        report.append("")

    # Summary
    report.append("SUMMARY")
    report.append("-" * 40)
    report.append(f"Total clinical variants analyzed: {len(results)}")
    report.append(f"High priority findings: {len(categorized['high_priority'])}")
    report.append(f"Moderate priority findings: {len(categorized['moderate_priority'])}")
    report.append(f"Low priority/normal findings: {len(categorized['low_priority'])}")
    report.append("")

    # High priority
    if categorized['high_priority']:
        report.append("🚨 HIGH PRIORITY FINDINGS (Requires Attention)")
        report.append("=" * 60)
        for r in categorized['high_priority']:
            report.append(f"\n{r['gene']} - {r['name']}")
            report.append(f"  rsID: {r['rsid']}")
            report.append(f"  Your genotype: {r['your_genotype']}")
            report.append(f"  Risk level: {r['risk_level'].upper()}")
            report.append(f"  Interpretation: {r['description']}")
            report.append(f"  Effect: {r['effect']}")
            report.append(f"  Evidence: {r['evidence']}")
            report.append(f"  Recommended actions:")
            for intervention in r['interventions']:
                report.append(f"    • {intervention}")
        report.append("")

    # Moderate priority
    if categorized['moderate_priority']:
        report.append("⚠️  MODERATE PRIORITY FINDINGS (Monitor/Consider)")
        report.append("=" * 60)
        for r in categorized['moderate_priority']:
            report.append(f"\n{r['gene']} - {r['name']}")
            report.append(f"  rsID: {r['rsid']}")
            report.append(f"  Your genotype: {r['your_genotype']}")
            report.append(f"  Risk level: {r['risk_level']}")
            report.append(f"  Interpretation: {r['description']}")
            report.append(f"  Effect: {r['effect']}")
            report.append(f"  Recommended actions:")
            for intervention in r['interventions']:
                report.append(f"    • {intervention}")
        report.append("")

    # By category summary
    report.append("FINDINGS BY CATEGORY")
    report.append("=" * 60)

    categories = {}
    for r in results:
        cat = r['category']
        if cat not in categories:
            categories[cat] = []
        categories[cat].append(r)

    category_names = {
        'cardiovascular': '❤️  Cardiovascular',
        'metabolic': '🔥 Metabolic/Diabetes',
        'cancer': '🎗️  Cancer Risk',
        'pharmacogenomics': '💊 Drug Response',
        'nutrient': '🥗 Nutrient Metabolism',
        'neuro': '🧠 Neurological'
    }

    for cat, name in category_names.items():
        if cat in categories:
            report.append(f"\n{name}")
            report.append("-" * 40)
            for r in categories[cat]:
                risk_icon = "🔴" if r['risk_level'] in ['high', 'elevated'] else "🟡" if r['risk_level'] == 'moderate' else "🟢"
                report.append(f"  {risk_icon} {r['gene']} ({r['rsid']}): {r['your_genotype']} - {r['risk_level']}")

    report.append("")
    report.append("=" * 80)
    report.append("DISCLAIMER")
    report.append("=" * 80)
    report.append("""
This analysis is for educational purposes only and is NOT medical advice.
Genetic variants interact in complex ways, and environmental factors play major roles.
Many variants have incomplete penetrance - having a risk allele does not mean you will
develop a condition. Always consult healthcare providers for medical decisions.
Consider professional genetic counseling for significant findings.
""")

    return "\n".join(report)


def export_to_json(results: List[Dict], apoe_status: Optional[str], filepath: str):
    """Export results to JSON for further analysis."""
    output = {
        'apoe_status': apoe_status,
        'summary': {
            'total_variants': len(results),
            'high_priority': len([r for r in results if r['risk_level'] in ['high', 'elevated']]),
            'moderate_priority': len([r for r in results if r['risk_level'] == 'moderate'])
        },
        'variants': results
    }

    with open(filepath, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"Results exported to {filepath}")


def main():
    import sys

    # Default file path
    filepath = "AncestryDNA.txt"
    if len(sys.argv) > 1:
        filepath = sys.argv[1]

    # Load data
    df = load_ancestry_data(filepath)

    # Analyze clinical variants
    print("\nAnalyzing clinical variants...")
    results = analyze_clinical_variants(df)

    # Determine APOE status
    apoe_status = determine_apoe_status(results)

    # Generate report
    report = generate_report(results, apoe_status)
    print(report)

    # Export to JSON
    export_to_json(results, apoe_status, "clinical_variants_report.json")

    # Save report to file
    with open("clinical_variants_report.txt", 'w') as f:
        f.write(report)
    print("\nReport saved to clinical_variants_report.txt")

    return results, apoe_status


if __name__ == "__main__":
    main()
