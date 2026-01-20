"""
Prompt templates for AI-powered genetic report generation.

These prompts are carefully designed to generate accurate, helpful,
and appropriately cautious genetic health reports.
"""

SYSTEM_PROMPT = """You are a knowledgeable genetic counselor AI assistant helping users 
understand their personal genomics data. Your role is to:

1. Explain genetic findings in clear, accurate language
2. Provide actionable health insights where appropriate
3. Always include appropriate medical disclaimers
4. Emphasize that genetics is just one factor in health outcomes
5. Recommend professional genetic counseling for significant findings
6. Never make definitive medical diagnoses

You are analyzing results from GenomeInsight, a privacy-first personal genomics toolkit.
All analysis runs locally on the user's machine - their data never leaves their computer.

Important guidelines:
- Be scientifically accurate but accessible
- Highlight both risks and protective factors
- Provide context about penetrance and environmental factors
- Suggest evidence-based lifestyle interventions when relevant
- Always recommend consulting healthcare providers for medical decisions"""


TECHNICAL_TEMPLATE = """Generate a technical genetic analysis report for a researcher or healthcare provider.

## Analysis Summary
- Source file: {source_file}
- Total SNPs analyzed: {snp_count:,}
- Clinical variants database size: {variants_analyzed}
- Variants found in subject: {variants_found}

## APOE Status
{apoe_section}

## High Priority Findings ({high_count})
{high_priority_section}

## Moderate Priority Findings ({moderate_count})
{moderate_priority_section}

## Gene Interactions Detected ({interaction_count})
{interactions_section}

## Category Breakdown
{category_breakdown}

---

Please generate a comprehensive technical report including:
1. Executive summary of key findings
2. Detailed variant analysis with clinical significance
3. Gene-gene interaction effects
4. Pharmacogenomic considerations
5. Risk stratification assessment
6. Recommended follow-up testing
7. Relevant literature citations (use general references, not specific papers)

Format the report in markdown with clear section headers."""


CONSUMER_TEMPLATE = """Generate a plain-language genetic health report for a consumer with no genetics background.

## Your DNA Analysis Results
- We analyzed {snp_count:,} genetic markers from your DNA
- We checked {variants_analyzed} health-related genetic variants
- We found {variants_found} variants to discuss with you

## Your APOE Status (Important for Brain Health)
{apoe_section}

## Key Findings That Need Attention ({high_count})
{high_priority_section}

## Other Notable Findings ({moderate_count})
{moderate_priority_section}

## Important Gene Combinations Found ({interaction_count})
{interactions_section}

---

Please generate a friendly, easy-to-understand report that:
1. Starts with a warm, reassuring introduction
2. Explains what the key findings mean in everyday language
3. Avoids medical jargon (or explains it when necessary)
4. Provides practical, actionable lifestyle suggestions
5. Clearly states what factors you CAN control
6. Emphasizes that genes are not destiny
7. Recommends when to speak with a doctor
8. Ends with an encouraging, empowering message

Use simple language, short paragraphs, and bullet points for readability.
Include emojis sparingly for visual appeal (❤️ 🧬 💪 🥗 etc).
Format the report in markdown."""


def format_apoe_section(apoe_status) -> str:
    """Format APOE status for prompt."""
    if not apoe_status:
        return "APOE genotype: Not determined (required SNPs not found in data)"
    
    return f"""APOE Genotype: {apoe_status.genotype}
Risk Category: {apoe_status.risk_category}
Has ε4 allele: {"Yes" if apoe_status.has_e4 else "No"}
Has ε2 allele: {"Yes" if apoe_status.has_e2 else "No"}
ε4 homozygous: {"Yes" if apoe_status.is_e4_homozygous else "No"}
Interpretation: {apoe_status.interpretation}
Recommendations: {', '.join(apoe_status.recommendations)}"""


def format_variant_section(variants: list) -> str:
    """Format variant list for prompt."""
    if not variants:
        return "None found."
    
    lines = []
    for v in variants:
        lines.append(f"""
- **{v.gene}** ({v.rsid}): {v.name}
  - Genotype: {v.genotype}
  - Risk Level: {v.risk_level.value}
  - Category: {v.category.value}
  - Description: {v.description}
  - Evidence: {v.evidence.value}
  - Suggested interventions: {', '.join(v.interventions) if v.interventions else 'None specified'}""")
    
    return "\n".join(lines)


def format_interactions_section(interactions: list) -> str:
    """Format gene interactions for prompt."""
    if not interactions:
        return "No significant gene-gene interactions detected."
    
    lines = []
    for i in interactions:
        genes = " + ".join(i.get("genes", []))
        lines.append(f"""
- **{genes}** ({i.get('type', 'interaction')})
  - Severity: {i.get('severity', 'unknown')}
  - Note: {i.get('note', 'No details')}
  - Recommendations: {', '.join(i.get('recommendations', []))}""")
    
    return "\n".join(lines)


def format_category_breakdown(by_category: dict) -> str:
    """Format category breakdown for prompt."""
    if not by_category:
        return "No variants categorized."
    
    lines = []
    for category, variants in by_category.items():
        lines.append(f"- {category.value}: {len(variants)} variants")
    
    return "\n".join(lines)


def build_prompt(result, style: str = "consumer") -> str:
    """
    Build the complete prompt for report generation.
    
    Args:
        result: AnalysisResult from clinical analyzer
        style: 'technical' or 'consumer'
        
    Returns:
        Formatted prompt string
    """
    template = CONSUMER_TEMPLATE if style == "consumer" else TECHNICAL_TEMPLATE
    
    return template.format(
        source_file=result.source_file,
        snp_count=result.snp_count,
        variants_analyzed=result.variants_analyzed,
        variants_found=result.variants_found,
        apoe_section=format_apoe_section(result.apoe_status),
        high_count=len(result.high_priority),
        high_priority_section=format_variant_section(result.high_priority),
        moderate_count=len(result.moderate_priority),
        moderate_priority_section=format_variant_section(result.moderate_priority),
        interaction_count=len(result.gene_interactions),
        interactions_section=format_interactions_section(result.gene_interactions),
        category_breakdown=format_category_breakdown(result.by_category),
    )
