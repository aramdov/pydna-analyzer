#!/usr/bin/env python3
"""
compare_promethease.py

Compares pydna-analyzer output against Promethease to evaluate:
- RSID coverage overlap
- Risk classification accuracy (where both tools have data)
- Category coverage gaps
- Bad-repute variants that pydna-analyzer misses

Optionally uses Claude SDK (--ai) to:
- Suggest correct risk levels for pydna variants classified as "unknown"
- Assess whether Promethease-only bad-repute variants should be added to pydna-analyzer

Usage:
    python compare_promethease.py
    python compare_promethease.py --ai
    python compare_promethease.py --ai --model claude-haiku-4-5-20251001
    python compare_promethease.py --output report.md
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

# ─── Paths ────────────────────────────────────────────────────────────────────

BASE_DIR = Path(__file__).parent
PROMETHEASE_CSV = BASE_DIR.parent / "promethease" / "Promethease Table.csv"
DEFAULT_DNA_FILE = BASE_DIR / "AncestryDNA.txt"


# ─── Data models ──────────────────────────────────────────────────────────────

@dataclass
class PydnaVariant:
    rsid: str
    gene: str
    name: str
    category: str
    genotype: str
    risk_level: str       # unknown / normal / low / moderate / high
    risk_allele: str
    description: str
    effect: str
    evidence: str
    interventions: list[str] = field(default_factory=list)


@dataclass
class PrometheaseVariant:
    rsid: str
    name: str             # original e.g. "rs1801133(C;T)"
    magnitude: float
    repute: str           # Good / Bad / (blank)
    summary: str
    gene: str
    chrom: str
    pos: str
    publications: str
    freq: str


@dataclass
class ComparisonRow:
    rsid: str
    gene: str
    pydna_risk_level: str
    pydna_genotype: str
    pydna_description: str
    prom_magnitude: float
    prom_repute: str
    prom_summary: str
    status: str           # overlap / pydna_only / prom_only
    mismatch: bool = False
    mismatch_note: str = ""
    ai_suggestion: str = ""


# ─── Loaders ──────────────────────────────────────────────────────────────────

def _variants_from_analysis_json(data: dict) -> list[PydnaVariant]:
    variants = []
    for v in data.get("variant_results", []):
        variants.append(PydnaVariant(
            rsid=v["rsid"],
            gene=v.get("gene", ""),
            name=v.get("name", ""),
            category=v.get("category", ""),
            genotype=v.get("genotype", ""),
            risk_level=v.get("risk_level", "unknown"),
            risk_allele=v.get("risk_allele", ""),
            description=v.get("description", ""),
            effect=v.get("effect", ""),
            evidence=v.get("evidence", ""),
            interventions=v.get("interventions", []),
        ))
    return variants


def load_pydna_from_json(path: Path) -> list[PydnaVariant]:
    with open(path) as f:
        data = json.load(f)
    return _variants_from_analysis_json(data)


def load_pydna_from_dna_file(path: Path) -> list[PydnaVariant]:
    from pydna_analyzer.clinical.analyzer import ClinicalAnalyzer
    from pydna_analyzer.core.data_loader import load_dna_data

    dataset = load_dna_data(path)
    result = ClinicalAnalyzer().analyze(dataset)

    return [
        PydnaVariant(
            rsid=v.rsid,
            gene=v.gene,
            name=v.name,
            category=v.category.value,
            genotype=v.genotype,
            risk_level=v.risk_level.value,
            risk_allele=v.risk_allele,
            description=v.description,
            effect=v.effect,
            evidence=v.evidence.value,
            interventions=v.interventions,
        )
        for v in result.variant_results
    ]


def load_promethease(path: Path) -> list[PrometheaseVariant]:
    variants = []
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get("Name", "")
            m = re.match(r"(rs\d+)", name)
            if not m:
                continue
            rsid = m.group(1)
            try:
                mag = float(row.get("Magnitude", "0") or "0")
            except ValueError:
                mag = 0.0
            variants.append(PrometheaseVariant(
                rsid=rsid,
                name=name,
                magnitude=mag,
                repute=row.get("Repute", "").strip(),
                summary=row.get("Summary", "").strip(),
                gene=row.get("Gene", "").strip(),
                chrom=row.get("Chrom", "").strip(),
                pos=row.get("Pos", "").strip(),
                publications=row.get("Publications", "").strip(),
                freq=row.get("Freq", "").strip(),
            ))
    return variants


# ─── Comparison logic ─────────────────────────────────────────────────────────

# Rough magnitude → risk_level mapping for mismatch detection
def magnitude_to_risk(mag: float) -> str:
    if mag < 1.0:
        return "normal"
    elif mag < 2.0:
        return "low"
    elif mag < 3.0:
        return "moderate"
    else:
        return "high"


def is_mismatch(pydna_risk: str, prom_repute: str, prom_mag: float) -> tuple[bool, str]:
    """
    Flag mismatches between pydna risk classification and Promethease repute/magnitude.

    Key rule: if Promethease says Bad with magnitude >= 2, pydna should NOT be
    "unknown" or "normal". That gap is the signal.
    """
    if prom_repute == "Bad" and prom_mag >= 2.0:
        if pydna_risk in ("unknown", "normal"):
            return True, f"Promethease: Bad/mag={prom_mag}, pydna: {pydna_risk} — likely under-classified"
    if prom_repute == "Good" and pydna_risk in ("moderate", "high"):
        return True, f"Promethease: Good/mag={prom_mag}, pydna: {pydna_risk} — possible over-classification"
    return False, ""


def build_comparison(
    pydna_variants: list[PydnaVariant],
    prom_variants: list[PrometheaseVariant],
) -> list[ComparisonRow]:
    pydna_map = {v.rsid: v for v in pydna_variants}
    prom_map = {v.rsid: v for v in prom_variants}

    rows: list[ComparisonRow] = []

    # Overlapping RSIDs
    overlap = set(pydna_map) & set(prom_map)
    for rsid in sorted(overlap):
        p = pydna_map[rsid]
        pr = prom_map[rsid]
        mismatch, note = is_mismatch(p.risk_level, pr.repute, pr.magnitude)
        rows.append(ComparisonRow(
            rsid=rsid,
            gene=p.gene,
            pydna_risk_level=p.risk_level,
            pydna_genotype=p.genotype,
            pydna_description=p.description,
            prom_magnitude=pr.magnitude,
            prom_repute=pr.repute,
            prom_summary=pr.summary,
            status="overlap",
            mismatch=mismatch,
            mismatch_note=note,
        ))

    # pydna-only
    for rsid in sorted(set(pydna_map) - set(prom_map)):
        p = pydna_map[rsid]
        rows.append(ComparisonRow(
            rsid=rsid,
            gene=p.gene,
            pydna_risk_level=p.risk_level,
            pydna_genotype=p.genotype,
            pydna_description=p.description,
            prom_magnitude=0.0,
            prom_repute="N/A",
            prom_summary="Not in Promethease output",
            status="pydna_only",
        ))

    # Promethease Bad-repute only (not tracked by pydna)
    for rsid in sorted(set(prom_map) - set(pydna_map)):
        pr = prom_map[rsid]
        if pr.repute == "Bad" and pr.magnitude >= 1.0:
            rows.append(ComparisonRow(
                rsid=rsid,
                gene=pr.gene,
                pydna_risk_level="not tracked",
                pydna_genotype="N/A",
                pydna_description="pydna-analyzer does not cover this variant",
                prom_magnitude=pr.magnitude,
                prom_repute=pr.repute,
                prom_summary=pr.summary,
                status="prom_only",
            ))

    return rows


# ─── Summary stats ────────────────────────────────────────────────────────────

def print_summary(
    pydna_variants: list[PydnaVariant],
    prom_variants: list[PrometheaseVariant],
    rows: list[ComparisonRow],
) -> None:
    overlap = [r for r in rows if r.status == "overlap"]
    pydna_only = [r for r in rows if r.status == "pydna_only"]
    prom_only_bad = [r for r in rows if r.status == "prom_only"]
    mismatches = [r for r in overlap if r.mismatch]
    unknowns = [r for r in pydna_only if r.pydna_risk_level == "unknown"]

    print("=" * 60)
    print("PYDNA-ANALYZER vs. PROMETHEASE — COMPARISON SUMMARY")
    print("=" * 60)
    print(f"\n{'pydna-analyzer variants found:':<40} {len(pydna_variants)}")
    print(f"{'Promethease SNPs in table:':<40} {len(prom_variants)}")
    print(f"\n{'Overlapping RSIDs:':<40} {len(overlap)}")
    print(f"{'  → Classification mismatches:':<40} {len(mismatches)}")
    print(f"{'pydna-only (not in Promethease):':<40} {len(pydna_only)}")
    print(f"{'  → risk_level=unknown in pydna:':<40} {len(unknowns)}")
    print(f"{'Promethease Bad-repute gaps:':<40} {len(prom_only_bad)}")

    if mismatches:
        print("\n── CLASSIFICATION MISMATCHES ─────────────────────────")
        for r in mismatches:
            print(f"  {r.rsid} ({r.gene}): {r.mismatch_note}")

    if prom_only_bad:
        print("\n── PROMETHEASE BAD-REPUTE VARIANTS NOT IN PYDNA ─────")
        for r in sorted(prom_only_bad, key=lambda x: -x.prom_magnitude):
            print(f"  {r.rsid} ({r.gene}) mag={r.prom_magnitude}: {r.prom_summary[:70]}")

    if unknowns:
        print("\n── PYDNA 'UNKNOWN' RISK LEVELS (candidates for AI fix) ─")
        for r in unknowns:
            print(f"  {r.rsid} ({r.gene}) genotype={r.pydna_genotype}: {r.pydna_description[:60]}")

    print()


# ─── AI enrichment ────────────────────────────────────────────────────────────

AI_SYSTEM_PROMPT = """You are a genomics expert assistant. You help evaluate genetic variant classifications
from a personal genomics tool called pydna-analyzer.

When given a variant, you assess whether the tool's risk classification is correct based on
established literature (SNPedia, ClinVar, published GWAS). Be concise and actionable.
Do not provide medical advice — frame everything as informational/educational."""


def build_risk_classification_prompt(row: ComparisonRow, pydna_v: PydnaVariant) -> str:
    return f"""Variant: {row.rsid} | Gene: {row.gene}
Genotype observed: {row.pydna_genotype}
Risk allele: {pydna_v.risk_allele}
Known effect: {pydna_v.effect}
Current pydna-analyzer classification: risk_level = "{row.pydna_risk_level}" ("{row.pydna_description}")

Task: Based on established literature (SNPedia, ClinVar, GWAS), what risk_level should be assigned?
Choose from: normal / low / moderate / high

Reply in this format:
SUGGESTED_RISK_LEVEL: <level>
RATIONALE: <1-2 sentence explanation>
"""


def build_gap_assessment_prompt(row: ComparisonRow) -> str:
    return f"""Variant: {row.rsid} | Gene: {row.gene}
Promethease magnitude: {row.prom_magnitude} | Repute: {row.repute if hasattr(row, 'repute') else row.prom_repute}
Summary: {row.prom_summary}

This variant appears in Promethease as "Bad" repute but is NOT tracked by pydna-analyzer.
pydna-analyzer focuses on cardiovascular, metabolic, cancer, pharmacogenomics, nutrient, and neuro categories.

Task: Should pydna-analyzer add this variant?
Reply in this format:
SHOULD_ADD: yes / no / maybe
CATEGORY: <which pydna category it would fit>
RATIONALE: <1-2 sentences>
"""


def _call_claude_cli(prompt: str, model: str, system: str = "") -> str:
    """
    Call Claude headlessly via `claude -p` (uses your Claude subscription).
    No API key needed — authenticates via your local Claude login.
    """
    import subprocess
    full_prompt = f"{system}\n\n{prompt}" if system else prompt
    cmd = ["claude", "-p", full_prompt, "--model", model]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or f"claude -p exited {result.returncode}")
    return result.stdout.strip()


def _call_sdk(prompt: str, model: str, system: str = "") -> str:
    """
    Call Claude via Anthropic SDK (requires ANTHROPIC_API_KEY env var).
    """
    import anthropic
    import httpx
    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        raise RuntimeError("ANTHROPIC_API_KEY not set")
    http_proxy = os.environ.get("HTTPS_PROXY") or os.environ.get("HTTP_PROXY")
    http_client = httpx.Client(proxy=http_proxy) if http_proxy else httpx.Client()
    client = anthropic.Anthropic(api_key=api_key, http_client=http_client)
    messages = [{"role": "user", "content": prompt}]
    kwargs: dict = {"model": model, "max_tokens": 512, "messages": messages}
    if system:
        kwargs["system"] = system
    response = client.messages.create(**kwargs)
    return response.content[0].text.strip()


def _call_ai(prompt: str, model: str, system: str = "", use_cli: bool = True) -> str:
    """Dispatch to claude CLI or SDK based on use_cli flag."""
    if use_cli:
        return _call_claude_cli(prompt, model, system)
    return _call_sdk(prompt, model, system)


def enrich_with_ai(
    rows: list[ComparisonRow],
    pydna_map: dict[str, PydnaVariant],
    model: str = "claude-haiku-4-5-20251001",
    use_cli: bool = True,
) -> list[ComparisonRow]:
    """Enrich comparison rows with Claude suggestions (via CLI or SDK)."""

    # 1. Fix "unknown" risk levels in pydna-only variants
    unknowns = [r for r in rows if r.status == "pydna_only" and r.pydna_risk_level == "unknown"]
    if unknowns:
        print(f"\n[AI] Assessing {len(unknowns)} 'unknown' risk level variants...")
        for r in unknowns:
            pydna_v = pydna_map.get(r.rsid)
            if not pydna_v:
                continue
            try:
                r.ai_suggestion = _call_ai(
                    build_risk_classification_prompt(r, pydna_v), model, AI_SYSTEM_PROMPT, use_cli
                )
                print(f"  {r.rsid} ({r.gene}): {r.ai_suggestion.splitlines()[0]}")
            except Exception as e:
                print(f"  {r.rsid}: error — {e}", file=sys.stderr)

    # 2. Assess Promethease-only Bad-repute gaps
    prom_gaps = [r for r in rows if r.status == "prom_only"]
    if prom_gaps:
        print(f"\n[AI] Assessing {len(prom_gaps)} Promethease Bad-repute gaps...")
        for r in prom_gaps:
            try:
                r.ai_suggestion = _call_ai(
                    build_gap_assessment_prompt(r), model, AI_SYSTEM_PROMPT, use_cli
                )
                print(f"  {r.rsid} ({r.gene}): {r.ai_suggestion.splitlines()[0]}")
            except Exception as e:
                print(f"  {r.rsid}: error — {e}", file=sys.stderr)

    # 3. Adjudicate overlapping mismatches
    mismatches = [r for r in rows if r.status == "overlap" and r.mismatch]
    if mismatches:
        print(f"\n[AI] Adjudicating {len(mismatches)} classification mismatches...")
        for r in mismatches:
            pydna_v = pydna_map.get(r.rsid)
            if not pydna_v:
                continue
            prompt = (
                build_risk_classification_prompt(r, pydna_v)
                + f"\nAdditional context from Promethease:\n"
                + f"  Magnitude: {r.prom_magnitude} | Repute: {r.prom_repute}\n"
                + f"  Summary: {r.prom_summary}\n"
            )
            try:
                r.ai_suggestion = _call_ai(prompt, model, AI_SYSTEM_PROMPT, use_cli)
                print(f"  {r.rsid} ({r.gene}): {r.ai_suggestion.splitlines()[0]}")
            except Exception as e:
                print(f"  {r.rsid}: error — {e}", file=sys.stderr)

    return rows


def generate_ai_narrative(
    rows: list[ComparisonRow],
    pydna_variants: list[PydnaVariant],
    prom_variants: list[PrometheaseVariant],
    model: str = "claude-sonnet-4-6",
    use_cli: bool = True,
) -> str:
    """Ask Claude to write a natural-language comparison narrative."""
    overlap = [r for r in rows if r.status == "overlap"]
    pydna_only = [r for r in rows if r.status == "pydna_only"]
    prom_gaps = [r for r in rows if r.status == "prom_only"]
    mismatches = [r for r in overlap if r.mismatch]
    unknowns_with_ai = [r for r in pydna_only if r.ai_suggestion and "SUGGESTED_RISK_LEVEL" in r.ai_suggestion]

    context = (
        f"pydna-analyzer vs Promethease Comparison Data:\n\n"
        f"Total pydna variants found: {len(pydna_variants)}\n"
        f"Total Promethease SNPs: {len(prom_variants)}\n"
        f"Overlapping RSIDs: {len(overlap)}\n"
        f"Classification mismatches: {len(mismatches)}\n"
        f"pydna-only variants: {len(pydna_only)}\n"
        f"Promethease Bad-repute gaps (not in pydna): {len(prom_gaps)}\n\n"
        f"Mismatches:\n"
        + (
            "\n".join(f"- {r.rsid} ({r.gene}): {r.mismatch_note}" for r in mismatches)
            or "none"
        )
        + f"\n\nPromethease Bad-repute not in pydna:\n"
        + (
            "\n".join(f"- {r.rsid} ({r.gene}) mag={r.prom_magnitude}: {r.prom_summary}" for r in prom_gaps)
            or "none"
        )
        + f"\n\nAI suggestions for unknown risk levels:\n"
        + (
            "\n".join(f"- {r.rsid} ({r.gene}): {r.ai_suggestion}" for r in unknowns_with_ai)
            or "none"
        )
    )

    prompt = (
        "Based on the comparison data below, write a concise (300-400 word) executive summary "
        "evaluating pydna-analyzer's genetic/health analysis quality relative to Promethease.\n\n"
        "Cover: (1) overall coverage assessment, (2) where pydna excels vs lags, "
        "(3) the most important bugs or gaps to fix, (4) recommendations.\n\n"
        + context
    )

    try:
        return _call_ai(prompt, model, AI_SYSTEM_PROMPT, use_cli)
    except Exception as e:
        print(f"  [narrative error: {e}]", file=sys.stderr)
        return ""


# ─── CSV export ───────────────────────────────────────────────────────────────

def export_csv(rows: list[ComparisonRow], output_path: Path) -> None:
    fields = [
        "rsid", "gene", "status", "pydna_risk_level", "pydna_genotype",
        "pydna_description", "prom_magnitude", "prom_repute", "prom_summary",
        "mismatch", "mismatch_note", "ai_suggestion",
    ]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for r in rows:
            writer.writerow({
                "rsid": r.rsid,
                "gene": r.gene,
                "status": r.status,
                "pydna_risk_level": r.pydna_risk_level,
                "pydna_genotype": r.pydna_genotype,
                "pydna_description": r.pydna_description[:100],
                "prom_magnitude": r.prom_magnitude,
                "prom_repute": r.prom_repute,
                "prom_summary": r.prom_summary[:100],
                "mismatch": r.mismatch,
                "mismatch_note": r.mismatch_note,
                "ai_suggestion": r.ai_suggestion.replace("\n", " | ") if r.ai_suggestion else "",
            })
    print(f"\nCSV exported to: {output_path}")


def export_markdown(
    rows: list[ComparisonRow],
    output_path: Path,
    narrative: str = "",
) -> None:
    overlap = [r for r in rows if r.status == "overlap"]
    pydna_only = [r for r in rows if r.status == "pydna_only"]
    prom_gaps = [r for r in rows if r.status == "prom_only"]
    mismatches = [r for r in overlap if r.mismatch]

    lines = ["# pydna-analyzer vs Promethease: Comparison Report\n"]

    if narrative:
        lines.append("## Executive Summary\n")
        lines.append(narrative + "\n")

    lines.append("## Stats\n")
    lines.append(f"| Metric | Count |")
    lines.append(f"|--------|-------|")
    lines.append(f"| pydna-analyzer variants | {len(pydna_only) + len(overlap)} |")
    lines.append(f"| Promethease SNPs | {sum(1 for r in rows)} |")
    lines.append(f"| RSID overlap | {len(overlap)} |")
    lines.append(f"| Classification mismatches | {len(mismatches)} |")
    lines.append(f"| Promethease Bad-repute gaps | {len(prom_gaps)} |")
    lines.append("")

    if mismatches:
        lines.append("## Classification Mismatches\n")
        lines.append("| RSID | Gene | pydna Risk | Prom Mag | Prom Repute | Note |")
        lines.append("|------|------|-----------|---------|------------|------|")
        for r in mismatches:
            lines.append(f"| {r.rsid} | {r.gene} | {r.pydna_risk_level} | {r.prom_magnitude} | {r.prom_repute} | {r.mismatch_note} |")
        lines.append("")

    if prom_gaps:
        lines.append("## Promethease Bad-repute Variants Not in pydna-analyzer\n")
        lines.append("| RSID | Gene | Magnitude | Summary | AI Assessment |")
        lines.append("|------|------|-----------|---------|---------------|")
        for r in sorted(prom_gaps, key=lambda x: -x.prom_magnitude):
            ai = r.ai_suggestion.replace("\n", " | ") if r.ai_suggestion else ""
            lines.append(f"| {r.rsid} | {r.gene} | {r.prom_magnitude} | {r.prom_summary[:60]} | {ai[:80]} |")
        lines.append("")

    unknowns = [r for r in pydna_only if r.pydna_risk_level == "unknown"]
    if unknowns:
        lines.append("## pydna 'Unknown' Risk Levels\n")
        lines.append("| RSID | Gene | Genotype | AI Suggestion |")
        lines.append("|------|------|----------|---------------|")
        for r in unknowns:
            ai = r.ai_suggestion.replace("\n", " | ") if r.ai_suggestion else "—"
            lines.append(f"| {r.rsid} | {r.gene} | {r.pydna_genotype} | {ai[:100]} |")
        lines.append("")

    output_path.write_text("\n".join(lines))
    print(f"Markdown report exported to: {output_path}")


# ─── CLI ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ai", action="store_true", help="Enable AI enrichment layer")

    # Backend selection — default is claude CLI (subscription), opt-in to SDK (API key)
    backend = parser.add_mutually_exclusive_group()
    backend.add_argument(
        "--claude-cli",
        action="store_true",
        default=True,
        help="Use `claude -p` (Claude subscription, no API key needed) — DEFAULT",
    )
    backend.add_argument(
        "--sdk",
        action="store_true",
        help="Use Anthropic SDK (requires ANTHROPIC_API_KEY env var)",
    )

    parser.add_argument(
        "--model",
        default="claude-sonnet-4-6",
        help="Model for per-variant AI calls (default: claude-sonnet-4-6)",
    )
    parser.add_argument(
        "--narrative-model",
        default="claude-sonnet-4-6",
        help="Model for the executive summary narrative (default: claude-sonnet-4-6)",
    )
    parser.add_argument(
        "--pydna-json",
        type=Path,
        default=None,
        help="Optional path to a precomputed pydna-analyzer analyze.json",
    )
    parser.add_argument(
        "--dna-file",
        type=Path,
        default=DEFAULT_DNA_FILE,
        help="DNA file to analyze live when --pydna-json is not provided",
    )
    parser.add_argument(
        "--promethease-csv",
        type=Path,
        default=PROMETHEASE_CSV,
        help="Path to Promethease Table.csv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output path (.csv or .md). Defaults to comparison_report.md in this directory.",
    )
    args = parser.parse_args()

    # Load data
    print("Loading pydna-analyzer data...")
    if args.pydna_json is not None:
        pydna_variants = load_pydna_from_json(args.pydna_json)
    else:
        pydna_variants = load_pydna_from_dna_file(args.dna_file)
    print(f"  → {len(pydna_variants)} variants")

    print("Loading Promethease data...")
    prom_variants = load_promethease(args.promethease_csv)
    print(f"  → {len(prom_variants)} SNPs")

    # Build comparison
    rows = build_comparison(pydna_variants, prom_variants)

    # Print summary
    print_summary(pydna_variants, prom_variants, rows)

    # AI enrichment
    narrative = ""
    if args.ai:
        use_cli = not args.sdk  # default to claude -p unless --sdk passed
        backend_label = "claude -p (subscription)" if use_cli else "Anthropic SDK"
        print(f"\n[AI enrichment enabled — {backend_label}, model={args.model}]")
        pydna_map = {v.rsid: v for v in pydna_variants}
        rows = enrich_with_ai(rows, pydna_map, model=args.model, use_cli=use_cli)

        print(f"\n[AI narrative — {backend_label}, model={args.narrative_model}]")
        narrative = generate_ai_narrative(
            rows, pydna_variants, prom_variants,
            model=args.narrative_model, use_cli=use_cli,
        )
        print("\n── NARRATIVE SUMMARY ─────────────────────────────────")
        print(narrative)

    # Export
    output = args.output or BASE_DIR / "comparison_report.md"
    suffix = output.suffix.lower()
    if suffix == ".csv":
        export_csv(rows, output)
    else:
        export_markdown(rows, output, narrative=narrative)


if __name__ == "__main__":
    main()
