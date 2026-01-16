"""
GenomeInsight CLI.

Modern command-line interface for personal genomics analysis.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text
from rich.progress import Progress, SpinnerColumn, TextColumn

from genomeinsight import __version__
from genomeinsight.core.data_loader import load_dna_data, DataFormat
from genomeinsight.clinical.analyzer import ClinicalAnalyzer, AnalysisResult
from genomeinsight.clinical.variants import RiskLevel, Category
from genomeinsight.reports.html_report import generate_html_report
from genomeinsight.reports.json_export import export_to_json

app = typer.Typer(
    name="genomeinsight",
    help="🧬 GenomeInsight: Personal genomics toolkit for DTC DNA analysis",
    add_completion=False,
)
console = Console()

# Risk level styling
RISK_STYLES = {
    RiskLevel.HIGH: ("🔴", "bold red"),
    RiskLevel.ELEVATED: ("🟠", "bold yellow"),
    RiskLevel.MODERATE: ("🟡", "yellow"),
    RiskLevel.LOW: ("🟢", "green"),
    RiskLevel.NORMAL: ("🟢", "green"),
    RiskLevel.UNKNOWN: ("⚪", "dim"),
}

# Category icons
CATEGORY_ICONS = {
    Category.CARDIOVASCULAR: "❤️",
    Category.METABOLIC: "🔥",
    Category.CANCER: "🎗️",
    Category.PHARMACOGENOMICS: "💊",
    Category.NUTRIENT: "🥗",
    Category.NEURO: "🧠",
}


def version_callback(value: bool):
    """Show version and exit."""
    if value:
        console.print(f"[bold blue]GenomeInsight[/] version {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        False,
        "--version",
        "-v",
        callback=version_callback,
        is_eager=True,
        help="Show version",
    ),
):
    """GenomeInsight: Comprehensive personal genomics toolkit."""
    pass


@app.command()
def analyze(
    filepath: Path = typer.Argument(
        ...,
        help="Path to DNA data file (AncestryDNA, 23andMe, MyHeritage, or VCF)",
        exists=True,
    ),
    format: Optional[str] = typer.Option(
        None,
        "--format",
        "-f",
        help="Data format (auto-detected if not specified): ancestrydna, 23andme, myheritage, vcf",
    ),
    output: Optional[Path] = typer.Option(
        None,
        "--output",
        "-o",
        help="Output file path (JSON or HTML based on extension)",
    ),
    html: bool = typer.Option(
        False,
        "--html",
        help="Generate interactive HTML report",
    ),
    json_output: bool = typer.Option(
        False,
        "--json",
        help="Output results as JSON",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet",
        "-q",
        help="Minimal output",
    ),
):
    """
    Analyze DNA data for clinical variants.
    
    Examples:
        genomeinsight analyze AncestryDNA.txt
        genomeinsight analyze 23andme.txt --html -o report.html
        genomeinsight analyze data.vcf --json -o results.json
    """
    # Determine format
    data_format = DataFormat.AUTO
    if format:
        try:
            data_format = DataFormat(format.lower())
        except ValueError:
            console.print(f"[red]Unknown format: {format}[/]")
            raise typer.Exit(1)
    
    # Load data
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
        transient=True,
    ) as progress:
        progress.add_task("Loading DNA data...", total=None)
        try:
            dataset = load_dna_data(filepath, data_format)
        except Exception as e:
            console.print(f"[red]Error loading file: {e}[/]")
            raise typer.Exit(1)
    
    if not quiet:
        console.print(f"[green]✓[/] Loaded {dataset.snp_count:,} SNPs from {filepath.name}")
    
    # Analyze
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
        transient=True,
    ) as progress:
        progress.add_task("Analyzing clinical variants...", total=None)
        analyzer = ClinicalAnalyzer()
        result = analyzer.analyze(dataset)
    
    if not quiet:
        console.print(f"[green]✓[/] Found {result.variants_found} clinical variants")
    
    # Output results
    if json_output or (output and output.suffix == '.json'):
        output_path = output or Path("genomeinsight_results.json")
        export_to_json(result, output_path)
        console.print(f"[green]✓[/] Results saved to {output_path}")
    
    if html or (output and output.suffix == '.html'):
        output_path = output if output and output.suffix == '.html' else Path("genomeinsight_report.html")
        generate_html_report(result, output_path)
        console.print(f"[green]✓[/] HTML report saved to {output_path}")
    
    # Print summary to console
    if not quiet and not json_output:
        _print_summary(result)


def _print_summary(result: AnalysisResult):
    """Print analysis summary to console."""
    console.print()
    
    # Header
    console.print(Panel.fit(
        "[bold blue]🧬 GenomeInsight Analysis Report[/]",
        border_style="blue",
    ))
    console.print()
    
    # APOE Status
    if result.apoe_status:
        apoe = result.apoe_status
        risk_color = {
            "high": "bold red",
            "elevated": "yellow",
            "average": "white",
            "reduced": "green",
        }.get(apoe.risk_category, "white")
        
        console.print(Panel(
            f"[bold]APOE Genotype:[/] {apoe.genotype}\n"
            f"[bold]Risk Category:[/] [{risk_color}]{apoe.risk_category.upper()}[/]\n"
            f"[dim]{apoe.interpretation}[/]",
            title="🧠 APOE Status",
            border_style="blue",
        ))
        console.print()
    
    # Summary table
    table = Table(title="📊 Analysis Summary", show_header=True)
    table.add_column("Metric", style="cyan")
    table.add_column("Count", justify="right")
    table.add_row("Total SNPs analyzed", f"{result.snp_count:,}")
    table.add_row("Clinical variants found", str(result.variants_found))
    table.add_row("High priority findings", f"[red]{len(result.high_priority)}[/]")
    table.add_row("Moderate priority", f"[yellow]{len(result.moderate_priority)}[/]")
    table.add_row("Gene interactions detected", str(len(result.gene_interactions)))
    console.print(table)
    console.print()
    
    # High priority findings
    if result.high_priority:
        console.print("[bold red]🚨 HIGH PRIORITY FINDINGS[/]")
        console.print()
        for r in result.high_priority[:5]:  # Show top 5
            icon, style = RISK_STYLES.get(r.risk_level, ("⚪", "white"))
            console.print(f"  {icon} [bold]{r.gene}[/] - {r.name}")
            console.print(f"     Genotype: [cyan]{r.genotype}[/] | {r.description}")
            if r.interaction_notes:
                console.print(f"     [yellow]⚠ {r.interaction_notes[0]}[/]")
            console.print()
    
    # Gene interactions
    if result.gene_interactions:
        console.print("[bold yellow]⚡ GENE INTERACTIONS DETECTED[/]")
        console.print()
        for interaction in result.gene_interactions:
            console.print(Panel(
                interaction['note'],
                title=f"🔗 {' + '.join(interaction['genes'])}",
                border_style="yellow",
            ))
            console.print()
    
    # Category summary
    console.print("[bold]📁 Findings by Category[/]")
    for category, results in result.by_category.items():
        icon = CATEGORY_ICONS.get(category, "📌")
        high_risk = sum(1 for r in results if r.risk_level in [RiskLevel.HIGH, RiskLevel.ELEVATED])
        console.print(f"  {icon} {category.value.title()}: {len(results)} variants ({high_risk} elevated risk)")
    
    console.print()
    console.print("[dim]Note: This analysis is for educational purposes only. Consult healthcare providers for medical decisions.[/]")


@app.command()
def info(
    filepath: Path = typer.Argument(
        ...,
        help="Path to DNA data file",
        exists=True,
    ),
):
    """
    Show information about a DNA data file.
    """
    try:
        dataset = load_dna_data(filepath)
    except Exception as e:
        console.print(f"[red]Error: {e}[/]")
        raise typer.Exit(1)
    
    console.print(Panel.fit(
        f"[bold]File:[/] {filepath.name}\n"
        f"[bold]Format:[/] {dataset.format.value}\n"
        f"[bold]SNPs:[/] {dataset.snp_count:,}",
        title="📄 File Information",
    ))


@app.command()
def variants():
    """
    List all clinical variants in the database.
    """
    from genomeinsight.clinical.variants import CLINICAL_VARIANTS
    
    table = Table(title="📚 Clinical Variants Database", show_lines=True)
    table.add_column("rsID", style="cyan")
    table.add_column("Gene", style="green")
    table.add_column("Name")
    table.add_column("Category")
    table.add_column("Evidence")
    
    for rsid, variant in CLINICAL_VARIANTS.items():
        icon = CATEGORY_ICONS.get(variant.category, "📌")
        evidence_style = "green" if variant.evidence.value == "strong" else "yellow"
        table.add_row(
            rsid,
            variant.gene,
            variant.name,
            f"{icon} {variant.category.value}",
            f"[{evidence_style}]{variant.evidence.value}[/]",
        )
    
    console.print(table)
    console.print(f"\n[dim]Total: {len(CLINICAL_VARIANTS)} variants[/]")


@app.command()
def prs(
    filepath: Path = typer.Argument(
        ...,
        help="Path to DNA data file",
        exists=True,
    ),
    weights: Path = typer.Option(
        ...,
        "--weights",
        "-w",
        help="Path to weights file (CSV or PGS Catalog format)",
        exists=True,
    ),
    name: Optional[str] = typer.Option(
        None,
        "--name",
        "-n",
        help="Name for this score (defaults to filename)",
    ),
    output: Optional[Path] = typer.Option(
        None,
        "--output",
        "-o",
        help="Output JSON file path",
    ),
):
    """
    Calculate polygenic risk score from weights file.
    
    Examples:
        genomeinsight prs my_dna.txt --weights cad_prs.csv
        genomeinsight prs data.txt -w heart_disease.csv -n "Heart Disease Risk"
    """
    from genomeinsight.polygenic import PRSCalculator
    
    # Load DNA data
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
        transient=True,
    ) as progress:
        progress.add_task("Loading DNA data...", total=None)
        try:
            dataset = load_dna_data(filepath)
        except Exception as e:
            console.print(f"[red]Error loading file: {e}[/]")
            raise typer.Exit(1)
    
    console.print(f"[green]✓[/] Loaded {dataset.snp_count:,} SNPs")
    
    # Calculate PRS
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
        transient=True,
    ) as progress:
        progress.add_task("Calculating polygenic risk score...", total=None)
        try:
            calculator = PRSCalculator()
            result = calculator.calculate_from_file(dataset, weights, score_name=name)
        except Exception as e:
            console.print(f"[red]Error calculating PRS: {e}[/]")
            raise typer.Exit(1)
    
    # Display results
    console.print()
    console.print(Panel.fit(
        f"[bold blue]📊 Polygenic Risk Score: {result.score_name}[/]",
        border_style="blue",
    ))
    console.print()
    
    # Results table
    table = Table(show_header=False, box=None)
    table.add_column("Metric", style="cyan")
    table.add_column("Value", justify="right")
    
    table.add_row("Raw Score", f"{result.raw_score:.4f}")
    if result.normalized_score is not None:
        table.add_row("Z-Score", f"{result.normalized_score:.2f}")
    if result.percentile is not None:
        table.add_row("Percentile", f"{result.percentile:.1f}th")
    table.add_row("Risk Category", f"[bold]{result.risk_category.upper()}[/]")
    table.add_row("SNPs Used", f"{result.snps_used} / {result.snps_available}")
    table.add_row("Coverage", result.coverage_percent)
    
    console.print(table)
    console.print()
    
    # Interpretation
    if result.interpretation:
        console.print(Panel(
            result.interpretation,
            title="📝 Interpretation",
            border_style="dim",
        ))
    
    # Coverage warning
    if result.coverage < 0.5:
        console.print()
        console.print("[yellow]⚠ Low coverage warning: Less than 50% of SNPs were found.[/]")
        console.print("[dim]Results may be less accurate. Consider using a more comprehensive genotyping service.[/]")
    
    # Save to JSON if requested
    if output:
        import json
        output_data = {
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
        with open(output, "w") as f:
            json.dump(output_data, f, indent=2)
        console.print(f"\n[green]✓[/] Results saved to {output}")
    
    console.print()
    console.print("[dim]Note: Polygenic scores are for research/educational purposes only.[/]")


if __name__ == "__main__":
    app()
