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
    ai_report: bool = typer.Option(
        False,
        "--ai",
        help="Generate AI-powered natural language report",
    ),
    ai_style: str = typer.Option(
        "consumer",
        "--ai-style",
        help="AI report style: 'technical', 'consumer', or 'both'",
    ),
    ai_provider: Optional[str] = typer.Option(
        None,
        "--ai-provider",
        help="LLM provider: 'openai' or 'anthropic' (auto-detected if not set)",
    ),
    ai_output: Optional[Path] = typer.Option(
        None,
        "--ai-output",
        help="Output path for AI report (defaults to console)",
    ),
):
    """
    Analyze DNA data for clinical variants.
    
    Examples:
        genomeinsight analyze AncestryDNA.txt
        genomeinsight analyze 23andme.txt --html -o report.html
        genomeinsight analyze data.vcf --json -o results.json
        genomeinsight analyze data.txt --ai --ai-style consumer
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
    
    # Generate AI report if requested
    if ai_report:
        _generate_ai_report(
            result=result,
            style=ai_style,
            provider=ai_provider,
            output_path=ai_output,
            quiet=quiet,
        )
    
    # Print summary to console
    if not quiet and not json_output and not ai_report:
        _print_summary(result)


def _generate_ai_report(
    result: AnalysisResult,
    style: str,
    provider: Optional[str],
    output_path: Optional[Path],
    quiet: bool,
):
    """Generate AI-powered natural language report."""
    from genomeinsight.ai import AIReportGenerator, ReportStyle, LLMProvider, get_client
    
    # Validate style
    valid_styles = ["technical", "consumer", "both"]
    if style not in valid_styles:
        console.print(f"[red]Invalid AI style: {style}. Choose from: {', '.join(valid_styles)}[/]")
        raise typer.Exit(1)
    
    # Get LLM client
    try:
        llm_provider = None
        if provider:
            try:
                llm_provider = LLMProvider(provider.lower())
            except ValueError:
                console.print(f"[red]Unknown provider: {provider}. Choose 'openai' or 'anthropic'[/]")
                raise typer.Exit(1)
        
        client = get_client(provider=llm_provider)
        
        if not quiet:
            console.print(f"[green]✓[/] Using {client.provider.value.title()} for AI report")
    
    except ValueError as e:
        console.print(f"[red]AI Error: {e}[/]")
        console.print("[dim]Set OPENAI_API_KEY or ANTHROPIC_API_KEY environment variable[/]")
        raise typer.Exit(1)
    except ImportError as e:
        console.print(f"[red]Missing dependency: {e}[/]")
        console.print("[dim]Install with: uv add openai  or  uv add anthropic[/]")
        raise typer.Exit(1)
    
    generator = AIReportGenerator(client)
    
    # Generate report(s)
    if style == "both":
        styles_to_generate = [ReportStyle.TECHNICAL, ReportStyle.CONSUMER]
    else:
        styles_to_generate = [ReportStyle.TECHNICAL if style == "technical" else ReportStyle.CONSUMER]
    
    for report_style in styles_to_generate:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
            transient=True,
        ) as progress:
            progress.add_task(f"Generating {report_style.value} AI report...", total=None)
            try:
                report = generator.generate(result, style=report_style)
            except Exception as e:
                console.print(f"[red]Error generating AI report: {e}[/]")
                raise typer.Exit(1)
        
        if output_path:
            # Append style suffix if generating both
            if style == "both":
                stem = output_path.stem
                suffix = output_path.suffix or ".md"
                final_path = output_path.parent / f"{stem}_{report_style.value}{suffix}"
            else:
                final_path = output_path
            
            with open(final_path, "w") as f:
                f.write(report)
            console.print(f"[green]✓[/] {report_style.value.title()} AI report saved to {final_path}")
        else:
            # Print to console
            console.print()
            console.print(Panel.fit(
                f"[bold blue]🤖 AI-Generated {report_style.value.title()} Report[/]",
                border_style="blue",
            ))
            console.print()
            from rich.markdown import Markdown
            console.print(Markdown(report))


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


@app.command()
def pgx(
    filepath: Path = typer.Argument(
        ...,
        help="Path to DNA data file (AncestryDNA, 23andMe, MyHeritage, or VCF)",
        exists=True,
    ),
    gene: Optional[str] = typer.Option(
        None,
        "--gene",
        "-g",
        help="Analyze a single gene (e.g. CYP2C19)",
    ),
    output: Optional[Path] = typer.Option(
        None,
        "--output",
        "-o",
        help="Output JSON file path",
    ),
):
    """
    Pharmacogenomics analysis: star alleles, metabolizer phenotypes, drug recommendations.

    Examples:
        genomeinsight pgx AncestryDNA.txt
        genomeinsight pgx data.txt --gene CYP2C19
        genomeinsight pgx data.txt -o pgx_results.json
    """
    from genomeinsight.pharmacogenomics import PGxAnalyzer, MetabolizerPhenotype

    # Load data
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

    console.print(f"[green]✓[/] Loaded {dataset.snp_count:,} SNPs from {filepath.name}")

    # Run PGx analysis
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
        transient=True,
    ) as progress:
        progress.add_task("Running pharmacogenomics analysis...", total=None)
        analyzer = PGxAnalyzer()
        if gene:
            single_result = analyzer.analyze_gene(dataset, gene.upper())
            if single_result is None:
                console.print(f"[red]Unknown gene: {gene}[/]")
                console.print(
                    f"[dim]Available genes: {', '.join(sorted(analyzer.gene_definitions))}[/]"
                )
                raise typer.Exit(1)
            # Wrap in PGxResult for consistent display
            from genomeinsight.pharmacogenomics import PGxResult

            result = PGxResult(
                gene_results=[single_result],
                total_genes_tested=1,
                total_snps_tested=single_result.snps_tested,
                total_snps_missing=single_result.snps_missing,
            )
        else:
            result = analyzer.analyze(dataset)

    _print_pgx_results(result)

    # Save JSON if requested
    if output:
        import json

        output_data = {
            "genes": [
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
        with output.open("w") as f:
            json.dump(output_data, f, indent=2)
        console.print(f"\n[green]✓[/] Results saved to {output}")


def _print_pgx_results(result):
    """Print pharmacogenomics results with Rich tables."""
    from genomeinsight.pharmacogenomics import MetabolizerPhenotype

    console.print()
    console.print(
        Panel.fit(
            "[bold blue]💊 Pharmacogenomics Report[/]",
            border_style="blue",
        )
    )
    console.print()

    # Gene summary table
    gene_table = Table(title="Gene Summary", show_lines=True)
    gene_table.add_column("Gene", style="cyan")
    gene_table.add_column("Diplotype", style="green")
    gene_table.add_column("Phenotype")
    gene_table.add_column("Activity Score", justify="right")
    gene_table.add_column("Confidence")

    PHENOTYPE_STYLES = {
        MetabolizerPhenotype.POOR: ("🔴", "bold red"),
        MetabolizerPhenotype.INTERMEDIATE: ("🟠", "yellow"),
        MetabolizerPhenotype.NORMAL: ("🟢", "green"),
        MetabolizerPhenotype.RAPID: ("🔵", "blue"),
        MetabolizerPhenotype.ULTRARAPID: ("⚡", "bold blue"),
        MetabolizerPhenotype.HIGH_SENSITIVITY: ("🔴", "bold red"),
        MetabolizerPhenotype.INTERMEDIATE_SENSITIVITY: ("🟠", "yellow"),
        MetabolizerPhenotype.NORMAL_SENSITIVITY: ("🟢", "green"),
        MetabolizerPhenotype.HIGH_ACTIVITY: ("🟢", "green"),
        MetabolizerPhenotype.INTERMEDIATE_ACTIVITY: ("🟠", "yellow"),
        MetabolizerPhenotype.LOW_ACTIVITY: ("🔴", "bold red"),
    }

    for gr in result.gene_results:
        icon, style = PHENOTYPE_STYLES.get(gr.phenotype, ("⚪", "white"))
        score_str = f"{gr.activity_score:.1f}" if gr.activity_score is not None else "—"
        conf_style = {"high": "green", "moderate": "yellow", "low": "red"}.get(
            gr.confidence, "white"
        )
        gene_table.add_row(
            gr.gene,
            gr.diplotype,
            f"{icon} [{style}]{gr.phenotype.value}[/]",
            score_str,
            f"[{conf_style}]{gr.confidence}[/]",
        )

    console.print(gene_table)
    console.print()

    # Drug recommendations
    all_recs = result.all_drug_recommendations
    if all_recs:
        rec_table = Table(title="💊 Drug Recommendations", show_lines=True)
        rec_table.add_column("Drug", style="bold")
        rec_table.add_column("Gene", style="cyan")
        rec_table.add_column("CPIC", justify="center")
        rec_table.add_column("Recommendation")
        rec_table.add_column("Alternatives", style="dim")

        for rec in all_recs:
            cpic_style = {
                "strong": "bold red",
                "moderate": "yellow",
                "optional": "dim",
            }.get(rec.cpic_level.value, "white")
            rec_table.add_row(
                rec.drug_name,
                rec.gene,
                f"[{cpic_style}]{rec.cpic_level.value.upper()}[/]",
                rec.recommendation,
                ", ".join(rec.alternatives) if rec.alternatives else "—",
            )

        console.print(rec_table)
    else:
        console.print("[green]No actionable drug-gene interactions detected.[/]")

    console.print()
    console.print(f"[dim]{result.summary}[/]")
    console.print()
    console.print(
        "[dim]⚠ Disclaimer: This analysis is for educational purposes only. "
        "Do not adjust medications without consulting a healthcare provider and "
        "clinical-grade pharmacogenomic testing.[/]"
    )


@app.command()
def ancestry(
    filepath: Path = typer.Argument(
        ...,
        help="Path to DNA data file (AncestryDNA, 23andMe, MyHeritage, or VCF)",
        exists=True,
    ),
    output: Optional[Path] = typer.Option(
        None,
        "--output",
        "-o",
        help="Output JSON file path",
    ),
    bootstrap: int = typer.Option(
        100,
        "--bootstrap",
        "-b",
        help="Number of bootstrap iterations for confidence intervals",
    ),
):
    """
    Estimate ancestry composition from DNA data.

    Uses maximum likelihood estimation with ancestry-informative markers
    to estimate population proportions with confidence intervals.

    Examples:
        genomeinsight ancestry AncestryDNA.txt
        genomeinsight ancestry data.txt -o ancestry.json
        genomeinsight ancestry data.txt --bootstrap 200
    """
    from genomeinsight.ancestry import AncestryAnalyzer

    # Load data
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

    console.print(f"[green]✓[/] Loaded {dataset.snp_count:,} SNPs from {filepath.name}")

    # Run ancestry estimation
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
        transient=True,
    ) as progress:
        progress.add_task("Estimating ancestry composition...", total=None)
        analyzer = AncestryAnalyzer(n_bootstrap=bootstrap)
        result = analyzer.analyze(dataset)

    coverage_pct = f"{result.coverage * 100:.1f}%"
    console.print(
        f"[green]✓[/] Matched {result.snps_used} / {result.snps_available} "
        f"ancestry markers ({coverage_pct})"
    )

    if result.snps_used == 0:
        console.print(f"\n[yellow]{result.interpretation}[/]")
        raise typer.Exit(0)

    _print_ancestry_results(result)

    # Save JSON if requested
    if output:
        import json

        output_data = {
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
            "interpretation": result.interpretation,
        }
        with output.open("w") as f:
            json.dump(output_data, f, indent=2)
        console.print(f"\n[green]✓[/] Results saved to {output}")


def _print_ancestry_results(result: "AncestryResult") -> None:  # noqa: F821
    """Print ancestry estimation results with Rich tables."""
    console.print()
    console.print(
        Panel.fit(
            "[bold blue]🌍 Ancestry Composition[/]",
            border_style="blue",
        )
    )
    console.print()

    # Population table with visual bars
    table = Table(show_lines=True)
    table.add_column("Population", style="bold")
    table.add_column("%", justify="right", style="cyan")
    table.add_column("", min_width=25)
    table.add_column("95% CI", style="dim")

    for pop in result.populations:
        if pop.proportion < 0.02:
            continue
        pct = f"{pop.proportion * 100:.1f}%"
        bar_len = int(pop.proportion * 30)
        bar = "█" * bar_len
        ci = f"{pop.confidence_low * 100:.1f} - {pop.confidence_high * 100:.1f}%"
        table.add_row(pop.population, pct, f"[blue]{bar}[/]", ci)

    console.print(table)
    console.print()

    # Regional summary
    console.print("[bold]📊 Regional Summary[/]")
    for region, proportion in sorted(
        result.top_regions.items(), key=lambda x: x[1], reverse=True
    ):
        if proportion >= 0.02:
            console.print(f"  {region}: {proportion * 100:.1f}%")

    # Interpretation
    console.print()
    console.print(
        Panel(
            result.interpretation,
            title="📝 Interpretation",
            border_style="dim",
        )
    )

    console.print()
    console.print(
        "[dim]⚠ Note: Ancestry estimates are based on a curated set of "
        "ancestry-informative markers. Commercial services use significantly "
        "more data for higher resolution. Results are for educational purposes only.[/]"
    )


@app.command()
def serve(
    host: str = typer.Option("127.0.0.1", "--host", help="Host to bind to"),
    port: int = typer.Option(8000, "--port", "-p", help="Port to bind to"),
    reload: bool = typer.Option(False, "--reload", help="Enable auto-reload for development"),
):
    """
    Start the GenomeInsight REST API server.

    Examples:
        genomeinsight serve
        genomeinsight serve --port 9000
        genomeinsight serve --reload
    """
    try:
        import uvicorn
    except ImportError:
        console.print("[red]Missing dependency: uvicorn[/]")
        console.print("[dim]Install with: uv sync --extra api[/]")
        raise typer.Exit(1)

    try:
        from genomeinsight.api import create_app  # noqa: F401
    except ImportError:
        console.print("[red]Missing dependency: fastapi[/]")
        console.print("[dim]Install with: uv sync --extra api[/]")
        raise typer.Exit(1)

    console.print(f"[bold blue]🧬 GenomeInsight API[/] starting on http://{host}:{port}")
    console.print("[dim]API docs: http://{}:{}/docs[/]".format(host, port))

    uvicorn.run(
        "genomeinsight.api:create_app",
        host=host,
        port=port,
        reload=reload,
        factory=True,
    )


if __name__ == "__main__":
    app()
