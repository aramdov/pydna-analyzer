"""
Interactive HTML report generator.
"""

from __future__ import annotations

from pathlib import Path
from datetime import datetime

from jinja2 import Template

from pydna_analyzer.clinical.analyzer import AnalysisResult
from pydna_analyzer.clinical.variants import RiskLevel, Category


HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PyDNA Analyzer Report</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        :root {
            --primary: #3b82f6;
            --danger: #ef4444;
            --warning: #f59e0b;
            --success: #22c55e;
            --bg: #0f172a;
            --card: #1e293b;
            --text: #f8fafc;
            --muted: #94a3b8;
        }
        
        * { box-sizing: border-box; margin: 0; padding: 0; }
        
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: var(--bg);
            color: var(--text);
            line-height: 1.6;
        }
        
        .container { max-width: 1200px; margin: 0 auto; padding: 2rem; }
        
        header {
            text-align: center;
            padding: 3rem 0;
            border-bottom: 1px solid rgba(255,255,255,0.1);
            margin-bottom: 2rem;
        }
        
        header h1 { font-size: 2.5rem; margin-bottom: 0.5rem; }
        header p { color: var(--muted); }
        
        .grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 1.5rem; }
        
        .card {
            background: var(--card);
            border-radius: 12px;
            padding: 1.5rem;
            border: 1px solid rgba(255,255,255,0.05);
        }
        
        .card h3 {
            font-size: 0.875rem;
            color: var(--muted);
            text-transform: uppercase;
            letter-spacing: 0.05em;
            margin-bottom: 0.5rem;
        }
        
        .card .value {
            font-size: 2rem;
            font-weight: 700;
        }
        
        .stat-danger .value { color: var(--danger); }
        .stat-warning .value { color: var(--warning); }
        .stat-success .value { color: var(--success); }
        
        .section { margin: 3rem 0; }
        .section h2 { font-size: 1.5rem; margin-bottom: 1.5rem; display: flex; align-items: center; gap: 0.5rem; }
        
        .apoe-card {
            background: linear-gradient(135deg, var(--card) 0%, rgba(59, 130, 246, 0.1) 100%);
            border-left: 4px solid var(--primary);
        }
        
        .apoe-genotype { font-size: 2.5rem; font-weight: 700; color: var(--primary); }
        .apoe-risk { margin-top: 0.5rem; padding: 0.25rem 0.75rem; border-radius: 999px; display: inline-block; font-size: 0.875rem; }
        .apoe-risk.high { background: rgba(239, 68, 68, 0.2); color: var(--danger); }
        .apoe-risk.elevated { background: rgba(245, 158, 11, 0.2); color: var(--warning); }
        .apoe-risk.average { background: rgba(148, 163, 184, 0.2); color: var(--muted); }
        .apoe-risk.reduced { background: rgba(34, 197, 94, 0.2); color: var(--success); }
        
        .variant-table {
            width: 100%;
            border-collapse: collapse;
            margin: 1rem 0;
        }
        
        .variant-table th, .variant-table td {
            text-align: left;
            padding: 1rem;
            border-bottom: 1px solid rgba(255,255,255,0.05);
        }
        
        .variant-table th { color: var(--muted); font-weight: 500; font-size: 0.875rem; }
        .variant-table tr:hover { background: rgba(255,255,255,0.02); }
        
        .risk-badge {
            padding: 0.25rem 0.75rem;
            border-radius: 999px;
            font-size: 0.75rem;
            font-weight: 600;
        }
        
        .risk-high { background: rgba(239, 68, 68, 0.2); color: var(--danger); }
        .risk-elevated { background: rgba(245, 158, 11, 0.2); color: var(--warning); }
        .risk-moderate { background: rgba(251, 191, 36, 0.2); color: #fbbf24; }
        .risk-low, .risk-normal { background: rgba(34, 197, 94, 0.2); color: var(--success); }
        
        .category-badge {
            padding: 0.25rem 0.5rem;
            border-radius: 6px;
            font-size: 0.75rem;
            background: rgba(255,255,255,0.1);
        }
        
        .interaction-alert {
            background: rgba(245, 158, 11, 0.1);
            border-left: 4px solid var(--warning);
            padding: 1.5rem;
            margin: 1rem 0;
            border-radius: 0 12px 12px 0;
        }
        
        .interaction-alert h4 { color: var(--warning); margin-bottom: 0.5rem; }
        
        .chart-container { height: 300px; margin: 1rem 0; }
        
        .tabs {
            display: flex;
            gap: 0.5rem;
            margin-bottom: 1rem;
            border-bottom: 1px solid rgba(255,255,255,0.1);
            padding-bottom: 0.5rem;
        }
        
        .tab {
            padding: 0.5rem 1rem;
            background: none;
            border: none;
            color: var(--muted);
            cursor: pointer;
            border-radius: 6px;
            transition: all 0.2s;
        }
        
        .tab:hover { color: var(--text); }
        .tab.active { background: var(--primary); color: white; }
        
        .tab-content { display: none; }
        .tab-content.active { display: block; }
        
        .disclaimer {
            background: rgba(239, 68, 68, 0.1);
            border: 1px solid rgba(239, 68, 68, 0.3);
            border-radius: 12px;
            padding: 1.5rem;
            margin-top: 3rem;
        }
        
        .disclaimer h4 { color: var(--danger); margin-bottom: 0.5rem; }
        
        footer {
            text-align: center;
            padding: 2rem;
            color: var(--muted);
            font-size: 0.875rem;
            margin-top: 3rem;
            border-top: 1px solid rgba(255,255,255,0.1);
        }
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>🧬 PyDNA Analyzer Report</h1>
            <p>Generated {{ generated_at }}</p>
        </header>
        
        <!-- Summary Stats -->
        <div class="grid">
            <div class="card">
                <h3>SNPs Analyzed</h3>
                <div class="value">{{ result.snp_count | format_number }}</div>
            </div>
            <div class="card">
                <h3>Clinical Variants Found</h3>
                <div class="value">{{ result.variants_found }}</div>
            </div>
            <div class="card stat-danger">
                <h3>High Priority</h3>
                <div class="value">{{ result.high_priority | length }}</div>
            </div>
            <div class="card stat-warning">
                <h3>Moderate Priority</h3>
                <div class="value">{{ result.moderate_priority | length }}</div>
            </div>
        </div>
        
        {% if result.apoe_status %}
        <!-- APOE Status -->
        <div class="section">
            <h2>🧠 APOE Status</h2>
            <div class="card apoe-card">
                <div class="apoe-genotype">{{ result.apoe_status.genotype }}</div>
                <div class="apoe-risk {{ result.apoe_status.risk_category }}">
                    {{ result.apoe_status.risk_category | upper }} RISK
                </div>
                <p style="margin-top: 1rem; color: var(--muted);">
                    {{ result.apoe_status.interpretation }}
                </p>
                <h4 style="margin-top: 1.5rem; margin-bottom: 0.5rem;">Recommendations:</h4>
                <ul style="padding-left: 1.5rem;">
                    {% for rec in result.apoe_status.recommendations %}
                    <li>{{ rec }}</li>
                    {% endfor %}
                </ul>
            </div>
        </div>
        {% endif %}
        
        {% if result.gene_interactions %}
        <!-- Gene Interactions -->
        <div class="section">
            <h2>⚡ Gene Interactions Detected</h2>
            {% for interaction in result.gene_interactions %}
            <div class="interaction-alert">
                <h4>🔗 {{ interaction.genes | join(' + ') }}</h4>
                <p>{{ interaction.note }}</p>
                {% if interaction.recommendations %}
                <h5 style="margin-top: 1rem;">Recommendations:</h5>
                <ul style="padding-left: 1.5rem;">
                    {% for rec in interaction.recommendations %}
                    <li>{{ rec }}</li>
                    {% endfor %}
                </ul>
                {% endif %}
            </div>
            {% endfor %}
        </div>
        {% endif %}
        
        <!-- Charts -->
        <div class="section">
            <h2>📊 Analysis Overview</h2>
            <div class="grid">
                <div class="card">
                    <div id="category-chart" class="chart-container"></div>
                </div>
                <div class="card">
                    <div id="risk-chart" class="chart-container"></div>
                </div>
            </div>
        </div>
        
        <!-- Variant Tables by Priority -->
        <div class="section">
            <h2>🔬 Detailed Findings</h2>
            
            <div class="tabs">
                <button class="tab active" onclick="showTab('high')">High Priority ({{ result.high_priority | length }})</button>
                <button class="tab" onclick="showTab('moderate')">Moderate ({{ result.moderate_priority | length }})</button>
                <button class="tab" onclick="showTab('all')">All Variants</button>
            </div>
            
            <div id="high" class="tab-content active">
                {% if result.high_priority %}
                <table class="variant-table">
                    <thead>
                        <tr>
                            <th>Gene</th>
                            <th>Variant</th>
                            <th>Your Genotype</th>
                            <th>Risk</th>
                            <th>Interpretation</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for v in result.high_priority %}
                        <tr>
                            <td><strong>{{ v.gene }}</strong><br><span style="color: var(--muted); font-size: 0.75rem;">{{ v.rsid }}</span></td>
                            <td>{{ v.name }}</td>
                            <td style="font-family: monospace;">{{ v.genotype }}</td>
                            <td><span class="risk-badge risk-{{ v.risk_level.value }}">{{ v.risk_level.value | upper }}</span></td>
                            <td>{{ v.description }}</td>
                        </tr>
                        {% endfor %}
                    </tbody>
                </table>
                {% else %}
                <p style="color: var(--muted); padding: 2rem; text-align: center;">No high priority findings</p>
                {% endif %}
            </div>
            
            <div id="moderate" class="tab-content">
                {% if result.moderate_priority %}
                <table class="variant-table">
                    <thead>
                        <tr>
                            <th>Gene</th>
                            <th>Variant</th>
                            <th>Your Genotype</th>
                            <th>Risk</th>
                            <th>Interpretation</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for v in result.moderate_priority %}
                        <tr>
                            <td><strong>{{ v.gene }}</strong><br><span style="color: var(--muted); font-size: 0.75rem;">{{ v.rsid }}</span></td>
                            <td>{{ v.name }}</td>
                            <td style="font-family: monospace;">{{ v.genotype }}</td>
                            <td><span class="risk-badge risk-{{ v.risk_level.value }}">{{ v.risk_level.value | upper }}</span></td>
                            <td>{{ v.description }}</td>
                        </tr>
                        {% endfor %}
                    </tbody>
                </table>
                {% else %}
                <p style="color: var(--muted); padding: 2rem; text-align: center;">No moderate priority findings</p>
                {% endif %}
            </div>
            
            <div id="all" class="tab-content">
                <table class="variant-table">
                    <thead>
                        <tr>
                            <th>Gene</th>
                            <th>Variant</th>
                            <th>Category</th>
                            <th>Your Genotype</th>
                            <th>Risk</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for v in result.variant_results %}
                        <tr>
                            <td><strong>{{ v.gene }}</strong><br><span style="color: var(--muted); font-size: 0.75rem;">{{ v.rsid }}</span></td>
                            <td>{{ v.name }}</td>
                            <td><span class="category-badge">{{ category_icons[v.category.value] }} {{ v.category.value }}</span></td>
                            <td style="font-family: monospace;">{{ v.genotype }}</td>
                            <td><span class="risk-badge risk-{{ v.risk_level.value }}">{{ v.risk_level.value }}</span></td>
                        </tr>
                        {% endfor %}
                    </tbody>
                </table>
            </div>
        </div>
        
        <div class="disclaimer">
            <h4>⚠️ Important Disclaimer</h4>
            <p>This analysis is for <strong>educational purposes only</strong> and is NOT medical advice. 
            Genetic variants interact in complex ways, and environmental factors play major roles. 
            Many variants have incomplete penetrance - having a risk allele does not mean you will 
            develop a condition. Always consult healthcare providers for medical decisions. 
            Consider professional genetic counseling for significant findings.</p>
        </div>
        
        <footer>
            <p>Generated by PyDNA Analyzer v{{ version }}</p>
        </footer>
    </div>
    
    <script>
        // Category distribution chart
        const categoryData = {{ category_data | safe }};
        Plotly.newPlot('category-chart', [{
            values: categoryData.values,
            labels: categoryData.labels,
            type: 'pie',
            hole: 0.4,
            marker: {
                colors: ['#ef4444', '#f59e0b', '#3b82f6', '#22c55e', '#8b5cf6', '#ec4899']
            }
        }], {
            title: { text: 'Variants by Category', font: { color: '#f8fafc' } },
            paper_bgcolor: 'transparent',
            plot_bgcolor: 'transparent',
            font: { color: '#94a3b8' },
            showlegend: true,
            legend: { orientation: 'h', y: -0.1 },
            margin: { t: 40, b: 40, l: 20, r: 20 }
        }, { responsive: true });
        
        // Risk distribution chart
        const riskData = {{ risk_data | safe }};
        Plotly.newPlot('risk-chart', [{
            x: riskData.labels,
            y: riskData.values,
            type: 'bar',
            marker: {
                color: ['#ef4444', '#f59e0b', '#fbbf24', '#22c55e', '#22c55e']
            }
        }], {
            title: { text: 'Variants by Risk Level', font: { color: '#f8fafc' } },
            paper_bgcolor: 'transparent',
            plot_bgcolor: 'transparent',
            font: { color: '#94a3b8' },
            xaxis: { tickangle: -45 },
            margin: { t: 40, b: 80, l: 40, r: 20 }
        }, { responsive: true });
        
        // Tab functionality
        function showTab(tabId) {
            document.querySelectorAll('.tab-content').forEach(el => el.classList.remove('active'));
            document.querySelectorAll('.tab').forEach(el => el.classList.remove('active'));
            document.getElementById(tabId).classList.add('active');
            event.target.classList.add('active');
        }
    </script>
</body>
</html>
"""


CATEGORY_ICONS = {
    "cardiovascular": "❤️",
    "metabolic": "🔥",
    "cancer": "🎗️",
    "pharmacogenomics": "💊",
    "nutrient": "🥗",
    "neuro": "🧠",
}


def generate_html_report(result: AnalysisResult, output_path: Path):
    """Generate an interactive HTML report."""
    from jinja2 import Environment
    from pydna_analyzer import __version__
    
    # Prepare category data for chart
    category_counts = {}
    for v in result.variant_results:
        cat = v.category.value
        category_counts[cat] = category_counts.get(cat, 0) + 1
    
    category_data = {
        "labels": [f"{CATEGORY_ICONS.get(k, '')} {k.title()}" for k in category_counts.keys()],
        "values": list(category_counts.values()),
    }
    
    # Prepare risk data for chart
    risk_counts = {}
    for v in result.variant_results:
        risk = v.risk_level.value
        risk_counts[risk] = risk_counts.get(risk, 0) + 1
    
    risk_order = ["high", "elevated", "moderate", "low", "normal"]
    risk_data = {
        "labels": [r.title() for r in risk_order if r in risk_counts],
        "values": [risk_counts[r] for r in risk_order if r in risk_counts],
    }
    
    # Create Jinja environment with custom filter
    env = Environment()
    env.filters['format_number'] = lambda x: f"{x:,}"
    template = env.from_string(HTML_TEMPLATE)
    
    import json
    html = template.render(
        result=result,
        generated_at=datetime.now().strftime("%Y-%m-%d %H:%M"),
        version=__version__,
        category_data=json.dumps(category_data),
        risk_data=json.dumps(risk_data),
        category_icons=CATEGORY_ICONS,
    )
    
    output_path = Path(output_path)
    output_path.write_text(html)

