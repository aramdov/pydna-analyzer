"""
AI-powered report generation module.

Generates natural language reports from genetic analysis results
using OpenAI or Anthropic LLMs.
"""

from pydna_analyzer.ai.client import (
    LLMClient,
    LLMProvider,
    OpenAIClient,
    AnthropicClient,
    get_client,
)
from pydna_analyzer.ai.report_generator import (
    AIReportGenerator,
    ReportStyle,
)

__all__ = [
    "LLMClient",
    "LLMProvider",
    "OpenAIClient",
    "AnthropicClient",
    "get_client",
    "AIReportGenerator",
    "ReportStyle",
]
