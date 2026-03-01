"""
Tests for AI report generation module.

Uses mocked API responses to test without real API keys.
"""

import pytest
from unittest.mock import Mock, patch, MagicMock
from dataclasses import dataclass, field
from typing import Optional
from enum import Enum


# Create mock types for testing without importing real modules
class MockRiskLevel(Enum):
    HIGH = "high"
    ELEVATED = "elevated"
    MODERATE = "moderate"
    LOW = "low"
    NORMAL = "normal"
    UNKNOWN = "unknown"


class MockCategory(Enum):
    CARDIOVASCULAR = "cardiovascular"
    METABOLIC = "metabolic"
    CANCER = "cancer"
    PHARMACOGENOMICS = "pharmacogenomics"
    NUTRIENT = "nutrient"
    NEURO = "neuro"


class MockEvidenceLevel(Enum):
    STRONG = "strong"
    MODERATE = "moderate"
    LIMITED = "limited"


@dataclass
class MockVariantResult:
    rsid: str = "rs1234"
    gene: str = "TEST"
    name: str = "Test Variant"
    category: MockCategory = MockCategory.CARDIOVASCULAR
    genotype: str = "AG"
    risk_level: MockRiskLevel = MockRiskLevel.MODERATE
    description: str = "Test description"
    evidence: MockEvidenceLevel = MockEvidenceLevel.STRONG
    interventions: list = field(default_factory=list)


@dataclass
class MockAPOEResult:
    genotype: str = "ε3/ε3"
    risk_category: str = "average"
    has_e4: bool = False
    has_e2: bool = False
    is_e4_homozygous: bool = False
    interpretation: str = "Average risk"
    recommendations: list = field(default_factory=lambda: ["Lifestyle measures"])


@dataclass
class MockAnalysisResult:
    source_file: str = "test.txt"
    snp_count: int = 100000
    variants_analyzed: int = 35
    variants_found: int = 15
    apoe_status: Optional[MockAPOEResult] = None
    high_priority: list = field(default_factory=list)
    moderate_priority: list = field(default_factory=list)
    low_priority: list = field(default_factory=list)
    gene_interactions: list = field(default_factory=list)
    by_category: dict = field(default_factory=dict)


class TestLLMClients:
    """Tests for LLM client classes."""
    
    def test_openai_client_requires_api_key(self):
        """Test OpenAI client raises error without API key."""
        from pydna_analyzer.ai.client import OpenAIClient
        
        with patch.dict("os.environ", {}, clear=True):
            with pytest.raises(ValueError, match="OpenAI API key required"):
                OpenAIClient()
    
    def test_anthropic_client_requires_api_key(self):
        """Test Anthropic client raises error without API key."""
        from pydna_analyzer.ai.client import AnthropicClient
        
        with patch.dict("os.environ", {}, clear=True):
            with pytest.raises(ValueError, match="Anthropic API key required"):
                AnthropicClient()
    
    def test_openai_client_accepts_api_key_param(self):
        """Test OpenAI client accepts API key as parameter."""
        from pydna_analyzer.ai.client import OpenAIClient, LLMProvider
        
        client = OpenAIClient(api_key="test-key-123")
        assert client.api_key == "test-key-123"
        assert client.provider == LLMProvider.OPENAI
    
    def test_anthropic_client_accepts_api_key_param(self):
        """Test Anthropic client accepts API key as parameter."""
        from pydna_analyzer.ai.client import AnthropicClient, LLMProvider
        
        client = AnthropicClient(api_key="test-key-123")
        assert client.api_key == "test-key-123"
        assert client.provider == LLMProvider.ANTHROPIC
    
    def test_openai_client_from_env(self):
        """Test OpenAI client reads from environment."""
        from pydna_analyzer.ai.client import OpenAIClient
        
        with patch.dict("os.environ", {"OPENAI_API_KEY": "env-key-456"}):
            client = OpenAIClient()
            assert client.api_key == "env-key-456"
    
    def test_anthropic_client_from_env(self):
        """Test Anthropic client reads from environment."""
        from pydna_analyzer.ai.client import AnthropicClient
        
        with patch.dict("os.environ", {"ANTHROPIC_API_KEY": "env-key-789"}):
            client = AnthropicClient()
            assert client.api_key == "env-key-789"


class TestGetClient:
    """Tests for get_client auto-detection."""
    
    def test_get_client_no_keys_raises_error(self):
        """Test get_client raises error when no keys available."""
        from pydna_analyzer.ai.client import get_client
        
        with patch.dict("os.environ", {}, clear=True):
            with pytest.raises(ValueError, match="No API key found"):
                get_client()
    
    def test_get_client_prefers_openai(self):
        """Test get_client prefers OpenAI when both keys present."""
        from pydna_analyzer.ai.client import get_client, LLMProvider
        
        with patch.dict("os.environ", {
            "OPENAI_API_KEY": "openai-key",
            "ANTHROPIC_API_KEY": "anthropic-key",
        }):
            client = get_client()
            assert client.provider == LLMProvider.OPENAI
    
    def test_get_client_uses_anthropic_if_only_anthropic(self):
        """Test get_client uses Anthropic when only that key exists."""
        from pydna_analyzer.ai.client import get_client, LLMProvider
        
        with patch.dict("os.environ", {"ANTHROPIC_API_KEY": "anthropic-key"}, clear=True):
            client = get_client()
            assert client.provider == LLMProvider.ANTHROPIC
    
    def test_get_client_explicit_provider(self):
        """Test get_client respects explicit provider."""
        from pydna_analyzer.ai.client import get_client, LLMProvider
        
        with patch.dict("os.environ", {"ANTHROPIC_API_KEY": "anthropic-key"}):
            client = get_client(provider=LLMProvider.ANTHROPIC)
            assert client.provider == LLMProvider.ANTHROPIC


class TestPrompts:
    """Tests for prompt building."""
    
    def test_build_prompt_consumer_style(self):
        """Test consumer style prompt generation."""
        from pydna_analyzer.ai.prompts import build_prompt
        
        result = MockAnalysisResult()
        prompt = build_prompt(result, style="consumer")
        
        assert "100,000" in prompt  # SNP count formatted
        assert "plain-language" in prompt.lower() or "consumer" in prompt.lower()
    
    def test_build_prompt_technical_style(self):
        """Test technical style prompt generation."""
        from pydna_analyzer.ai.prompts import build_prompt
        
        result = MockAnalysisResult()
        prompt = build_prompt(result, style="technical")
        
        assert "100,000" in prompt
        assert "technical" in prompt.lower() or "researcher" in prompt.lower()
    
    def test_format_apoe_section_with_status(self):
        """Test APOE section formatting with status."""
        from pydna_analyzer.ai.prompts import format_apoe_section
        
        apoe = MockAPOEResult()
        section = format_apoe_section(apoe)
        
        assert "ε3/ε3" in section
        assert "average" in section.lower()
    
    def test_format_apoe_section_without_status(self):
        """Test APOE section formatting without status."""
        from pydna_analyzer.ai.prompts import format_apoe_section
        
        section = format_apoe_section(None)
        assert "not determined" in section.lower()
    
    def test_format_variant_section_empty(self):
        """Test empty variant section."""
        from pydna_analyzer.ai.prompts import format_variant_section
        
        section = format_variant_section([])
        assert "none" in section.lower()
    
    def test_format_variant_section_with_variants(self):
        """Test variant section with data."""
        from pydna_analyzer.ai.prompts import format_variant_section
        
        variants = [MockVariantResult()]
        section = format_variant_section(variants)
        
        assert "TEST" in section
        assert "rs1234" in section
        assert "AG" in section


class TestReportGenerator:
    """Tests for AIReportGenerator."""
    
    def test_generator_init_with_client(self):
        """Test generator initialization with client."""
        from pydna_analyzer.ai.report_generator import AIReportGenerator
        
        mock_client = Mock()
        generator = AIReportGenerator(client=mock_client)
        assert generator._client == mock_client
    
    def test_generator_generate_consumer_report(self):
        """Test consumer report generation."""
        from pydna_analyzer.ai.report_generator import AIReportGenerator, ReportStyle
        from pydna_analyzer.ai.client import LLMProvider
        
        mock_client = Mock()
        mock_client.provider = LLMProvider.OPENAI
        mock_client.generate.return_value = "# AI Generated Report\n\nThis is a test report."
        
        generator = AIReportGenerator(client=mock_client)
        result = MockAnalysisResult()
        
        report = generator.generate(result, style=ReportStyle.CONSUMER)
        
        assert "PyDNA Analyzer" in report
        assert "Disclaimer" in report
        mock_client.generate.assert_called_once()
    
    def test_generator_generate_technical_report(self):
        """Test technical report generation."""
        from pydna_analyzer.ai.report_generator import AIReportGenerator, ReportStyle
        from pydna_analyzer.ai.client import LLMProvider
        
        mock_client = Mock()
        mock_client.provider = LLMProvider.ANTHROPIC
        mock_client.generate.return_value = "## Technical Analysis\n\nDetailed findings..."
        
        generator = AIReportGenerator(client=mock_client)
        result = MockAnalysisResult()
        
        report = generator.generate(result, style=ReportStyle.TECHNICAL)
        
        assert "Technical" in report
        assert "Disclaimer" in report
        mock_client.generate.assert_called_once()
    
    def test_generator_generate_both(self):
        """Test generating both report styles."""
        from pydna_analyzer.ai.report_generator import AIReportGenerator
        from pydna_analyzer.ai.client import LLMProvider
        
        mock_client = Mock()
        mock_client.provider = LLMProvider.OPENAI
        mock_client.generate.return_value = "Test report content"
        
        generator = AIReportGenerator(client=mock_client)
        result = MockAnalysisResult()
        
        reports = generator.generate_both(result)
        
        assert "technical" in reports
        assert "consumer" in reports
        assert mock_client.generate.call_count == 2


class TestOpenAIClientGeneration:
    """Tests for OpenAI client generation with mocked SDK."""
    
    def test_openai_generate_with_system_prompt(self):
        """Test OpenAI generation with system prompt."""
        from pydna_analyzer.ai.client import OpenAIClient
        
        # Create mock OpenAI SDK
        mock_openai_module = MagicMock()
        mock_client_instance = MagicMock()
        mock_response = MagicMock()
        mock_response.choices = [MagicMock()]
        mock_response.choices[0].message.content = "Generated response"
        mock_client_instance.chat.completions.create.return_value = mock_response
        mock_openai_module.OpenAI.return_value = mock_client_instance
        
        with patch.dict("sys.modules", {"openai": mock_openai_module}):
            client = OpenAIClient(api_key="test-key")
            result = client.generate(
                prompt="Test prompt",
                system="Test system prompt"
            )
        
        assert result == "Generated response"
        mock_client_instance.chat.completions.create.assert_called_once()
        call_args = mock_client_instance.chat.completions.create.call_args
        messages = call_args.kwargs["messages"]
        assert len(messages) == 2
        assert messages[0]["role"] == "system"
        assert messages[1]["role"] == "user"


class TestAnthropicClientGeneration:
    """Tests for Anthropic client generation with mocked SDK."""
    
    def test_anthropic_generate_with_system_prompt(self):
        """Test Anthropic generation with system prompt."""
        from pydna_analyzer.ai.client import AnthropicClient
        
        # Create mock Anthropic SDK
        mock_anthropic_module = MagicMock()
        mock_client_instance = MagicMock()
        mock_response = MagicMock()
        mock_content_block = MagicMock()
        mock_content_block.text = "Generated response"
        mock_response.content = [mock_content_block]
        mock_client_instance.messages.create.return_value = mock_response
        mock_anthropic_module.Anthropic.return_value = mock_client_instance
        
        with patch.dict("sys.modules", {"anthropic": mock_anthropic_module}):
            client = AnthropicClient(api_key="test-key")
            result = client.generate(
                prompt="Test prompt",
                system="Test system prompt"
            )
        
        assert result == "Generated response"
        mock_client_instance.messages.create.assert_called_once()
        call_args = mock_client_instance.messages.create.call_args
        assert call_args.kwargs["system"] == "Test system prompt"
