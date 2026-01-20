"""
LLM client abstraction for multiple providers.

Supports OpenAI and Anthropic with a unified interface.
"""

from __future__ import annotations

import os
from abc import ABC, abstractmethod
from enum import Enum
from typing import Optional


class LLMProvider(Enum):
    """Supported LLM providers."""
    OPENAI = "openai"
    ANTHROPIC = "anthropic"


class LLMClient(ABC):
    """Abstract base class for LLM clients."""
    
    @abstractmethod
    def generate(
        self,
        prompt: str,
        system: Optional[str] = None,
        max_tokens: int = 4096,
    ) -> str:
        """
        Generate text from a prompt.
        
        Args:
            prompt: The user prompt to send
            system: Optional system prompt
            max_tokens: Maximum tokens in response
            
        Returns:
            Generated text response
        """
        pass
    
    @property
    @abstractmethod
    def provider(self) -> LLMProvider:
        """Return the provider type."""
        pass


class OpenAIClient(LLMClient):
    """OpenAI GPT client."""
    
    def __init__(
        self,
        api_key: Optional[str] = None,
        model: str = "gpt-4o",
    ):
        """
        Initialize OpenAI client.
        
        Args:
            api_key: OpenAI API key (defaults to OPENAI_API_KEY env var)
            model: Model to use (default: gpt-4o)
        """
        self.api_key = api_key or os.environ.get("OPENAI_API_KEY")
        if not self.api_key:
            raise ValueError(
                "OpenAI API key required. Set OPENAI_API_KEY environment variable "
                "or pass api_key parameter."
            )
        self.model = model
        self._client = None
    
    @property
    def provider(self) -> LLMProvider:
        return LLMProvider.OPENAI
    
    def _get_client(self):
        """Lazy-load the OpenAI client."""
        if self._client is None:
            try:
                from openai import OpenAI
                self._client = OpenAI(api_key=self.api_key)
            except ImportError:
                raise ImportError(
                    "OpenAI package not installed. Install with: "
                    "uv add openai  or  pip install openai"
                )
        return self._client
    
    def generate(
        self,
        prompt: str,
        system: Optional[str] = None,
        max_tokens: int = 4096,
    ) -> str:
        """Generate text using OpenAI API."""
        client = self._get_client()
        
        messages = []
        if system:
            messages.append({"role": "system", "content": system})
        messages.append({"role": "user", "content": prompt})
        
        response = client.chat.completions.create(
            model=self.model,
            messages=messages,
            max_tokens=max_tokens,
            temperature=0.7,
        )
        
        return response.choices[0].message.content


class AnthropicClient(LLMClient):
    """Anthropic Claude client."""
    
    def __init__(
        self,
        api_key: Optional[str] = None,
        model: str = "claude-sonnet-4-20250514",
    ):
        """
        Initialize Anthropic client.
        
        Args:
            api_key: Anthropic API key (defaults to ANTHROPIC_API_KEY env var)
            model: Model to use (default: claude-sonnet-4-20250514)
        """
        self.api_key = api_key or os.environ.get("ANTHROPIC_API_KEY")
        if not self.api_key:
            raise ValueError(
                "Anthropic API key required. Set ANTHROPIC_API_KEY environment variable "
                "or pass api_key parameter."
            )
        self.model = model
        self._client = None
    
    @property
    def provider(self) -> LLMProvider:
        return LLMProvider.ANTHROPIC
    
    def _get_client(self):
        """Lazy-load the Anthropic client."""
        if self._client is None:
            try:
                from anthropic import Anthropic
                self._client = Anthropic(api_key=self.api_key)
            except ImportError:
                raise ImportError(
                    "Anthropic package not installed. Install with: "
                    "uv add anthropic  or  pip install anthropic"
                )
        return self._client
    
    def generate(
        self,
        prompt: str,
        system: Optional[str] = None,
        max_tokens: int = 4096,
    ) -> str:
        """Generate text using Anthropic API."""
        client = self._get_client()
        
        kwargs = {
            "model": self.model,
            "max_tokens": max_tokens,
            "messages": [{"role": "user", "content": prompt}],
        }
        
        if system:
            kwargs["system"] = system
        
        response = client.messages.create(**kwargs)
        
        return response.content[0].text


def get_client(
    provider: Optional[LLMProvider] = None,
    api_key: Optional[str] = None,
) -> LLMClient:
    """
    Get an LLM client, auto-detecting provider if not specified.
    
    Args:
        provider: Specific provider to use (auto-detected if None)
        api_key: API key (uses environment variable if None)
        
    Returns:
        Configured LLMClient instance
        
    Raises:
        ValueError: If no API key is available for any provider
    """
    if provider == LLMProvider.OPENAI:
        return OpenAIClient(api_key=api_key)
    elif provider == LLMProvider.ANTHROPIC:
        return AnthropicClient(api_key=api_key)
    
    # Auto-detect from environment
    openai_key = api_key or os.environ.get("OPENAI_API_KEY")
    anthropic_key = api_key or os.environ.get("ANTHROPIC_API_KEY")
    
    if openai_key:
        return OpenAIClient(api_key=openai_key)
    elif anthropic_key:
        return AnthropicClient(api_key=anthropic_key)
    else:
        raise ValueError(
            "No API key found. Set one of:\n"
            "  - OPENAI_API_KEY (for GPT-4)\n"
            "  - ANTHROPIC_API_KEY (for Claude)\n"
            "Or pass --ai-provider and api_key explicitly."
        )
