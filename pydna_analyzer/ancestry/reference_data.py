"""Reference allele frequency data for ancestry-informative markers."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass
class AIMDatabase:
    """Database of ancestry-informative markers with population frequencies."""

    metadata: dict[str, Any]
    populations: dict[str, dict[str, Any]]
    aims: dict[str, dict[str, Any]]

    @classmethod
    def load(cls, filepath: Path | None = None) -> AIMDatabase:
        """Load AIM database from bundled JSON file."""
        if filepath is None:
            filepath = Path(__file__).parent / "data" / "aim_frequencies.json"
        with filepath.open() as f:
            data = json.load(f)
        return cls(
            metadata=data["metadata"],
            populations=data["populations"],
            aims=data["aims"],
        )

    def get_aim_rsids(self) -> set[str]:
        """Get set of all AIM rsIDs."""
        return set(self.aims.keys())

    def get_population_names(self) -> list[str]:
        """Get list of population names."""
        return list(self.populations.keys())

    def get_frequency(self, rsid: str, population: str) -> float | None:
        """Get allele frequency for a specific AIM and population."""
        aim = self.aims.get(rsid)
        if aim is None:
            return None
        return aim["frequencies"].get(population)

    def get_effect_allele(self, rsid: str) -> str | None:
        """Get the effect allele for an AIM."""
        aim = self.aims.get(rsid)
        if aim is None:
            return None
        return aim["effect_allele"]
