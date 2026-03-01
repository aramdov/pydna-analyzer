"""
Multi-format DNA data loaders.

Supports:
- AncestryDNA (.txt)
- 23andMe (.txt)
- MyHeritage (.csv)
- VCF files (.vcf, .vcf.gz)
"""

from __future__ import annotations

import gzip
from abc import ABC, abstractmethod
from enum import Enum
from pathlib import Path
from typing import Optional

import pandas as pd
from pydantic import BaseModel, ConfigDict


class DataFormat(str, Enum):
    """Supported DNA data file formats."""
    ANCESTRY_DNA = "ancestrydna"
    TWENTY_THREE_ME = "23andme"
    MY_HERITAGE = "myheritage"
    VCF = "vcf"
    AUTO = "auto"


class GenotypeData(BaseModel):
    """Container for loaded genotype data."""
    
    rsid: str
    chromosome: str
    position: int
    allele1: str
    allele2: str
    genotype: str
    
    model_config = ConfigDict(frozen=True)


class DNADataset(BaseModel):
    """Container for a complete DNA dataset."""
    
    format: DataFormat
    source_file: str
    snp_count: int
    dataframe: object  # pd.DataFrame - using object for Pydantic compatibility
    
    model_config = ConfigDict(arbitrary_types_allowed=True)
    
    def get_genotype(self, rsid: str) -> Optional[str]:
        """Get genotype for a specific rsID."""
        df = self.dataframe
        match = df[df['rsid'] == rsid]
        if len(match) > 0:
            return match.iloc[0]['genotype']
        return None
    
    def get_rsids(self) -> set[str]:
        """Get set of all rsIDs in the dataset."""
        return set(self.dataframe['rsid'].values)
    
    def filter_by_rsids(self, rsids: list[str]) -> pd.DataFrame:
        """Filter dataset to specific rsIDs."""
        return self.dataframe[self.dataframe['rsid'].isin(rsids)]


class BaseLoader(ABC):
    """Abstract base class for DNA data loaders."""
    
    @abstractmethod
    def load(self, filepath: Path) -> pd.DataFrame:
        """Load DNA data from file."""
        pass
    
    @abstractmethod
    def detect(self, filepath: Path) -> bool:
        """Check if this loader can handle the file."""
        pass
    
    def _normalize_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """Normalize column names and create genotype column."""
        # Ensure standard column names
        df.columns = df.columns.str.lower().str.strip()
        
        # Rename common variations
        column_mapping = {
            'snpid': 'rsid',
            'snp': 'rsid',
            'marker': 'rsid',
            'chr': 'chromosome',
            'chrom': 'chromosome',
            'pos': 'position',
            'alleles': 'genotype',
        }
        df = df.rename(columns=column_mapping)
        
        # Create genotype if separate allele columns exist
        if 'genotype' not in df.columns:
            if 'allele1' in df.columns and 'allele2' in df.columns:
                df['genotype'] = df['allele1'].astype(str) + df['allele2'].astype(str)
        
        return df


class AncestryDNALoader(BaseLoader):
    """Loader for AncestryDNA raw data files."""
    
    def detect(self, filepath: Path) -> bool:
        """Detect AncestryDNA format by checking file structure."""
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        if 'AncestryDNA' in line:
                            return True
                    else:
                        # Check column structure
                        parts = line.strip().split('\t')
                        if len(parts) == 5:  # rsid, chr, pos, allele1, allele2
                            return True
                        break
        except Exception:
            pass
        return False
    
    def load(self, filepath: Path) -> pd.DataFrame:
        """Load AncestryDNA format file."""
        # Read file, skipping comment lines and the header row
        # AncestryDNA files have: comments (#), then a header (rsid chromosome...), then data
        lines = []
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                # Skip header row
                if line.strip().startswith('rsid'):
                    continue
                lines.append(line)
        
        # Parse the data
        import io
        df = pd.read_csv(
            io.StringIO(''.join(lines)),
            sep='\t',
            names=['rsid', 'chromosome', 'position', 'allele1', 'allele2'],
            dtype={
                'rsid': str,
                'chromosome': str,
                'position': int,
                'allele1': str,
                'allele2': str,
            }
        )
        df['genotype'] = df['allele1'] + df['allele2']
        return df


class TwentyThreeMeLoader(BaseLoader):
    """Loader for 23andMe raw data files."""
    
    def detect(self, filepath: Path) -> bool:
        """Detect 23andMe format."""
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        if '23andMe' in line:
                            return True
                    else:
                        # 23andMe uses: rsid, chr, pos, genotype (single column)
                        parts = line.strip().split('\t')
                        if len(parts) == 4:
                            return True
                        break
        except Exception:
            pass
        return False
    
    def load(self, filepath: Path) -> pd.DataFrame:
        """Load 23andMe format file."""
        df = pd.read_csv(
            filepath,
            sep='\t',
            comment='#',
            names=['rsid', 'chromosome', 'position', 'genotype'],
            dtype={
                'rsid': str,
                'chromosome': str,
                'position': int,
                'genotype': str,
            }
        )
        # Split genotype into alleles
        df['allele1'] = df['genotype'].str[0]
        df['allele2'] = df['genotype'].str[1:].str.strip()
        df.loc[df['allele2'] == '', 'allele2'] = df.loc[df['allele2'] == '', 'allele1']
        return df


class MyHeritageLoader(BaseLoader):
    """Loader for MyHeritage raw data files."""
    
    def detect(self, filepath: Path) -> bool:
        """Detect MyHeritage format."""
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    if 'MyHeritage' in line:
                        return True
                    if not line.startswith('#'):
                        # Check for CSV format with RSID header
                        if 'RSID' in line.upper():
                            return True
                        break
        except Exception:
            pass
        return False
    
    def load(self, filepath: Path) -> pd.DataFrame:
        """Load MyHeritage format file."""
        df = pd.read_csv(
            filepath,
            comment='#',
        )
        df = self._normalize_dataframe(df)
        return df


class VCFLoader(BaseLoader):
    """Loader for VCF (Variant Call Format) files."""

    def __init__(
        self,
        sample: str | int | None = None,
        include_filtered: bool = False,
    ) -> None:
        self.sample = sample
        self.include_filtered = include_filtered
        self.warnings: list[str] = []

    def detect(self, filepath: Path) -> bool:
        """Detect VCF format."""
        try:
            opener = gzip.open if str(filepath).endswith(".gz") else open
            with opener(filepath, "rt") as f:
                for line in f:
                    if line.startswith("##fileformat=VCF"):
                        return True
                    if not line.startswith("#"):
                        break
        except Exception:
            pass
        return False

    def load(self, filepath: Path) -> pd.DataFrame:
        """Load VCF format file.

        Raises:
            ValueError: If the file is missing a ##fileformat=VCF header or #CHROM line.
        """
        self.warnings = []
        opener = gzip.open if str(filepath).endswith(".gz") else open

        has_fileformat = False
        header_cols: list[str] = []
        sample_idx: int | None = None
        rows: list[dict[str, object]] = []

        with opener(filepath, "rt") as f:
            for line in f:
                if line.startswith("##"):
                    if line.startswith("##fileformat=VCF"):
                        has_fileformat = True
                    continue
                if line.startswith("#CHROM"):
                    header_cols = line.strip().split("\t")
                    continue

                # Data row — but first check we saw both required headers
                if not has_fileformat:
                    raise ValueError(
                        f"VCF file missing ##fileformat=VCF header: {filepath}"
                    )
                if not header_cols:
                    raise ValueError(
                        f"VCF file missing #CHROM header line: {filepath}"
                    )

                # Resolve sample index once on the first data row
                if sample_idx is None:
                    sample_idx = self._resolve_sample_index(header_cols)

                parsed = self._parse_data_row(line, sample_idx)
                if parsed is not None:
                    rows.append(parsed)

        # Files with headers but no data rows still need validation
        if not has_fileformat:
            raise ValueError(
                f"VCF file missing ##fileformat=VCF header: {filepath}"
            )
        if not header_cols:
            raise ValueError(
                f"VCF file missing #CHROM header line: {filepath}"
            )

        return pd.DataFrame(rows)

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _resolve_sample_index(self, header_cols: list[str]) -> int:
        """Return the column index for the chosen sample.

        The first 9 VCF columns are fixed (CHROM..FORMAT).
        Sample columns start at index 9.
        """
        sample_names = header_cols[9:]
        if not sample_names:
            raise ValueError("VCF file has no sample columns")

        if self.sample is None:
            return 9  # first sample

        if isinstance(self.sample, int):
            if self.sample < 0 or self.sample >= len(sample_names):
                raise ValueError(
                    f"Sample index {self.sample} out of range "
                    f"(file has {len(sample_names)} samples)"
                )
            return 9 + self.sample

        # sample is a name
        if self.sample in sample_names:
            return 9 + sample_names.index(self.sample)
        raise ValueError(
            f"Sample '{self.sample}' not found. "
            f"Available samples: {sample_names}"
        )

    def _parse_data_row(
        self,
        line: str,
        sample_idx: int,
    ) -> dict[str, object] | None:
        """Parse a single VCF data row, returning None if the row should be skipped."""
        parts = line.strip().split("\t")

        # Need at least up to the sample column
        if len(parts) <= sample_idx:
            self.warnings.append(f"Row too short, skipping: {line.strip()!r}")
            return None

        chrom = parts[0]

        # Position must be a valid integer
        try:
            pos = int(parts[1])
        except ValueError:
            self.warnings.append(
                f"Invalid position '{parts[1]}', skipping row"
            )
            return None

        # FILTER check
        filter_val = parts[6]
        if not self.include_filtered and filter_val not in ("PASS", "."):
            return None

        # rsid
        rsid = parts[2] if parts[2] != "." else f"{chrom}_{pos}"
        ref = parts[3]
        alts = parts[4].split(",")
        alleles_list = [ref] + alts

        # FORMAT / sample genotype
        format_fields = parts[8].split(":")
        if "GT" not in format_fields:
            self.warnings.append(
                f"No GT in FORMAT for row at {chrom}:{pos}, skipping"
            )
            return None
        gt_field_idx = format_fields.index("GT")

        sample_fields = parts[sample_idx].split(":")
        gt = sample_fields[gt_field_idx]
        gt_parts = gt.replace("|", "/").split("/")

        resolved: list[str] = []
        for g in gt_parts:
            try:
                idx = int(g)
            except ValueError:
                resolved.append("?")
                self.warnings.append(
                    f"Non-integer GT allele '{g}' at {chrom}:{pos}"
                )
                continue
            if idx < 0 or idx >= len(alleles_list):
                resolved.append("?")
                self.warnings.append(
                    f"GT index {idx} out of range at {chrom}:{pos} "
                    f"(alleles: {alleles_list})"
                )
            else:
                resolved.append(alleles_list[idx])

        allele1 = resolved[0] if resolved else "?"
        allele2 = resolved[1] if len(resolved) > 1 else allele1

        return {
            "rsid": rsid,
            "chromosome": chrom,
            "position": pos,
            "allele1": allele1,
            "allele2": allele2,
            "genotype": allele1 + allele2,
        }


# Registry of loaders
LOADERS: dict[DataFormat, BaseLoader] = {
    DataFormat.ANCESTRY_DNA: AncestryDNALoader(),
    DataFormat.TWENTY_THREE_ME: TwentyThreeMeLoader(),
    DataFormat.MY_HERITAGE: MyHeritageLoader(),
    DataFormat.VCF: VCFLoader(),
}


def detect_format(filepath: Path) -> Optional[DataFormat]:
    """Auto-detect the format of a DNA data file."""
    for format_type, loader in LOADERS.items():
        if loader.detect(filepath):
            return format_type
    return None


def load_dna_data(
    filepath: str | Path,
    format: DataFormat = DataFormat.AUTO,
) -> DNADataset:
    """
    Load DNA data from various file formats.
    
    Args:
        filepath: Path to the DNA data file
        format: Data format (auto-detected if not specified)
    
    Returns:
        DNADataset containing the loaded genotypes
    
    Raises:
        ValueError: If format cannot be detected or is unsupported
    """
    filepath = Path(filepath)
    
    if not filepath.exists():
        raise FileNotFoundError(f"DNA data file not found: {filepath}")
    
    # Auto-detect format if needed
    if format == DataFormat.AUTO:
        detected = detect_format(filepath)
        if detected is None:
            raise ValueError(
                f"Could not auto-detect format for {filepath}. "
                "Please specify format explicitly."
            )
        format = detected
    
    # Load data
    loader = LOADERS.get(format)
    if loader is None:
        raise ValueError(f"Unsupported format: {format}")
    
    df = loader.load(filepath)
    
    return DNADataset(
        format=format,
        source_file=str(filepath),
        snp_count=len(df),
        dataframe=df,
    )
