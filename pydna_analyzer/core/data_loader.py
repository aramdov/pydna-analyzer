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
    
    def detect(self, filepath: Path) -> bool:
        """Detect VCF format."""
        try:
            opener = gzip.open if str(filepath).endswith('.gz') else open
            with opener(filepath, 'rt') as f:
                for line in f:
                    if line.startswith('##fileformat=VCF'):
                        return True
                    if not line.startswith('#'):
                        break
        except Exception:
            pass
        return False
    
    def load(self, filepath: Path) -> pd.DataFrame:
        """Load VCF format file."""
        opener = gzip.open if str(filepath).endswith('.gz') else open
        
        rows = []
        with opener(filepath, 'rt') as f:
            for line in f:
                if line.startswith('##'):
                    continue
                if line.startswith('#CHROM'):
                    headers = line.strip().split('\t')
                    sample_col = headers[-1] if len(headers) > 9 else None
                    continue
                
                parts = line.strip().split('\t')
                chrom = parts[0]
                pos = int(parts[1])
                rsid = parts[2] if parts[2] != '.' else f"{chrom}_{pos}"
                ref = parts[3]
                alt = parts[4].split(',')[0]  # Take first alt allele
                
                # Parse genotype from sample column
                if len(parts) > 9:
                    format_fields = parts[8].split(':')
                    sample_fields = parts[9].split(':')
                    gt_idx = format_fields.index('GT') if 'GT' in format_fields else 0
                    gt = sample_fields[gt_idx]
                    
                    # Convert 0/1 format to alleles
                    gt_parts = gt.replace('|', '/').split('/')
                    alleles = []
                    for g in gt_parts:
                        if g == '0':
                            alleles.append(ref)
                        elif g == '1':
                            alleles.append(alt)
                        else:
                            alleles.append('.')
                    
                    allele1, allele2 = alleles[0], alleles[1] if len(alleles) > 1 else alleles[0]
                else:
                    allele1, allele2 = ref, ref
                
                rows.append({
                    'rsid': rsid,
                    'chromosome': chrom,
                    'position': pos,
                    'allele1': allele1,
                    'allele2': allele2,
                    'genotype': allele1 + allele2,
                })
        
        return pd.DataFrame(rows)


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
