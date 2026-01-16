"""Tests for DNA data loaders."""

import pytest
from pathlib import Path

from genomeinsight.core.data_loader import (
    load_dna_data,
    detect_format,
    DataFormat,
    DNADataset,
    AncestryDNALoader,
    TwentyThreeMeLoader,
)


class TestDataFormat:
    """Tests for format detection."""

    def test_detect_ancestrydna(self, sample_ancestrydna_file):
        """Should detect AncestryDNA format."""
        fmt = detect_format(sample_ancestrydna_file)
        assert fmt == DataFormat.ANCESTRY_DNA

    def test_detect_23andme(self, sample_23andme_file):
        """Should detect 23andMe format."""
        fmt = detect_format(sample_23andme_file)
        assert fmt == DataFormat.TWENTY_THREE_ME


class TestAncestryDNALoader:
    """Tests for AncestryDNA format loading."""

    def test_load_returns_dataset(self, sample_ancestrydna_file):
        """Should return a DNADataset."""
        dataset = load_dna_data(sample_ancestrydna_file)
        assert isinstance(dataset, DNADataset)
        assert dataset.format == DataFormat.ANCESTRY_DNA

    def test_snp_count(self, sample_ancestrydna_file):
        """Should count SNPs correctly."""
        dataset = load_dna_data(sample_ancestrydna_file)
        assert dataset.snp_count == 10  # 10 SNPs in sample data

    def test_get_genotype(self, sample_ancestrydna_file):
        """Should retrieve genotypes by rsID."""
        dataset = load_dna_data(sample_ancestrydna_file)

        # Test specific genotypes
        assert dataset.get_genotype("rs429358") == "TT"
        assert dataset.get_genotype("rs7412") == "CC"
        assert dataset.get_genotype("rs1801133") == "CT"  # Heterozygous

    def test_get_genotype_missing(self, sample_ancestrydna_file):
        """Should return None for missing rsIDs."""
        dataset = load_dna_data(sample_ancestrydna_file)
        assert dataset.get_genotype("rs999999999") is None

    def test_get_rsids(self, sample_ancestrydna_file):
        """Should return all rsIDs."""
        dataset = load_dna_data(sample_ancestrydna_file)
        rsids = dataset.get_rsids()

        assert "rs429358" in rsids
        assert "rs7412" in rsids
        assert "rs1801133" in rsids
        assert len(rsids) == 10


class TestTwentyThreeMeLoader:
    """Tests for 23andMe format loading."""

    def test_load_returns_dataset(self, sample_23andme_file):
        """Should return a DNADataset."""
        dataset = load_dna_data(sample_23andme_file)
        assert isinstance(dataset, DNADataset)
        assert dataset.format == DataFormat.TWENTY_THREE_ME

    def test_snp_count(self, sample_23andme_file):
        """Should count SNPs correctly."""
        dataset = load_dna_data(sample_23andme_file)
        assert dataset.snp_count == 5  # 5 SNPs in 23andMe sample

    def test_genotype_parsing(self, sample_23andme_file):
        """23andMe uses combined genotype column (e.g., 'CT' not 'C' 'T')."""
        dataset = load_dna_data(sample_23andme_file)

        # Should split correctly
        assert dataset.get_genotype("rs429358") == "TT"
        assert dataset.get_genotype("rs1801133") == "CT"


class TestLoadDNAData:
    """Tests for the main load_dna_data function."""

    def test_auto_detect_ancestry(self, sample_ancestrydna_file):
        """Should auto-detect AncestryDNA format."""
        dataset = load_dna_data(sample_ancestrydna_file)
        assert dataset.format == DataFormat.ANCESTRY_DNA

    def test_auto_detect_23andme(self, sample_23andme_file):
        """Should auto-detect 23andMe format."""
        dataset = load_dna_data(sample_23andme_file)
        assert dataset.format == DataFormat.TWENTY_THREE_ME

    def test_explicit_format(self, sample_ancestrydna_file):
        """Should respect explicit format specification."""
        dataset = load_dna_data(
            sample_ancestrydna_file,
            format=DataFormat.ANCESTRY_DNA
        )
        assert dataset.format == DataFormat.ANCESTRY_DNA

    def test_source_file_stored(self, sample_ancestrydna_file):
        """Should store source file path."""
        dataset = load_dna_data(sample_ancestrydna_file)
        assert sample_ancestrydna_file.name in dataset.source_file

    def test_missing_file_raises(self):
        """Should raise error for missing file."""
        with pytest.raises((FileNotFoundError, ValueError)):
            load_dna_data(Path("/nonexistent/file.txt"))
