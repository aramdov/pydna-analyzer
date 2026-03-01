"""Tests for VCF file loader."""

from __future__ import annotations

import gzip

from pydna_analyzer.core.data_loader import (
    DataFormat,
    VCFLoader,
    detect_format,
)


class TestVCFDetection:
    """Tests for VCF format detection."""

    def test_detect_vcf_file(self, sample_vcf_file):
        fmt = detect_format(sample_vcf_file)
        assert fmt == DataFormat.VCF

    def test_detect_vcf_gz(self, tmp_path, sample_vcf_content):
        fp = tmp_path / "test.vcf.gz"
        with gzip.open(fp, "wt") as f:
            f.write(sample_vcf_content)
        fmt = detect_format(fp)
        assert fmt == DataFormat.VCF

    def test_rejects_non_vcf(self, sample_ancestrydna_file):
        loader = VCFLoader()
        assert loader.detect(sample_ancestrydna_file) is False


class TestVCFBasicParsing:
    """Tests for basic VCF parsing correctness."""

    def test_load_returns_dataframe_columns(self, sample_vcf_file):
        loader = VCFLoader()
        df = loader.load(sample_vcf_file)
        assert set(df.columns) >= {
            "rsid",
            "chromosome",
            "position",
            "allele1",
            "allele2",
            "genotype",
        }

    def test_snp_count(self, sample_vcf_file):
        loader = VCFLoader()
        df = loader.load(sample_vcf_file)
        assert len(df) == 5

    def test_rsid_extraction(self, sample_vcf_file):
        loader = VCFLoader()
        df = loader.load(sample_vcf_file)
        rsids = set(df["rsid"])
        assert "rs429358" in rsids
        assert "rs7412" in rsids
        assert "rs1801133" in rsids

    def test_synthetic_rsid_for_missing(self, sample_vcf_file):
        """Row with '.' rsid should get chrom_pos synthetic ID."""
        loader = VCFLoader()
        df = loader.load(sample_vcf_file)
        rsids = set(df["rsid"])
        assert "6_26091179" in rsids

    def test_genotype_resolution(self, sample_vcf_file):
        """0/1 with REF=T ALT=C should give allele1=T, allele2=C."""
        loader = VCFLoader()
        df = loader.load(sample_vcf_file)
        row = df[df["rsid"] == "rs429358"].iloc[0]
        assert row["allele1"] == "T"
        assert row["allele2"] == "C"
        assert row["genotype"] == "TC"

    def test_homozygous_ref(self, sample_vcf_file):
        """0/0 should give REF/REF."""
        loader = VCFLoader()
        df = loader.load(sample_vcf_file)
        row = df[df["rsid"] == "rs7412"].iloc[0]
        assert row["allele1"] == "C"
        assert row["allele2"] == "C"

    def test_homozygous_alt(self, sample_vcf_file):
        """1/1 should give ALT/ALT."""
        loader = VCFLoader()
        df = loader.load(sample_vcf_file)
        row = df[df["rsid"] == "rs1801133"].iloc[0]
        assert row["allele1"] == "T"
        assert row["allele2"] == "T"

    def test_chromosome_and_position(self, sample_vcf_file):
        loader = VCFLoader()
        df = loader.load(sample_vcf_file)
        row = df[df["rsid"] == "rs429358"].iloc[0]
        assert row["chromosome"] == "19"
        assert row["position"] == 45411941
