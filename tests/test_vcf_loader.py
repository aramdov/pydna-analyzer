"""Tests for VCF file loader."""

from __future__ import annotations

import gzip

import pytest

from pydna_analyzer.core.data_loader import (
    DataFormat,
    VCFLoader,
    detect_format,
    load_dna_data,
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


class TestVCFMultiAllelic:
    """Tests for multi-allelic site handling."""

    def test_triallelic_het(self, tmp_path, multi_allelic_vcf_content):
        """GT 1/2 with ALT=G,T should give G/T."""
        fp = tmp_path / "multi.vcf"
        fp.write_text(multi_allelic_vcf_content)
        loader = VCFLoader()
        df = loader.load(fp)
        row = df[df["rsid"] == "rs111"].iloc[0]
        assert row["allele1"] == "G"
        assert row["allele2"] == "T"
        assert row["genotype"] == "GT"

    def test_triallelic_ref_alt3(self, tmp_path, multi_allelic_vcf_content):
        """GT 0/3 with ALT=G,T,C should give A/C."""
        fp = tmp_path / "multi.vcf"
        fp.write_text(multi_allelic_vcf_content)
        loader = VCFLoader()
        df = loader.load(fp)
        row = df[df["rsid"] == "rs222"].iloc[0]
        assert row["allele1"] == "A"
        assert row["allele2"] == "C"

    def test_biallelic_still_works(self, tmp_path, multi_allelic_vcf_content):
        """Standard biallelic 0/1 should still work."""
        fp = tmp_path / "multi.vcf"
        fp.write_text(multi_allelic_vcf_content)
        loader = VCFLoader()
        df = loader.load(fp)
        row = df[df["rsid"] == "rs333"].iloc[0]
        assert row["allele1"] == "A"
        assert row["allele2"] == "G"

    def test_hom_ref_multiallelic(self, tmp_path, multi_allelic_vcf_content):
        """GT 0/0 at multi-allelic site should give REF/REF."""
        fp = tmp_path / "multi.vcf"
        fp.write_text(multi_allelic_vcf_content)
        loader = VCFLoader()
        df = loader.load(fp)
        row = df[df["rsid"] == "rs444"].iloc[0]
        assert row["allele1"] == "A"
        assert row["allele2"] == "A"

    def test_gt_index_out_of_range(self, tmp_path):
        """GT index beyond alleles list should produce '?' and a warning."""
        content = (
            "##fileformat=VCFv4.3\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
            "1\t100\trs999\tA\tG\t50\tPASS\t.\tGT\t0/5\n"
        )
        fp = tmp_path / "bad_gt.vcf"
        fp.write_text(content)
        loader = VCFLoader()
        df = loader.load(fp)
        row = df.iloc[0]
        assert row["allele2"] == "?"
        assert len(loader.warnings) == 1
        assert "out of range" in loader.warnings[0]


class TestVCFFilterHandling:
    """Tests for FILTER column handling."""

    def test_default_skips_non_pass(self, tmp_path, filtered_vcf_content):
        """Default: only PASS and '.' rows loaded."""
        fp = tmp_path / "filtered.vcf"
        fp.write_text(filtered_vcf_content)
        loader = VCFLoader()
        df = loader.load(fp)
        assert len(df) == 2
        rsids = set(df["rsid"])
        assert "rs100" in rsids
        assert "rs300" in rsids
        assert "rs200" not in rsids

    def test_include_filtered_keeps_all(self, tmp_path, filtered_vcf_content):
        """include_filtered=True loads everything."""
        fp = tmp_path / "filtered.vcf"
        fp.write_text(filtered_vcf_content)
        loader = VCFLoader(include_filtered=True)
        df = loader.load(fp)
        assert len(df) == 4

    def test_dot_filter_treated_as_pass(self, tmp_path, filtered_vcf_content):
        """'.' in FILTER column should be treated like PASS."""
        fp = tmp_path / "filtered.vcf"
        fp.write_text(filtered_vcf_content)
        loader = VCFLoader()
        df = loader.load(fp)
        rsids = set(df["rsid"])
        assert "rs300" in rsids


class TestVCFMultiSample:
    """Tests for multi-sample VCF handling."""

    def test_default_selects_first_sample(self, tmp_path, multi_sample_vcf_content):
        """Default (sample=None) should select the first sample (NA12878)."""
        fp = tmp_path / "multi.vcf"
        fp.write_text(multi_sample_vcf_content)
        loader = VCFLoader()
        df = loader.load(fp)
        row = df[df["rsid"] == "rs429358"].iloc[0]
        assert row["genotype"] == "TC"

    def test_select_by_name(self, tmp_path, multi_sample_vcf_content):
        """Select sample by name."""
        fp = tmp_path / "multi.vcf"
        fp.write_text(multi_sample_vcf_content)
        loader = VCFLoader(sample="HG00096")
        df = loader.load(fp)
        row = df[df["rsid"] == "rs429358"].iloc[0]
        assert row["genotype"] == "CC"

    def test_select_by_index(self, tmp_path, multi_sample_vcf_content):
        """Select sample by 0-based index."""
        fp = tmp_path / "multi.vcf"
        fp.write_text(multi_sample_vcf_content)
        loader = VCFLoader(sample=2)
        df = loader.load(fp)
        row = df[df["rsid"] == "rs1801133"].iloc[0]
        assert row["genotype"] == "TT"

    def test_invalid_sample_name_raises(self, tmp_path, multi_sample_vcf_content):
        """Non-existent sample name should raise ValueError."""
        fp = tmp_path / "multi.vcf"
        fp.write_text(multi_sample_vcf_content)
        loader = VCFLoader(sample="DOES_NOT_EXIST")
        with pytest.raises(ValueError, match="not found"):
            loader.load(fp)

    def test_invalid_sample_index_raises(self, tmp_path, multi_sample_vcf_content):
        """Out-of-range sample index should raise ValueError."""
        fp = tmp_path / "multi.vcf"
        fp.write_text(multi_sample_vcf_content)
        loader = VCFLoader(sample=99)
        with pytest.raises(ValueError, match="out of range"):
            loader.load(fp)


class TestVCFValidation:
    """Tests for VCF structural validation."""

    def test_missing_fileformat_raises(self, tmp_path, malformed_vcf_content_no_header):
        fp = tmp_path / "bad.vcf"
        fp.write_text(malformed_vcf_content_no_header)
        loader = VCFLoader()
        with pytest.raises(ValueError, match="missing ##fileformat"):
            loader.load(fp)

    def test_missing_chrom_header_raises(self, tmp_path, malformed_vcf_content_no_chrom):
        fp = tmp_path / "bad.vcf"
        fp.write_text(malformed_vcf_content_no_chrom)
        loader = VCFLoader()
        with pytest.raises(ValueError, match="missing #CHROM"):
            loader.load(fp)

    def test_empty_file_raises(self, tmp_path):
        fp = tmp_path / "empty.vcf"
        fp.write_text("")
        loader = VCFLoader()
        with pytest.raises(ValueError, match="missing ##fileformat"):
            loader.load(fp)

    def test_headers_only_returns_empty_df(self, tmp_path):
        """VCF with headers but no data rows should return empty DataFrame."""
        content = (
            "##fileformat=VCFv4.3\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
        )
        fp = tmp_path / "headers_only.vcf"
        fp.write_text(content)
        loader = VCFLoader()
        df = loader.load(fp)
        assert len(df) == 0


class TestVCFEdgeCases:
    """Tests for edge cases and warning collection."""

    def test_bad_rows_skipped_with_warnings(self, tmp_path, malformed_vcf_content_bad_rows):
        """Malformed rows should be skipped, good rows kept, warnings collected."""
        fp = tmp_path / "bad_rows.vcf"
        fp.write_text(malformed_vcf_content_bad_rows)
        loader = VCFLoader()
        df = loader.load(fp)
        assert len(df) == 3
        assert len(loader.warnings) >= 3

    def test_gzipped_vcf(self, tmp_path, sample_vcf_content):
        """Gzipped VCF should load identically to uncompressed."""
        fp = tmp_path / "test.vcf.gz"
        with gzip.open(fp, "wt") as f:
            f.write(sample_vcf_content)
        loader = VCFLoader()
        df = loader.load(fp)
        assert len(df) == 5

    def test_phased_genotype(self, tmp_path):
        """Phased genotypes (|) should be handled same as unphased (/)."""
        content = (
            "##fileformat=VCFv4.3\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
            "1\t100\trs100\tA\tG\t50\tPASS\t.\tGT\t0|1\n"
        )
        fp = tmp_path / "phased.vcf"
        fp.write_text(content)
        loader = VCFLoader()
        df = loader.load(fp)
        row = df.iloc[0]
        assert row["allele1"] == "A"
        assert row["allele2"] == "G"

    def test_warnings_reset_on_each_load(
        self, tmp_path, malformed_vcf_content_bad_rows, sample_vcf_content
    ):
        """Warnings should be cleared on each call to load()."""
        loader = VCFLoader()
        fp1 = tmp_path / "bad.vcf"
        fp1.write_text(malformed_vcf_content_bad_rows)
        loader.load(fp1)
        assert len(loader.warnings) > 0

        fp2 = tmp_path / "good.vcf"
        fp2.write_text(sample_vcf_content)
        loader.load(fp2)
        assert len(loader.warnings) == 0


class TestVCFIntegration:
    """Integration tests: VCF through load_dna_data()."""

    def test_auto_detect_vcf(self, sample_vcf_file):
        """load_dna_data should auto-detect VCF and return DNADataset."""
        dataset = load_dna_data(sample_vcf_file)
        assert dataset.format == DataFormat.VCF
        assert dataset.snp_count == 5

    def test_vcf_dataset_get_genotype(self, sample_vcf_file):
        """DNADataset.get_genotype() should work with VCF-loaded data."""
        dataset = load_dna_data(sample_vcf_file)
        assert dataset.get_genotype("rs429358") == "TC"
        assert dataset.get_genotype("rs7412") == "CC"
