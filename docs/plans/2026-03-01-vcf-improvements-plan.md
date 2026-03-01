# VCF File Support Improvements — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make VCFLoader robust enough to handle real-world VCF files with proper validation, multi-allelic support, FILTER-based quality filtering, multi-sample selection, and comprehensive test coverage.

**Architecture:** Rewrite `VCFLoader` in `pydna_analyzer/core/data_loader.py` to add `__init__()` parameters (`sample`, `include_filtered`), a validation layer, per-row error handling with warning collection, and correct multi-allelic genotype resolution. Add ~25-30 tests in a new `tests/test_vcf_loader.py`. Update the `LOADERS` registry to use a factory function for configured loaders.

**Tech Stack:** Python 3.11+, pandas, gzip, pytest, pydantic (existing DNADataset model)

---

### Task 1: Add VCF test fixtures to conftest.py

**Files:**
- Modify: `tests/conftest.py` (append after line 251)

**Step 1: Add VCF fixture content strings**

Add these fixtures at the end of `tests/conftest.py`:

```python
# ===== VCF loader fixtures =====

VCF_HEADER = (
    "##fileformat=VCFv4.3\n"
    "##INFO=<ID=DP,Number=1,Type=Integer>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
)


@pytest.fixture
def sample_vcf_content():
    """Simple valid VCF: single sample, biallelic, all PASS."""
    return (
        VCF_HEADER
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
        "19\t45411941\trs429358\tT\tC\t50\tPASS\t.\tGT\t0/1\n"
        "19\t45412079\trs7412\tC\tT\t50\tPASS\t.\tGT\t0/0\n"
        "1\t11856378\trs1801133\tC\tT\t50\tPASS\t.\tGT\t1/1\n"
        "1\t11854476\trs1801131\tA\tC\t50\tPASS\t.\tGT\t0/1\n"
        "6\t26091179\t.\tC\tG\t50\tPASS\t.\tGT\t0/0\n"
    )


@pytest.fixture
def sample_vcf_file(tmp_path, sample_vcf_content):
    """Write sample VCF to a temp file."""
    fp = tmp_path / "test.vcf"
    fp.write_text(sample_vcf_content)
    return fp


@pytest.fixture
def multi_sample_vcf_content():
    """VCF with 3 samples."""
    return (
        VCF_HEADER
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\tHG00096\tHG00097\n"
        "19\t45411941\trs429358\tT\tC\t50\tPASS\t.\tGT\t0/1\t1/1\t0/0\n"
        "1\t11856378\trs1801133\tC\tT\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1\n"
    )


@pytest.fixture
def multi_allelic_vcf_content():
    """VCF with multi-allelic sites (2-3 ALT alleles)."""
    return (
        VCF_HEADER
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
        "1\t100\trs111\tA\tG,T\t50\tPASS\t.\tGT\t1/2\n"     # het: G/T
        "1\t200\trs222\tA\tG,T,C\t50\tPASS\t.\tGT\t0/3\n"   # het: A/C
        "1\t300\trs333\tA\tG\t50\tPASS\t.\tGT\t0/1\n"        # biallelic: A/G
        "1\t400\trs444\tA\tG,T\t50\tPASS\t.\tGT\t0/0\n"      # hom ref: A/A
    )


@pytest.fixture
def filtered_vcf_content():
    """VCF with mix of PASS, ., and non-PASS filters."""
    return (
        VCF_HEADER
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
        "1\t100\trs100\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
        "1\t200\trs200\tC\tT\t10\tLowQual\t.\tGT\t0/1\n"
        "1\t300\trs300\tG\tA\t50\t.\t.\tGT\t1/1\n"
        "1\t400\trs400\tT\tC\t5\tLowQual;LowDP\t.\tGT\t0/1\n"
    )


@pytest.fixture
def malformed_vcf_content_no_header():
    """VCF missing ##fileformat header."""
    return (
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
        "1\t100\trs100\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
    )


@pytest.fixture
def malformed_vcf_content_no_chrom():
    """VCF with ##fileformat but missing #CHROM line."""
    return (
        "##fileformat=VCFv4.3\n"
        "1\t100\trs100\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
    )


@pytest.fixture
def malformed_vcf_content_bad_rows():
    """VCF with valid headers but some broken data rows."""
    return (
        VCF_HEADER
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
        "1\t100\trs100\tA\tG\t50\tPASS\t.\tGT\t0/1\n"    # good
        "1\tnot_a_number\trs200\tA\tG\t50\tPASS\t.\tGT\t0/1\n"  # bad position
        "1\t300\n"                                          # too few columns
        "1\t400\trs400\tA\tG\t50\tPASS\t.\tDP\t42\n"       # no GT in FORMAT
        "1\t500\trs500\tA\tG\t50\tPASS\t.\tGT\t0/9\n"     # GT index out of range
        "1\t600\trs600\tC\tT\t50\tPASS\t.\tGT\t0/1\n"     # good
    )
```

**Step 2: Run existing tests to verify nothing is broken**

Run: `uv run pytest tests/conftest.py --collect-only -q`
Expected: Fixtures collected, no errors

**Step 3: Commit**

```bash
git add tests/conftest.py
git commit -m "test(vcf): add VCF test fixtures to conftest"
```

---

### Task 2: Create test file with basic parsing tests

**Files:**
- Create: `tests/test_vcf_loader.py`

**Step 1: Write the basic parsing tests**

Create `tests/test_vcf_loader.py`:

```python
"""Tests for VCF file loader."""

from __future__ import annotations

import gzip
from pathlib import Path

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
        assert set(df.columns) >= {"rsid", "chromosome", "position", "allele1", "allele2", "genotype"}

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
```

**Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_vcf_loader.py -v`
Expected: Tests FAIL because `VCFLoader()` constructor doesn't accept no-args cleanly yet (or some tests pass with current code — that's fine for tests that match current behavior). The key point is these tests exist and define the expected behavior.

**Step 3: Commit**

```bash
git add tests/test_vcf_loader.py
git commit -m "test(vcf): add basic VCF parsing tests"
```

---

### Task 3: Rewrite VCFLoader with validation and __init__ parameters

This is the core implementation task. Rewrite the `VCFLoader` class in `pydna_analyzer/core/data_loader.py` (lines 234-303).

**Files:**
- Modify: `pydna_analyzer/core/data_loader.py:234-312`

**Step 1: Replace VCFLoader class with the new implementation**

Replace lines 234-303 (`class VCFLoader` through end of `load` method) with:

```python
class VCFLoader(BaseLoader):
    """Loader for VCF (Variant Call Format) files.

    Args:
        sample: Sample name (str) or 0-based index (int) to extract.
                None means first sample.
        include_filtered: If True, include rows where FILTER is not PASS.
                         Default False: only PASS and '.' rows are loaded.
    """

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
            ValueError: If VCF structure is invalid (missing headers, bad sample).
        """
        self.warnings = []
        opener = gzip.open if str(filepath).endswith(".gz") else open

        sample_names: list[str] = []
        sample_col_idx: int | None = None
        found_fileformat = False
        found_chrom = False
        rows: list[dict] = []

        with opener(filepath, "rt") as f:
            for line_num, line in enumerate(f, 1):
                line = line.rstrip("\n\r")

                # Meta-information lines
                if line.startswith("##"):
                    if line.startswith("##fileformat=VCF"):
                        found_fileformat = True
                    continue

                # Header line
                if line.startswith("#CHROM"):
                    found_chrom = True
                    headers = line.split("\t")
                    if len(headers) < 8:
                        raise ValueError(
                            f"Invalid VCF: #CHROM header has {len(headers)} columns, "
                            "expected at least 8"
                        )
                    # Samples start at column 9 (0-indexed)
                    sample_names = headers[9:] if len(headers) > 9 else []
                    sample_col_idx = self._resolve_sample_index(sample_names)
                    continue

                # Validate we saw required headers
                if not found_fileformat:
                    raise ValueError(
                        "Not a valid VCF file: missing ##fileformat header"
                    )
                if not found_chrom:
                    raise ValueError(
                        "Invalid VCF: missing #CHROM header line"
                    )

                # Data row
                row = self._parse_data_row(
                    line, line_num, sample_col_idx, sample_names
                )
                if row is not None:
                    rows.append(row)

        # Edge case: headers only, no data rows — still need validation
        if not found_fileformat:
            raise ValueError(
                "Not a valid VCF file: missing ##fileformat header"
            )
        if not found_chrom:
            raise ValueError("Invalid VCF: missing #CHROM header line")

        return pd.DataFrame(rows)

    def _resolve_sample_index(self, sample_names: list[str]) -> int | None:
        """Resolve sample parameter to a column index (0-based within sample columns)."""
        if not sample_names:
            return None

        if self.sample is None:
            return 0  # default: first sample

        if isinstance(self.sample, int):
            if self.sample < 0 or self.sample >= len(sample_names):
                raise ValueError(
                    f"Sample index {self.sample} out of range. "
                    f"VCF has {len(sample_names)} samples: {', '.join(sample_names)}"
                )
            return self.sample

        # String: lookup by name
        if self.sample not in sample_names:
            raise ValueError(
                f"Sample '{self.sample}' not found in VCF. "
                f"Available: {', '.join(sample_names)}"
            )
        return sample_names.index(self.sample)

    def _parse_data_row(
        self,
        line: str,
        line_num: int,
        sample_col_idx: int | None,
        sample_names: list[str],
    ) -> dict | None:
        """Parse a single VCF data row. Returns None if row should be skipped."""
        parts = line.split("\t")

        # Validate column count (need at least 8 fixed columns)
        if len(parts) < 8:
            self.warnings.append(
                f"Line {line_num}: skipped — only {len(parts)} columns (expected 8+)"
            )
            return None

        # FILTER check
        filter_val = parts[6]
        if not self.include_filtered and filter_val not in ("PASS", "."):
            return None

        # Parse position
        try:
            pos = int(parts[1])
        except ValueError:
            self.warnings.append(
                f"Line {line_num}: skipped — invalid position '{parts[1]}'"
            )
            return None

        chrom = parts[0]
        rsid = parts[2] if parts[2] != "." else f"{chrom}_{pos}"
        ref = parts[3]
        alts = parts[4].split(",")
        alleles_list = [ref] + alts  # index 0 = ref, 1+ = alts

        # Parse genotype from sample column
        if sample_col_idx is not None and len(parts) > 9:
            actual_col = 9 + sample_col_idx
            if actual_col >= len(parts):
                self.warnings.append(
                    f"Line {line_num}: skipped — sample column {sample_col_idx} missing"
                )
                return None

            format_fields = parts[8].split(":")
            if "GT" not in format_fields:
                self.warnings.append(
                    f"Line {line_num}: skipped — no GT field in FORMAT '{parts[8]}'"
                )
                return None

            gt_field_idx = format_fields.index("GT")
            sample_fields = parts[actual_col].split(":")

            if gt_field_idx >= len(sample_fields):
                self.warnings.append(
                    f"Line {line_num}: skipped — GT field index out of range in sample data"
                )
                return None

            gt = sample_fields[gt_field_idx]
            gt_indices = gt.replace("|", "/").split("/")

            resolved = []
            for idx_str in gt_indices:
                if idx_str == ".":
                    resolved.append(".")
                    continue
                try:
                    idx = int(idx_str)
                except ValueError:
                    resolved.append("?")
                    self.warnings.append(
                        f"Line {line_num}: invalid GT index '{idx_str}'"
                    )
                    continue
                if idx < 0 or idx >= len(alleles_list):
                    resolved.append("?")
                    self.warnings.append(
                        f"Line {line_num}: GT index {idx} out of range "
                        f"(have {len(alleles_list)} alleles)"
                    )
                else:
                    resolved.append(alleles_list[idx])

            allele1 = resolved[0] if resolved else "."
            allele2 = resolved[1] if len(resolved) > 1 else allele1
        else:
            # No sample columns — default to homozygous reference
            allele1, allele2 = ref, ref

        return {
            "rsid": rsid,
            "chromosome": chrom,
            "position": pos,
            "allele1": allele1,
            "allele2": allele2,
            "genotype": allele1 + allele2,
        }
```

**Step 2: Update the LOADERS registry**

The global `LOADERS` dict currently creates `VCFLoader()` at module level. Since `VCFLoader.__init__` now has default parameters, this still works. No change needed to the registry for the default case. The line at ~311:

```python
DataFormat.VCF: VCFLoader(),
```

This already works — `VCFLoader()` creates a loader with `sample=None, include_filtered=False`.

**Step 3: Run tests**

Run: `uv run pytest tests/test_vcf_loader.py -v`
Expected: All `TestVCFDetection` and `TestVCFBasicParsing` tests PASS

**Step 4: Run full suite to check for regressions**

Run: `uv run pytest -v`
Expected: All 188+ tests pass

**Step 5: Commit**

```bash
git add pydna_analyzer/core/data_loader.py
git commit -m "feat(vcf): rewrite VCFLoader with validation and multi-allelic support"
```

---

### Task 4: Add multi-allelic tests

**Files:**
- Modify: `tests/test_vcf_loader.py`

**Step 1: Add multi-allelic test class**

Append to `tests/test_vcf_loader.py`:

```python
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
```

**Step 2: Run tests**

Run: `uv run pytest tests/test_vcf_loader.py::TestVCFMultiAllelic -v`
Expected: All PASS

**Step 3: Commit**

```bash
git add tests/test_vcf_loader.py
git commit -m "test(vcf): add multi-allelic genotype resolution tests"
```

---

### Task 5: Add FILTER handling tests

**Files:**
- Modify: `tests/test_vcf_loader.py`

**Step 1: Add filter test class**

Append to `tests/test_vcf_loader.py`:

```python
class TestVCFFilterHandling:
    """Tests for FILTER column handling."""

    def test_default_skips_non_pass(self, tmp_path, filtered_vcf_content):
        """Default: only PASS and '.' rows loaded."""
        fp = tmp_path / "filtered.vcf"
        fp.write_text(filtered_vcf_content)
        loader = VCFLoader()
        df = loader.load(fp)
        # rs100 (PASS) and rs300 (.) should be loaded; rs200 and rs400 (LowQual) skipped
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
        assert "rs300" in rsids  # has '.' filter
```

**Step 2: Run tests**

Run: `uv run pytest tests/test_vcf_loader.py::TestVCFFilterHandling -v`
Expected: All PASS

**Step 3: Commit**

```bash
git add tests/test_vcf_loader.py
git commit -m "test(vcf): add FILTER column handling tests"
```

---

### Task 6: Add multi-sample tests

**Files:**
- Modify: `tests/test_vcf_loader.py`

**Step 1: Add multi-sample test class**

Append to `tests/test_vcf_loader.py`:

```python
class TestVCFMultiSample:
    """Tests for multi-sample VCF handling."""

    def test_default_selects_first_sample(self, tmp_path, multi_sample_vcf_content):
        """Default (sample=None) should select the first sample (NA12878)."""
        fp = tmp_path / "multi.vcf"
        fp.write_text(multi_sample_vcf_content)
        loader = VCFLoader()
        df = loader.load(fp)
        # NA12878 has rs429358 0/1 → T/C; HG00096 has 1/1 → C/C
        row = df[df["rsid"] == "rs429358"].iloc[0]
        assert row["genotype"] == "TC"

    def test_select_by_name(self, tmp_path, multi_sample_vcf_content):
        """Select sample by name."""
        fp = tmp_path / "multi.vcf"
        fp.write_text(multi_sample_vcf_content)
        loader = VCFLoader(sample="HG00096")
        df = loader.load(fp)
        row = df[df["rsid"] == "rs429358"].iloc[0]
        assert row["genotype"] == "CC"  # HG00096 is 1/1

    def test_select_by_index(self, tmp_path, multi_sample_vcf_content):
        """Select sample by 0-based index."""
        fp = tmp_path / "multi.vcf"
        fp.write_text(multi_sample_vcf_content)
        loader = VCFLoader(sample=2)  # HG00097
        df = loader.load(fp)
        row = df[df["rsid"] == "rs1801133"].iloc[0]
        assert row["genotype"] == "TT"  # HG00097 is 1/1

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
```

**Step 2: Run tests**

Run: `uv run pytest tests/test_vcf_loader.py::TestVCFMultiSample -v`
Expected: All PASS

**Step 3: Commit**

```bash
git add tests/test_vcf_loader.py
git commit -m "test(vcf): add multi-sample selection tests"
```

---

### Task 7: Add validation and edge case tests

**Files:**
- Modify: `tests/test_vcf_loader.py`

**Step 1: Add validation and edge case test classes**

Append to `tests/test_vcf_loader.py`:

```python
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
        # 2 good rows (rs100, rs600), others have issues
        assert len(df) == 3  # rs100 (good), rs500 (loaded with '?'), rs600 (good)
        assert len(loader.warnings) >= 3  # bad position, too few cols, no GT

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

    def test_warnings_reset_on_each_load(self, tmp_path, malformed_vcf_content_bad_rows, sample_vcf_content):
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
```

**Step 2: Run all VCF tests**

Run: `uv run pytest tests/test_vcf_loader.py -v`
Expected: All ~27 tests PASS

**Step 3: Run full suite**

Run: `uv run pytest -v`
Expected: All tests pass (188 existing + ~27 new)

**Step 4: Commit**

```bash
git add tests/test_vcf_loader.py
git commit -m "test(vcf): add validation, edge case, and integration tests"
```

---

### Task 8: Update docs and README

**Files:**
- Modify: `README.md`
- Modify: `CLAUDE.md`

**Step 1: Update README VCF section**

In `README.md`, the features table (line 25) currently says `📁 **Multi-Format Support** | AncestryDNA, 23andMe, MyHeritage, and VCF files`. Update this to:

```
| 📁 **Multi-Format Support** | AncestryDNA, 23andMe, MyHeritage, and VCF files (multi-sample, quality filtering) |
```

In the roadmap section (line 227), change:
```
- [ ] VCF file support improvements
```
to:
```
- [x] VCF file support improvements
```

**Step 2: Update CLAUDE.md test count**

In `CLAUDE.md`, update the test count in the "Test Structure" section to reflect the new VCF tests (will be ~215 tests total — update with actual count from test run).

**Step 3: Update README test count**

In `README.md` features table, update the test count (line 35) to the new total.

**Step 4: Commit**

```bash
git add README.md CLAUDE.md
git commit -m "docs: update README and CLAUDE.md for VCF improvements"
```

---

### Task summary

| Task | Description | Tests |
|------|-------------|-------|
| 1 | Add VCF fixtures to conftest.py | — |
| 2 | Create test file with basic parsing tests | 11 |
| 3 | Rewrite VCFLoader (core implementation) | — |
| 4 | Multi-allelic tests | 5 |
| 5 | FILTER handling tests | 3 |
| 6 | Multi-sample tests | 5 |
| 7 | Validation and edge case tests | 11 |
| 8 | Update docs | — |
| **Total** | | **~35 tests** |
