# VCF File Support Improvements — Design

**Date**: 2026-03-01
**Status**: Design approved

## Overview

Make the VCF parser robust enough to handle real-world VCF files correctly, with proper
validation, multi-allelic support, basic quality filtering, and comprehensive test coverage.

## Scope

### In scope

- **Validation** — Verify VCF structure on load (header, column count, genotype fields). Clear
  `ValueError` messages instead of crashing on malformed input.
- **Multi-allelic support** — Sites with multiple ALT alleles resolve GT indices against the full
  alleles list. GT `0/2` with `ALT=G,T` → allele1=A (ref), allele2=T. One row per site.
- **FILTER column** — Default: only load PASS or `.` rows. `include_filtered=True` overrides.
- **Multi-sample awareness** — Default to first sample (not last). `sample` parameter selects by
  name or 0-based index. Error if sample not found.
- **Test coverage** — ~25-30 tests covering all VCF parsing paths.
- **Gzip support** — Already works, keep it.

### Out of scope

- INFO field extraction
- Structural variants / symbolic alleles (`<DEL>`, `<DUP>`)
- Streaming/chunked parsing for very large files
- VCF writing/export

## API Changes

Parameters added to `VCFLoader.__init__()`, keeping the `BaseLoader` interface intact:

```python
class VCFLoader(BaseLoader):
    def __init__(
        self,
        sample: str | int | None = None,  # name or 0-based index; None = first
        include_filtered: bool = False,     # if True, include non-PASS rows
    ):
        self.sample = sample
        self.include_filtered = include_filtered
```

`load_dna_data()` auto-detection uses defaults (first sample, PASS-only). Direct instantiation
allows configuration:

```python
loader = VCFLoader(sample="NA12878", include_filtered=True)
data = loader.load("multi_sample.vcf.gz")
```

## Parsing Logic

1. **Validation layer** (runs first):
   - Verify `##fileformat=VCF` header exists
   - Verify `#CHROM` header line has at least 8 required columns
   - Validate requested sample exists in the header
   - Failure → `ValueError` with clear message

2. **FILTER check** (per-row):
   - Check column 6 (FILTER)
   - Skip non-PASS/`.` rows unless `include_filtered=True`

3. **Multi-allelic resolution** (per-row):
   - Parse all ALT alleles: `alts = parts[4].split(",")`
   - Build alleles list: `[ref] + alts` (index 0 = ref, 1+ = alts)
   - Parse GT indices and map to alleles list
   - One row per site, correctly resolved

## Error Handling

| Scenario | Behavior |
|----------|----------|
| No `##fileformat=VCF` header | `ValueError` |
| Missing `#CHROM` header line | `ValueError` |
| Invalid sample name/index | `ValueError` |
| Row with wrong column count | Skip row, collect warning |
| GT field missing from FORMAT | Skip row, collect warning |
| GT index out of range | Set allele to `"?"`, collect warning |
| Empty file / headers only | Return empty DataFrame |
| Invalid position (non-integer) | Skip row, collect warning |

Warnings collected on `loader.warnings` list (keeps `load()` return type as DataFrame).

## Test Strategy

String fixtures in `conftest.py`, written to `tmp_path`. No binary VCF files committed.

| Test Class | Coverage |
|------------|----------|
| `TestVCFDetection` | `.vcf`, `.vcf.gz`, rejects non-VCF |
| `TestVCFBasicParsing` | Columns, rsID, synthetic rsID, genotype resolution |
| `TestVCFMultiAllelic` | Biallelic, tri-allelic GT indices, homozygous ref, out-of-range |
| `TestVCFFilterHandling` | Default skips non-PASS, include_filtered keeps all, `.` as PASS |
| `TestVCFMultiSample` | First sample default, select by name/index, invalid → ValueError |
| `TestVCFValidation` | Missing headers → error, empty file → empty DataFrame |
| `TestVCFEdgeCases` | Wrong column count, missing GT, warnings collected, gzip |
| `TestVCFIntegration` | Auto-detection, compatibility with ClinicalAnalyzer |
