# Dependency Policy

This document defines the supported versions for all pipeline dependencies.
The validation script (`scripts/validate_deps.sh`) enforces these requirements.

## Version Policy

We follow a **minimum version with range** policy:
- **Minimum** versions are tested and known to work
- **Maximum** versions are set conservatively to avoid untested breaking changes
- Patch/minor updates within ranges are generally safe

## Core Tools (Required)

| Tool | Minimum | Maximum | Notes |
|------|---------|---------|-------|
| `fastqc` | 0.11.9 | 0.12.x | QC reports |
| `fastp` | 0.20.0 | 1.x | Read trimming (v1.0 is latest major) |
| `bwa-mem2` | 2.2.1 | 2.x | Alignment (BWA-MEM2) |
| `samtools` | 1.15 | 1.x | BAM processing (htslib-based) |
| `bcftools` | 1.15 | 1.x | Variant calling (htslib-based) |

## Optional Tools

| Tool | Minimum | Maximum | Notes |
|------|---------|---------|-------|
| `vep` | 105 | latest | Annotation (optional) |

## Runtime Dependencies

| Tool | Minimum | Notes |
|------|---------|-------|
| `bash` | 3.2 | Shell (3.2 for macOS compat; 4.0+ preferred) |
| `python3` | 3.8 | Test scripts |
| `conda` | 4.10 | Environment management |

> **Note on version ranges**: Maximum versions are set permissively to allow minor/patch updates.
> The validation script will warn if a version is newer than tested. Review DEPENDENCIES.md
> when major versions are released.

## Docker Base Image

- **Base:** `ubuntu:22.04`
- **Miniconda:** `Miniconda3-py39_24.3.0-0`

## CI Environment Matrix

The CI tests against:
- Ubuntu 22.04 (primary)
- macOS 14 (ARM64)

## Version Pinning Strategy

### Conda Environment

The `.github/conda/integration-environment.yml` uses flexible pinning:
```yaml
dependencies:
  - python=3.11
  - fastqc>=0.11.9
  - fastp>=0.20.0
  - bwa-mem2>=2.2.1
  - samtools>=1.15
  - bcftools>=1.15
```

### Docker

The Dockerfile pins to known-good versions via conda ranges.
Build-time validation ensures no legacy versions slip through.

## Updating Dependencies

1. Update version ranges in this file
2. Update `scripts/validate_deps.sh` with new ranges
3. Update `.github/conda/integration-environment.yml`
4. Run `bash tests/run_tests.sh` to verify
5. Test Docker build: `docker build -t wgs-pipeline:test .`

## Version Validation

Run the validation script to check all dependencies:

```bash
# Full validation (requires all tools installed)
bash scripts/validate_deps.sh

# CI-friendly mode (skip unavailable tools)
bash scripts/validate_deps.sh --ci-mode

# Check specific tool
bash scripts/validate_deps.sh --tool samtools
```

The script exits with:
- `0` - All versions within acceptable ranges
- `1` - One or more versions out of range
- `2` - Required tool missing
