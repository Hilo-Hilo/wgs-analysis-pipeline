# WGS Pipeline Smoke Tests

End-to-end smoke testing framework for the WGS analysis pipeline.

## Quick Start

```bash
# Run fastest smoke test (mock tools, ~10 seconds)
make smoke

# Run with real bioinformatics tools
make smoke-real

# Test all quality profiles
make smoke-all
```

## Overview

The smoke test harness provides:

1. **Synthetic Reference Genome** - Deterministic ~100kb genome with known variant positions
2. **Synthetic FASTQ Generator** - Paired-end reads with various quality profiles
3. **Full Pipeline Execution** - QC → Cleaning → Alignment → Variant Calling
4. **Output Validation** - Automated pass/fail checks

## Components

### `generate_reference.py`

Generates a tiny synthetic reference genome:
- Multiple "chromosomes" (4 contigs, ~100kb total)
- Known variant positions for validation
- Deterministic output (fixed seed)
- Ready for BWA-MEM2 indexing

```bash
python3 tests/smoke/generate_reference.py \
    --output-dir /path/to/output \
    --name my_ref \
    --seed 42
```

### `generate_fastq.py`

Generates synthetic paired-end FASTQ files:
- Quality profiles: `good`, `poor`, `adapter`, `mixed`
- Reads align to the synthetic reference
- Configurable read count and length
- Deterministic output

```bash
python3 tests/smoke/generate_fastq.py \
    --reference /path/to/ref.fa \
    --output-dir /path/to/output \
    --profile good \
    --num-reads 2000
```

### `run_smoke.sh`

Main smoke test harness:

```bash
# Mock mode (fast, no dependencies)
./tests/smoke/run_smoke.sh --mock

# Real mode (requires tools)
./tests/smoke/run_smoke.sh --real --verbose

# Test edge cases
./tests/smoke/run_smoke.sh --mock --profile poor --keep
```

**Options:**
- `--mock` - Use mock tools (fast, CI-safe)
- `--real` - Use real bioinformatics tools
- `--profile PROFILE` - Quality profile (good/poor/adapter/mixed)
- `--keep` - Keep output artifacts after completion
- `--verbose` - Enable verbose output

## Quality Profiles

| Profile | Description | Use Case |
|---------|-------------|----------|
| `good` | High quality reads (Q30+) | Normal operation testing |
| `poor` | Low quality reads | Edge case / error handling |
| `adapter` | Heavy adapter contamination | Trimming validation |
| `mixed` | Mix of quality levels | Realistic scenario |

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | All tests passed |
| 2 | Invalid usage |
| 10 | Setup/data generation failed |
| 11 | Quality control step failed |
| 12 | Data cleaning step failed |
| 13 | Alignment step failed |
| 14 | Variant calling step failed |
| 15 | Output validation failed |

## CI Integration

For CI environments:

```bash
# Quick validation (recommended for PRs)
make ci-quick

# Full validation (recommended for main branch)
make ci-full
```

## Artifacts

When using `--keep`, outputs are preserved at `/tmp/wgs-smoke-$PID/`:

```
/tmp/wgs-smoke-12345/
├── data/
│   ├── raw/                 # Generated FASTQ files
│   ├── processed/           # Cleaned FASTQ files
│   └── reference/           # Synthetic reference + indexes
├── results/
│   ├── quality_control/     # FastQC reports
│   ├── alignment/           # BAM files
│   └── variants/            # VCF files
└── logs/                    # Pipeline logs
```

## Development

To add new test cases or profiles:

1. Edit `generate_fastq.py` to add quality profiles
2. Edit `run_smoke.sh` to add validation checks
3. Run `test_smoke_components.sh` to verify generators
