# Getting Started

## Quick Start (5 minutes)

```bash
# 1. Check system readiness
./scripts/check_requirements.sh

# 2. Validate dependency versions (optional but recommended)
./scripts/validate_deps.sh --ci-mode

# 3. Generate test data
python3 tests/generate_sample_data.py --output-dir data/raw --num-reads 2000

# 4. Dry-run quality control
conda activate wgs_analysis
./scripts/quality_control.sh --input-dir data/raw --dry-run
```

## Prerequisites

- Linux or macOS
- [Conda](https://docs.conda.io/en/latest/miniconda.html) (Miniconda or Mambaforge)
- 16 GB RAM, 400 GB free disk, 4+ CPU cores

### Create the environment

```bash
# Using version-pinned dependencies (see DEPENDENCIES.md for version policy)
conda create -n wgs_analysis -c bioconda -c conda-forge \
    python=3.11 \
    "fastqc>=0.11.9" \
    "fastp>=0.20.0" \
    "bwa>=0.7.17" \
    "samtools>=1.15" \
    "bcftools>=1.15" \
    -y
conda activate wgs_analysis

# Validate versions are within supported ranges
./scripts/validate_deps.sh

# Optional (annotation step only)
# conda install -n wgs_analysis -c bioconda ensembl-vep
```

> **Note**: See [DEPENDENCIES.md](DEPENDENCIES.md) for the full version policy and supported ranges.

## Prepare Your Data

Place paired-end FASTQ files in `data/raw/`:

| File pattern | Example |
|---|---|
| `*_R1.fastq.gz` or `*_1.fq.gz` | `MySample_R1.fastq.gz` |
| `*_R2.fastq.gz` or `*_2.fq.gz` | `MySample_R2.fastq.gz` |

For 30x human WGS, expect ~50-100 GB of raw FASTQ input.

## Run the Pipeline

### Step by step

```bash
conda activate wgs_analysis

# 1. Quality control
./scripts/quality_control.sh

# 2. Adapter/quality trimming
./scripts/data_cleaning.sh

# 3. Alignment (needs reference genome -- see below)
./scripts/alignment.sh

# 4. Variant calling
./scripts/variant_calling.sh

# 5. Annotation
./scripts/vep_annotation.sh
```

### Or all at once

```bash
bash run_pipeline.sh --sample-id MySample --input-dir data/raw --output-dir results
```

Every script accepts `--help`, `--dry-run`, `--threads N`, and `--verbose`.

### DGX GPU alignment (optional)

```bash
./scripts/alignment.sh \
  --input-dir data/processed \
  --reference data/reference/GRCh38/GRCh38_latest_genomic.fna \
  --output-dir results/alignment \
  --sample-id MySample \
  --threads 32 \
  --use-gpu --gpu-aligner parabricks --gpu-count 1
```

## Reference Genome

Download and index GRCh38 (one-time, ~30 min):

```bash
mkdir -p data/reference/GRCh38 && cd data/reference/GRCh38
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip *.fna.gz
mv *.fna GRCh38_latest_genomic.fna
samtools faidx GRCh38_latest_genomic.fna
bwa index GRCh38_latest_genomic.fna
cd ../../..
```

## Output

```
results/
  quality_control/   # FastQC HTML reports
  alignment/         # BAM files + mapping stats
  variants/          # Raw + filtered VCF
  annotation/        # VEP-annotated variants
logs/                # Per-step log files
```

### Key quality benchmarks

| Metric | Good | Investigate |
|--------|------|-------------|
| Mean quality score | >30 | <25 |
| Mapping rate | >85% | <75% |
| Variant count (30x WGS) | 4-5M | <1M or >10M |
| Ti/Tv ratio | 2.0-2.1 | <1.8 or >2.5 |

## Resource Estimates (30x human WGS, 4 cores)

| Stage | Time | Peak RAM |
|-------|------|----------|
| QC | 15-30 min | 2 GB |
| Cleaning | 30-60 min | 2 GB |
| Alignment | 6-10 h | 14 GB |
| Variant calling | 2-4 h | 3 GB |
| Annotation | 1-2 h | 6 GB |

## Next Steps

- **Tune for your hardware**: see profiles in `config/profiles/`
- **Run in Docker**: `docker compose up --build`
- **Troubleshoot**: [TROUBLESHOOTING.md](TROUBLESHOOTING.md)
