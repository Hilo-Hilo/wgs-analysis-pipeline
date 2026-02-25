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
conda create -n wgs_analysis -c bioconda -c conda-forge \
    python=3.11 \
    "fastqc>=0.11.9" \
    "fastp>=0.20.0" \
    "bwa-mem2>=2.2.1" \
    "samtools>=1.15" \
    "bcftools>=1.15" \
    -y
conda activate wgs_analysis

# Validate versions are within supported ranges
./scripts/validate_deps.sh

# Optional (annotation step only)
# conda install -n wgs_analysis -c bioconda ensembl-vep
```

See [Dependency Policy](../../DEPENDENCIES.md) for supported ranges.

## Prepare Data

Place paired-end FASTQ files in `data/raw/`:

| File pattern | Example |
|---|---|
| `*_R1.fastq.gz` or `*_1.fq.gz` | `MySample_R1.fastq.gz` |
| `*_R2.fastq.gz` or `*_2.fq.gz` | `MySample_R2.fastq.gz` |

## Run Pipeline

### Full run

```bash
bash run_pipeline.sh --sample-id MySample --input-dir data/raw --output-dir results
```

### Step-by-step

```bash
conda activate wgs_analysis
./scripts/quality_control.sh
./scripts/data_cleaning.sh
./scripts/alignment.sh
./scripts/variant_calling.sh
./scripts/vep_annotation.sh
```

### Optional GPU alignment (Parabricks)

```bash
./scripts/alignment.sh \
  --input-dir data/processed \
  --reference data/reference/GRCh38/GRCh38_latest_genomic.fna \
  --output-dir results/alignment \
  --sample-id MySample \
  --threads 32 \
  --use-gpu --gpu-aligner parabricks --gpu-count 1
```

## Reference Genome Setup

```bash
mkdir -p data/reference/GRCh38 && cd data/reference/GRCh38
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip *.fna.gz
mv *.fna GRCh38_latest_genomic.fna
samtools faidx GRCh38_latest_genomic.fna
bwa-mem2 index GRCh38_latest_genomic.fna
cd ../../..
```

## Output Layout

```text
results/
  quality_control/   # FastQC reports
  alignment/         # BAM + mapping stats
  variants/          # Raw + filtered VCF
  annotation/        # VEP outputs
logs/                # Per-step logs
```

## Related Docs

- [Troubleshooting](troubleshooting.md)
- [Sample Registry](../reference/sample-registry.md)
- [Codebase Overview](../architecture/codebase-overview.md)
