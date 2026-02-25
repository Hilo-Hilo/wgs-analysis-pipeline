# Troubleshooting

## Quick Diagnostics

```bash
./scripts/check_requirements.sh          # verify tools + system
conda activate wgs_analysis && conda list | grep -E "(fastqc|fastp|bwa|samtools|bcftools)"
df -h                                    # disk space
tail logs/*.log                          # recent errors
```

## Common Errors

### "Command not found" (fastqc, bwa, etc.)

Conda environment not activated, or tools not installed.

```bash
conda activate wgs_analysis
# If env missing:
conda create -n wgs_analysis -c bioconda -c conda-forge \
    python=3.11 fastqc fastp bwa samtools bcftools -y
```

### "No FASTQ files found"

Pipeline expects `*_R1.fastq.gz` / `*_R2.fastq.gz` (or `*_1.fq.gz` / `*_2.fq.gz`).

```bash
ls data/raw/           # check actual filenames
# generate test data if needed:
python3 tests/generate_sample_data.py --output-dir data/raw
```

### Docker build fails on conda channel/network errors

The Dockerfile now retries `conda/mamba` environment creation with exponential backoff (4 attempts). If your network is unstable, retrying the build usually succeeds:

```bash
docker build --build-arg MINICONDA_INSTALLER=Miniconda3-latest-Linux-aarch64.sh -t wgs-pipeline:latest .
```

If failures persist, verify DNS/network egress to `conda.anaconda.org`.

Note: the base Docker image intentionally excludes `vep` to avoid pulling legacy/incompatible toolchains on ARM. `vep` is optional and only needed for the annotation step.

### GPU alignment fails on DGX

GPU mode requires both NVIDIA runtime and Parabricks:

```bash
nvidia-smi
which pbrun
```

Run alignment in GPU mode explicitly:

```bash
./scripts/alignment.sh --use-gpu --gpu-aligner parabricks --gpu-count 1 --threads 32
```

If either `nvidia-smi` or `pbrun` is missing, alignment now fails fast with an actionable error and you can rerun in CPU mode by dropping `--use-gpu`.

### Out of memory / process killed

Reduce threads (each thread uses more RAM):

```bash
./scripts/alignment.sh --threads 2
```

Or use the laptop profile: `--profile config/profiles/laptop.conf`

### BWA index files missing

```bash
cd data/reference/GRCh38/
bwa index GRCh38_latest_genomic.fna
samtools faidx GRCh38_latest_genomic.fna
```

### Low mapping rate (<80%)

Check that your reference genome matches your sample species and build (GRCh38 for human). Run QC first to rule out contamination or adapter issues.

### Very few / too many variants

| Expected (30x WGS) | Concern |
|---------------------|---------|
| 4-5M variants       | Normal  |
| <1M variants        | Low coverage or wrong reference |
| >10M variants       | Contamination or lenient filters |

Adjust `VARIANT_MIN_DEPTH` and `VARIANT_MIN_QUALITY` in your config.

## VEP Annotation Errors

### Exit Codes Reference

| Code | Meaning | Action |
|------|---------|--------|
| 10 | VEP not installed | `conda install -c bioconda ensembl-vep` |
| 11 | bcftools not installed | `conda install -c bioconda bcftools` |
| 12 | Input VCF not found | Check file path, run variant calling first |
| 13 | Malformed VCF | See "Malformed VCF" section below |
| 14 | Empty VCF (no variants) | Use `--skip-empty-check` if expected |
| 15 | Invalid VEP cache | See "VEP cache issues" section below |
| 16 | Insufficient memory | Use `--buffer-size 500 --threads 1` |
| 17 | VEP execution failed | Check `logs/annotation.log` for details |
| 18 | Post-processing failed | Annotated VCF may be incomplete |

### VEP not found

VEP (Ensembl Variant Effect Predictor) is optional for the core pipeline but required for annotation.

```bash
# Install via conda
conda install -c bioconda ensembl-vep

# Or via Docker
docker pull ensemblorg/ensembl-vep

# Verify installation
vep --help
```

### VEP cache issues

VEP requires a local cache for offline annotation:

```bash
# Download GRCh38 cache (requires ~15GB disk space)
vep_install -a cf -s homo_sapiens -y GRCh38 -c /path/to/vep_cache

# Check cache directory
ls -la /path/to/vep_cache/homo_sapiens/

# Validate before running
./scripts/vep_annotation.sh --validate-only --cache-dir /path/to/vep_cache
```

### Malformed VCF errors

If VCF validation fails:

```bash
# Check VCF format
bcftools view -H input.vcf.gz | head

# Validate VCF
bcftools view -h input.vcf.gz

# Repair common issues
bcftools norm -c ws input.vcf.gz -O z -o fixed.vcf.gz
```

### Empty VCF (no variants)

All variants may have been filtered out during variant calling:

```bash
# Check variant count
bcftools view -H results/variants/*_filtered.vcf.gz | wc -l

# If intentional, proceed with:
./scripts/vep_annotation.sh --skip-empty-check
```

### Out of memory during annotation

VEP can be memory-intensive. For 16GB systems:

```bash
# Reduce buffer size and threads
./scripts/vep_annotation.sh --buffer-size 500 --threads 1 --min-memory-gb 4

# Monitor memory during annotation
watch -n 5 'free -h'
```

### Resuming interrupted annotation

If annotation is interrupted, use resume mode:

```bash
./scripts/vep_annotation.sh --resume

# Or force a fresh start
./scripts/vep_annotation.sh --force
```

### Dry run for validation

Before running annotation, validate setup:

```bash
# Check all prerequisites without running
./scripts/vep_annotation.sh --validate-only --verbose

# Preview what would happen
./scripts/vep_annotation.sh --dry-run --verbose
```

## Quality Benchmarks

| Metric | Good | Investigate |
|--------|------|-------------|
| Mean quality score | >30 | <25 |
| Adapter contamination | <5% | >10% |
| Duplication rate | <20% | >30% |
| Mapping rate | >85% | <75% |
| Ti/Tv ratio | 2.0-2.1 | <1.8 or >2.5 |

## Emergency Reset

```bash
rm -rf results/ logs/ temp/
mkdir -p results logs temp
./scripts/check_requirements.sh
```
