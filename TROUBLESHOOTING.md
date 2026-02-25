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
