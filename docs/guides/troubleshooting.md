# Troubleshooting

## Quick Diagnostics

```bash
./scripts/check_requirements.sh
./scripts/validate_deps.sh --ci-mode
conda activate wgs_analysis && conda list | grep -E "(fastqc|fastp|bwa-mem2|samtools|bcftools)"
df -h
tail logs/*.log
```

## Common Errors

### "Command not found" (fastqc, bwa-mem2, etc.)

```bash
conda activate wgs_analysis
conda install -n wgs_analysis -c bioconda -c conda-forge \
  python=3.11 fastqc fastp bwa-mem2 samtools bcftools -y
```

### "No FASTQ files found"

```bash
ls data/raw/
python3 tests/generate_sample_data.py --output-dir data/raw
```

### GPU alignment fails

```bash
nvidia-smi
which pbrun
./scripts/alignment.sh --use-gpu --gpu-aligner parabricks --gpu-count 1 --threads 32
```

If `nvidia-smi` or `pbrun` is missing, rerun in CPU mode (omit `--use-gpu`).

### CPU alignment killed / OOM

Lower thread count and keep sort memory conservative:

```bash
./scripts/alignment.sh --threads 2
```

### Reference index missing

```bash
cd data/reference/GRCh38/
bwa-mem2 index GRCh38_latest_genomic.fna
samtools faidx GRCh38_latest_genomic.fna
```

## VEP Annotation Issues

| Code | Meaning | Action |
|---|---|---|
| 10 | VEP not installed | `conda install -c bioconda ensembl-vep` |
| 11 | bcftools missing | `conda install -c bioconda bcftools` |
| 12 | Input VCF missing | verify variant-calling output path |
| 13 | Malformed VCF | validate with `bcftools view -h` |
| 14 | Empty VCF | use `--skip-empty-check` if expected |
| 15 | Invalid VEP cache | re-download cache |
| 16 | Insufficient memory | lower `--buffer-size` and `--threads` |
| 17 | VEP execution failed | inspect `logs/annotation.log` |
| 18 | Post-processing failed | inspect output + logs |

### VEP cache setup

```bash
vep_install -a cf -s homo_sapiens -y GRCh38 -c /path/to/vep_cache
./scripts/vep_annotation.sh --validate-only --cache-dir /path/to/vep_cache
```

## Emergency Reset

```bash
rm -rf results/ logs/ temp/
mkdir -p results logs temp
./scripts/check_requirements.sh
```
