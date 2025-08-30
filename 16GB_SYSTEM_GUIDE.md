# WGS Analysis on 16GB RAM Systems - Complete Guide

## What This Pipeline Does

**INPUT**: Paired-end FASTQ files (raw sequencing data, 50-100GB)
**OUTPUT**: Annotated genetic variants ready for analysis

See [INPUT_OUTPUT_SPECIFICATION.md](INPUT_OUTPUT_SPECIFICATION.md) for complete file list and sizes.

## System Requirements (16GB Optimized)

- **RAM**: 16GB (pipeline uses max 14GB, leaving 2GB for system)
- **CPU**: 4+ cores (optimized for 4-core processing)
- **Storage**: 400GB free space minimum
- **OS**: Linux or macOS
- **Time**: 8-15 hours for complete analysis

## Quick Setup

```bash
# 1. Check your system
./scripts/check_requirements.sh --min-ram 16 --min-disk 400

# 2. Install tools
conda create -n wgs_analysis -c bioconda -c conda-forge \
    python=3.9 fastqc fastp bwa samtools bcftools vep
conda activate wgs_analysis

# 3. Use 16GB configuration
source config/local_16gb.conf
```

## Complete Analysis Workflow

### Step 1: Prepare Your Data
```bash
# Create directories
mkdir -p data/raw data/reference results logs

# Copy your FASTQ files
cp your_sample_R1.fastq.gz data/raw/
cp your_sample_R2.fastq.gz data/raw/
```

### Step 2: Download Reference Genome
```bash
# Download GRCh38 (recommended for 16GB systems)
cd data/reference
mkdir GRCh38 && cd GRCh38
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip *.fna.gz
mv *.fna GRCh38_latest_genomic.fna

# Index reference (30-60 minutes)
samtools faidx GRCh38_latest_genomic.fna
bwa index GRCh38_latest_genomic.fna
cd ../../..
```

### Step 3: Run Analysis Pipeline
```bash
# Load 16GB optimized settings
source config/local_16gb.conf

# Run each step sequentially (prevents memory issues)
./scripts/quality_control.sh --threads 4              # 10-30 minutes
./scripts/data_cleaning.sh --threads 4                # 30-90 minutes
./scripts/alignment.sh --threads 4                    # 4-8 hours
./scripts/variant_calling.sh --threads 2              # 2-4 hours
./scripts/vep_annotation.sh --threads 2               # 1-2 hours

# Total time: 8-15 hours
```

## Expected Output Files

After successful completion, you will have:

```
results/
├── quality_control/
│   ├── sample_R1_fastqc.html          # View in browser
│   ├── sample_R2_fastqc.html          # View in browser
│   └── quality_summary.txt            # Quality metrics
├── alignment/
│   ├── sample_sorted.bam              # Main alignment file (80-120GB)
│   ├── sample_sorted.bam.bai          # Index file
│   └── alignment_stats.txt            # Mapping statistics
├── variants/
│   ├── sample_filtered.vcf.gz         # Quality-filtered variants
│   └── variant_stats.txt              # Variant statistics
└── annotation/
    ├── sample_annotated.vcf.gz        # MAIN OUTPUT - annotated variants
    ├── sample_high_impact.txt         # Clinically important variants
    └── sample_clinvar.txt             # Known disease variants
```

## Key Files for Analysis

1. **sample_annotated.vcf.gz** - Complete annotated variants for downstream analysis
2. **sample_high_impact.txt** - High-impact variants (likely functional effects)
3. **sample_clinvar.txt** - Variants with known clinical significance
4. **alignment_stats.txt** - Check mapping rate (should be >85%)

## Memory Management Features

### Automatic Resource Control
- **Sequential processing** - Only one memory-intensive step at a time
- **Memory monitoring** - Scripts stop if RAM usage exceeds 90%
- **Automatic cleanup** - Intermediate files removed to save space
- **Progress tracking** - Know current step and estimated completion

### 16GB Optimized Settings
- **BWA alignment**: 10GB max memory (vs 16GB+ in default configs)
- **Thread limits**: 4 threads max (vs 8+ in high-memory configs)
- **VEP annotation**: 6GB max memory (vs 16GB+ in default configs)
- **Compression**: Maximum compression to reduce file sizes

## Troubleshooting 16GB Systems

### "Out of Memory" Errors
```bash
# Reduce threads even further
./scripts/quality_control.sh --threads 2
./scripts/alignment.sh --threads 2

# Check memory usage
free -h
htop
```

### "No Space Left" Errors
```bash
# Clean up intermediate files
rm -rf temp/ logs/*.log
rm -rf results/*/temp/

# Check space usage
df -h
du -sh results/*
```

### Slow Performance
```bash
# Use fewer threads if system is struggling
export THREADS=2
source config/local_16gb.conf

# Monitor system resources
top
iostat 1
```

## Quality Validation

### Expected Results (30x WGS)
- **Mapping rate**: >85%
- **Total variants**: 4-5 million
- **High-impact variants**: 300-500
- **Ti/Tv ratio**: 2.0-2.1

### Check Results
```bash
# Mapping rate
grep "mapped (" results/alignment/alignment_stats.txt

# Variant count
bcftools stats results/variants/sample_filtered.vcf.gz | grep "number of records"

# High-impact count
wc -l results/annotation/sample_high_impact.txt
```

## Performance Tips for 16GB Systems

1. **Close other applications** during analysis
2. **Use SSD storage** for better I/O performance
3. **Run overnight** for long alignment steps
4. **Monitor temperature** - reduce threads if CPU overheats
5. **Use swap space** if available (but will slow down analysis)

## Next Steps After Analysis

### For Research
1. Filter variants by gene lists or pathways
2. Compare with public databases (gnomAD, UK Biobank)
3. Perform association studies or pathway analysis

### For Clinical Applications
1. Focus on high-impact variants in disease genes
2. Check ClinVar annotations for known pathogenic variants
3. Consult with genetic counselors for interpretation

### For Method Development
1. Benchmark against known reference samples
2. Optimize parameters for specific use cases
3. Integrate with downstream analysis pipelines

## Data Management

### Storage Recommendations
- **Keep**: annotated.vcf.gz, high_impact.txt, clinvar.txt
- **Archive**: sorted.bam, alignment_stats.txt
- **Delete**: intermediate files, logs (after confirming success)

### Backup Strategy
```bash
# Compress and backup essential results
tar -czf sample_analysis_results.tar.gz \
    results/annotation/sample_annotated.vcf.gz \
    results/annotation/sample_high_impact.txt \
    results/annotation/sample_clinvar.txt \
    results/alignment/alignment_stats.txt

# Store compressed backup safely
```

This guide provides everything needed to successfully run WGS analysis on a 16GB RAM system. The optimized configurations ensure reliable processing while the clear input/output specifications make it easy to understand exactly what files are created.