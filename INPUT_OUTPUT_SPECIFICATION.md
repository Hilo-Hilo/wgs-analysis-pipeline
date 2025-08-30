# WGS Analysis Pipeline - Input/Output Specification

## Overview

This pipeline transforms raw whole genome sequencing (WGS) data into annotated variants ready for analysis. It is optimized for local machines with 16GB RAM.

## Required Inputs

### Primary Input Files
- **Raw FASTQ Files** (paired-end, Illumina)
  - Forward reads: `*_R1.fastq.gz` or `*_1.fastq.gz`
  - Reverse reads: `*_R2.fastq.gz` or `*_2.fastq.gz`
  - Expected size: 50-100GB total for 30x WGS
  - Format: Gzip-compressed FASTQ

### Reference Genome
- **Human Reference Genome**:
  - GRCh38/hg38: 3.2GB download
- **BWA Indexes**: Generated automatically (adds ~3GB)
- **samtools Index**: Generated automatically (adds ~200MB)

### System Requirements
- **RAM**: 16GB (pipeline optimized for this constraint)
- **CPU**: 4+ cores (8+ recommended)
- **Storage**: 400GB free space minimum
- **OS**: Linux or macOS

## Complete Output Structure

```
results/
├── quality_control/
│   ├── fastqc_reports/
│   │   ├── sample_R1_fastqc.html           # Quality control report (web viewable)
│   │   ├── sample_R1_fastqc.zip            # Detailed QC data
│   │   ├── sample_R2_fastqc.html           # Quality control report (web viewable)
│   │   └── sample_R2_fastqc.zip            # Detailed QC data
│   └── quality_summary.txt                 # Overall quality metrics
│
├── cleaned_reads/
│   ├── sample_clean_R1.fastq.gz            # Adapter-trimmed forward reads
│   ├── sample_clean_R2.fastq.gz            # Adapter-trimmed reverse reads
│   └── cleaning_report.txt                 # Trimming statistics
│
├── alignment/
│   ├── sample_sorted.bam                   # Primary alignment file (80-120GB)
│   ├── sample_sorted.bam.bai               # BAM index file
│   ├── alignment_stats.txt                 # Mapping statistics
│   └── coverage_summary.txt                # Coverage depth metrics
│
├── variants/
│   ├── sample_raw.vcf.gz                   # Raw variant calls (1-5GB)
│   ├── sample_raw.vcf.gz.tbi               # VCF index
│   ├── sample_filtered.vcf.gz              # Quality-filtered variants
│   ├── sample_filtered.vcf.gz.tbi          # Filtered VCF index
│   └── variant_stats.txt                   # Variant calling statistics
│
└── annotation/
    ├── sample_annotated.vcf.gz              # Fully annotated variants
    ├── sample_annotated.vcf.gz.tbi          # Annotated VCF index
    ├── sample_high_impact.txt               # High-impact variants only
    ├── sample_clinvar.txt                   # ClinVar matches
    └── annotation_summary.txt               # Annotation statistics

logs/
├── quality_control.log                     # QC process logs
├── alignment.log                           # Alignment process logs
├── variant_calling.log                     # Variant calling logs
└── annotation.log                          # Annotation logs
```

## Expected File Sizes (30x WGS)

| File Type | Size Range | Description |
|-----------|------------|-------------|
| Raw FASTQ | 50-100GB | Input sequencing files |
| Cleaned FASTQ | 40-90GB | Adapter-trimmed reads |
| BAM file | 80-120GB | Aligned reads (primary output) |
| Raw VCF | 1-3GB | All variant calls |
| Filtered VCF | 500MB-1.5GB | Quality-filtered variants |
| Annotated VCF | 2-8GB | Variants with annotations |
| Logs | 10-50MB | Process logs |
| **Total Storage** | **200-350GB** | Complete analysis |

## Key Output Files for Analysis

### For Quality Assessment
- `results/quality_control/sample_R1_fastqc.html` - View in web browser
- `results/quality_control/sample_R2_fastqc.html` - View in web browser
- `results/alignment/alignment_stats.txt` - Check mapping rate

### For Variant Analysis
- `results/annotation/sample_annotated.vcf.gz` - Main analysis file
- `results/annotation/sample_high_impact.txt` - Clinically relevant variants
- `results/annotation/sample_clinvar.txt` - Known pathogenic variants

### For Further Processing
- `results/alignment/sample_sorted.bam` - For structural variants, CNVs
- `results/variants/sample_filtered.vcf.gz` - For custom analysis pipelines

## Memory-Optimized Processing Strategy

### Sequential Processing (16GB RAM)
1. **Quality Control**: 2GB peak memory
2. **Read Cleaning**: 4GB peak memory
3. **Alignment**: 12GB peak memory (optimized BWA settings)
4. **Variant Calling**: 6GB peak memory
5. **Annotation**: 8GB peak memory

### Automatic Cleanup
- Intermediate files removed after each step
- Only essential files retained
- Configurable cleanup behavior

## Expected Analysis Times (16GB RAM, 8 cores)

| Step | Time Range | Progress Indicators |
|------|------------|-------------------|
| Quality Control | 10-30 minutes | File-by-file progress |
| Read Cleaning | 30-90 minutes | Percentage complete |
| Alignment | 4-8 hours | Read pairs processed |
| Variant Calling | 2-4 hours | Chromosome progress |
| Annotation | 1-2 hours | Variant progress |
| **Total** | **8-15 hours** | Overall pipeline status |

## Validation Outputs

### Quality Metrics to Check
- **Mapping Rate**: Should be >85%
- **Mean Quality Score**: Should be >30
- **Coverage Depth**: Should be 25-35x for 30x sequencing
- **Variant Count**: 4-5 million for whole genome

### Warning Signs
- Mapping rate <80%: Possible contamination or wrong reference
- Quality score <25: Poor sequencing quality
- <3 million variants: Possible coverage issues
- >6 million variants: Possible contamination

## Resource Monitoring

The pipeline automatically monitors:
- Memory usage (prevents system crashes)
- Disk space (stops if <50GB remaining)
- CPU temperature (reduces threads if overheating)
- Process completion status

## Data Privacy and Security

### What's Stored
- Processed genomic data (alignments, variants)
- Quality metrics and statistics
- Analysis logs

### What's NOT Stored
- Raw personal identifiers
- Unencrypted sensitive data
- Cloud credentials or API keys

### Recommended Practices
- Use encrypted storage for sensitive data
- Regularly backup important results
- Follow institutional data policies
- Consider anonymizing sample IDs