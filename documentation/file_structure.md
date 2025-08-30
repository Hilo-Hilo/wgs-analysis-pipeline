# WGS Analysis File Structure

## Overview
This document describes the complete file and directory structure for the WGS analysis project after completing CHM13 T2T mapping, including purpose and expected contents of each directory.

## Root Directory Structure
```
wgs_analysis/
├── data/                    # Data files and reference materials
├── documentation/           # Project documentation
├── logs/                    # Analysis log files
├── results/                 # Analysis outputs and reports
├── scripts/                 # Custom analysis scripts
├── README.md                # Main project overview
└── CLAUDE.md                # Project-specific instructions
```

## Detailed Directory Contents

### `/data/` - Data Storage
Main directory for all input, processed, and reference data files.

```
data/
├── raw/                     # Raw sequencing data (original)
│   ├── SAMPLE001_L01_UDB-406_1.fq.gz  # Forward reads (54GB)
│   └── SAMPLE001_L01_UDB-406_2.fq.gz  # Reverse reads (77GB)
├── processed/               # Cleaned and processed data
│   ├── SAMPLE001_clean_R1.fq.gz       # Cleaned forward reads (41GB)
│   └── SAMPLE001_clean_R2.fq.gz       # Cleaned reverse reads (56GB)
└── reference/               # Reference genome files
    ├── CHM13/                           # CHM13 T2T v2.0 (PRIMARY)
    │   ├── chm13v2.0.fa                 # Reference sequence (2.9GB)
    │   ├── chm13v2.0.fa.fai             # FASTA index
    │   ├── chm13v2.0.fa.0123            # BWA-MEM2 index files
    │   ├── chm13v2.0.fa.amb
    │   ├── chm13v2.0.fa.ann
    │   ├── chm13v2.0.fa.bwt.2bit.64
    │   └── chm13v2.0.fa.pac
    └── GRCh38/                          # GRCh38 (DEPRECATED)
        ├── GRCh38_latest_genomic.fna    # Reference sequence
        ├── GRCh38_latest_genomic.fna.fai # FASTA index
        └── GRCh38_latest_genomic.fna.*  # BWA index files
```

### `/documentation/` - Project Documentation
Comprehensive project documentation and guides.

```
documentation/
├── README.md                        # Documentation index and navigation
├── projectplan.md                   # Project plan, phases, and progress tracking
├── environment_setup.md             # Environment configuration guide
├── analysis_workflow.md             # Step-by-step analysis procedures
├── file_structure.md                # This file - project structure guide
├── chm13_mapping_results.md         # Comprehensive CHM13 mapping results
├── variant_calling_guide.md         # Complete variant calling guide
├── quality_control_guide.md         # Quality control procedures
├── reference_genome_comparison.md   # CHM13 vs GRCh38 comparison
└── pipeline_parameters.md           # Optimized pipeline parameters
```

### `/logs/` - Analysis Logs
Log files from various analysis steps for troubleshooting and monitoring.

```
logs/
├── bam_indexing.log                 # BAM indexing process logs
├── bwa_mem2_alignment.log           # BWA-MEM2 alignment logs
├── chm13_mapping_optimized.log      # Optimized CHM13 mapping logs
├── data_cleaning.log                # Data cleaning process logs
├── fastp_cleaning.log               # fastp cleaning detailed logs
├── flagstat.log                     # BAM flagstat logs
├── idxstats.log                     # BAM index statistics logs
└── samtools_stats.log               # Detailed samtools statistics logs
```

### `/results/` - Analysis Results
Organized analysis outputs and reports by analysis stage.

```
results/
├── alignment/                       # Mapping and alignment results
│   ├── SAMPLE001_chm13_sorted.bam        # Final sorted alignment (73GB)
│   ├── SAMPLE001_chm13_sorted.bam.bai    # BAM index (8.9MB)
│   ├── SAMPLE001_chm13_stats.txt         # Basic alignment statistics
│   ├── SAMPLE001_chm13_detailed_stats.txt # Comprehensive statistics
│   ├── SAMPLE001_chm13_idxstats.txt      # Per-chromosome statistics
│   └── chm13_mapping_optimized_summary.txt # Mapping summary report
├── fastqc_raw/                      # Raw data quality reports
│   ├── SAMPLE001_L01_UDB-406_1_fastqc.html
│   ├── SAMPLE001_L01_UDB-406_1_fastqc.zip
│   ├── SAMPLE001_L01_UDB-406_1_fastqc/
│   │   ├── Icons/                       # FastQC report icons
│   │   ├── Images/                      # Quality metric plots
│   │   ├── fastqc_data.txt             # Raw quality data
│   │   ├── fastqc_report.html          # Interactive quality report
│   │   └── summary.txt                 # Quality summary
│   ├── SAMPLE001_L01_UDB-406_2_fastqc.html
│   ├── SAMPLE001_L01_UDB-406_2_fastqc.zip
│   ├── SAMPLE001_L01_UDB-406_2_fastqc/
│   │   └── [similar structure to R1]
│   └── quality_control_summary.txt     # Overall QC summary
├── fastqc_clean/                    # Cleaned data quality reports
│   └── [Post-cleaning FastQC reports - ready for generation]
├── quality_control/                 # Quality control and cleaning reports
│   ├── cleaning_summary.txt             # Data cleaning summary
│   ├── fastp_report.html               # Interactive cleaning report
│   └── fastp_report.json               # Machine-readable cleaning data
├── variants/                        # Variant calling results (ready for use)
│   └── [Variant files will be generated here]
└── reports/                         # Additional analysis reports
    └── [Future comprehensive reports]
```

### `/scripts/` - Custom Scripts
Custom analysis scripts and utilities for pipeline execution.

```
scripts/
├── quality_control.sh               # Quality control automation
├── data_cleaning.sh                 # Data cleaning pipeline
├── chm13_mapping.sh                 # Basic CHM13 mapping script
└── chm13_mapping_optimized.sh       # Optimized mapping pipeline
```

## File Naming Conventions

### Sample Files
- **Sample ID**: SAMPLE001 (consistent across all files)
- **Raw Data Pattern**: `{sample_id}_L01_UDB-406_{R1|R2}.fq.gz`
- **Cleaned Data Pattern**: `{sample_id}_clean_{R1|R2}.fq.gz`

### Analysis Results
- **Alignment Files**: `{sample_id}_chm13_{processing_stage}.{extension}`
- **Statistics Files**: `{sample_id}_chm13_{metric_type}.txt`
- **Report Files**: `{tool}_{report_type}.{html|json|txt}`

### Quality Control Reports
- **FastQC Raw**: `{sample_id}_L01_UDB-406_{read}_fastqc.{html|zip}`
- **FastQC Clean**: `{sample_id}_clean_{read}_fastqc.{html|zip}`
- **Processing Reports**: `{tool}_report.{html|json}`

## Current File Sizes and Storage

### Data Directory (203 GB total)
- **Raw FASTQ files**: 131 GB (54GB + 77GB)
- **Cleaned FASTQ files**: 97 GB (41GB + 56GB)
- **CHM13 reference**: 2.9 GB
- **GRCh38 reference**: 3.2 GB (deprecated, can be removed)

### Results Directory (73+ GB total)
- **Final BAM file**: 73 GB
- **BAM index**: 8.9 MB
- **Quality reports**: ~100 MB
- **Statistics files**: ~1 MB total

### Storage Optimization Completed
- ✅ **Temporary files removed**: Saved 1.5 GB (BAM temp file)
- ✅ **Direct BAM pipeline**: Avoided 330+ GB intermediate SAM
- ✅ **Compressed outputs**: All major files compressed

## Data Flow Summary

### Completed Pipeline
```
Raw FASTQ (131GB)
    ↓ [fastp cleaning]
Cleaned FASTQ (97GB) 
    ↓ [BWA-MEM2 + CHM13]
Sorted BAM (73GB)
    ↓ [samtools stats]
Quality Metrics & Statistics
```

### Ready for Next Steps
```
Sorted BAM + CHM13 Reference
    ↓ [bcftools/GATK]
Variant Calls (VCF)
    ↓ [annotation]
Annotated Variants
    ↓ [analysis]
Genomic Reports
```

## Quality Assessment Results

### Data Quality Metrics ✅
- **Raw data quality**: Excellent (Q39 average, 0 poor sequences)
- **Cleaning effectiveness**: 26% size reduction, quality preserved
- **Mapping quality**: 88.62% mapping rate to CHM13
- **Alignment quality**: 87.89% properly paired, 0.364% error rate

### File Integrity Status ✅
- **All critical files present**: BAM, index, statistics
- **No missing dependencies**: All required files available
- **Storage optimized**: Temporary files cleaned
- **Documentation complete**: All analyses documented

## Security and Backup Considerations

### File Permissions
```bash
# Data files (sensitive genomic data)
chmod 600 data/raw/*.fq.gz
chmod 600 data/processed/*.fq.gz
chmod 600 results/alignment/*.bam

# Scripts (executable)
chmod 755 scripts/*.sh

# Documentation (readable)
chmod 644 documentation/*.md
```

### Critical Files for Backup
1. **Raw data**: `data/raw/SAMPLE001_*.fq.gz` (131GB)
2. **Final alignment**: `results/alignment/SAMPLE001_chm13_sorted.bam*` (73GB)
3. **Quality reports**: `results/fastqc_*/` and `results/quality_control/`
4. **Documentation**: `documentation/` directory
5. **Analysis logs**: `logs/` directory

### Recommended Cleanup (Optional)
- **GRCh38 reference** (3.2GB): Can be removed as CHM13 is now primary
- **Intermediate statistics**: Some verbose log files can be archived

## Maintenance Tasks

### Regular Monitoring
- **Disk space**: Monitor `/results/` directory growth
- **File integrity**: Verify BAM file checksums periodically
- **Log rotation**: Archive old log files monthly

### Future Expansion
- **variants/**: Will grow with variant calling results
- **reports/**: Will contain comprehensive genomic analysis reports
- **Additional references**: May add population-specific references

## Troubleshooting File Issues

### Common Problems
1. **BAM file corruption**: Verify with `samtools quickcheck`
2. **Missing index files**: Regenerate with `samtools index`
3. **Permission errors**: Check file ownership and permissions
4. **Disk space**: Monitor available space regularly

### File Validation Commands
```bash
# Check BAM file integrity
samtools quickcheck results/alignment/SAMPLE001_chm13_sorted.bam

# Verify index exists and is current
ls -la results/alignment/SAMPLE001_chm13_sorted.bam.bai

# Check disk usage
du -sh results/
df -h .
```

## Project Status

### ✅ Completed Phases
- Environment setup and tool installation
- Raw data quality assessment
- Data cleaning and preprocessing
- CHM13 T2T reference genome setup
- Read mapping and alignment
- Post-mapping quality control
- File organization and cleanup

### 🎯 Ready for Execution
- Variant calling with bcftools or GATK
- Coverage analysis across chromosomes
- Structural variant detection
- Comprehensive genomic analysis

---

**File structure optimized**: July 20, 2025  
**Total project size**: ~276 GB (optimized from 330+ GB)  
**Status**: Production-ready for variant calling