# WGS Analysis Environment Setup - 16GB Optimized

## Overview
This document provides detailed instructions for setting up the conda environment for whole genome sequencing (WGS) analysis specifically optimized for 16GB RAM systems.

## System Requirements (16GB Optimized)
- **RAM**: 16GB (optimized configuration)
- **Available Storage**: 400GB+ for complete WGS analysis
- **CPU**: 4+ cores recommended
- **Operating System**: Linux, macOS, or Windows with WSL
- **Conda/Miniconda**: Version 4.0 or higher

## Environment Configuration

### Base Environment
- **Environment Name**: `wgs_analysis`
- **Python Version**: 3.11.13
- **R Version**: 4.3.3
- **Target Reference**: GRCh38 (hg38)
- **Architecture**: ARM64 optimized

### Installation Commands
```bash
# Create environment with base packages
conda create -n wgs_analysis python=3.11 r-base=4.3 -c conda-forge -y

# Activate environment
conda activate wgs_analysis

# Install bioinformatics tools (16GB optimized selection)
conda install -c bioconda fastp fastqc bwa samtools bcftools ensembl-vep -y
```

## Installed Software Stack

### Core Bioinformatics Tools (16GB Optimized)
| Tool | Version | Purpose | Channel | Memory Usage |
|------|---------|---------|---------|--------------|
| **fastp** | 1.0.1 | High-performance FASTQ preprocessing | bioconda | ~4GB |
| **fastqc** | 0.12.1 | Quality control for sequencing data | bioconda | ~2GB |
| **BWA** | 0.7.19 | Memory-efficient aligner for mapping | bioconda | ~10GB |
| **samtools** | 1.22.1 | SAM/BAM file manipulation | bioconda | ~3GB |
| **bcftools** | 1.22 | VCF file manipulation and variant calling | bioconda | ~3GB |
| **ensembl-vep** | Latest | Variant effect prediction and annotation | bioconda | ~6GB |
| **htslib** | 1.22.1 | High-throughput sequencing data formats | bioconda | N/A |

### Runtime Dependencies
| Tool | Version | Purpose |
|------|---------|---------|
| **OpenJDK** | 23.0.2 | Java runtime for FastQC |
| **Perl** | 5.32.1 | Perl runtime for FastQC |
| **ISA-L** | 2.31.1 | Compression library for fastp |

### Programming Languages
| Language | Version | Use Case |
|----------|---------|----------|
| **Python** | 3.11.13 | Scripting and data analysis |
| **R** | 4.3.3 | Statistical analysis and visualization |

## Project Directory Structure
```
wgs_analysis/
├── data/
│   ├── raw/              # Raw FASTQ files
│   ├── processed/        # Cleaned/processed FASTQ files
│   └── reference/        # Reference genome files
├── documentation/        # Project documentation
│   ├── README.md            # Project overview
│   ├── projectplan.md       # Detailed project plan
│   ├── environment_setup.md # This file
│   └── analysis_workflow.md # Analysis procedures
├── logs/                 # Analysis log files
├── results/              # Analysis outputs
│   ├── fastqc_raw/          # Raw data quality reports
│   ├── fastqc_clean/        # Cleaned data quality reports
│   ├── alignment/           # Mapping results
│   └── variants/            # Variant calling results
└── scripts/              # Custom analysis scripts
```

## Environment Verification

### Tool Verification Commands
```bash
# Activate environment
conda activate wgs_analysis

# Verify each tool installation
fastp --version
fastqc --version
bwa
samtools --version
bcftools --version
python --version
R --version
```

### Expected Output Examples
```
fastp 1.0.1
FastQC v0.12.1
BWA-0.7.19
samtools 1.22.1
bcftools 1.22
Python 3.11.13
R version 4.3.3
```

## Environment Management

### Activation/Deactivation
```bash
# Activate environment
conda activate wgs_analysis

# Deactivate environment
conda deactivate
```

### Environment Information
```bash
# List installed packages
conda list

# Export environment specification
conda env export > environment.yml

# Environment information
conda info --envs
```

## Troubleshooting

### Common Issues
1. **Package Conflicts**: Use `conda clean --all` to clear cache
2. **Architecture Mismatch**: Ensure ARM64 compatible packages
3. **Slow Installation**: Use `mamba` instead of conda for faster installs

### 16GB Performance Optimization
- **Memory Management**: Sequential processing to avoid memory conflicts
- **Thread Limits**: Conservative thread counts (4 max) to prevent memory overload
- **Storage**: SSD recommended for 400GB+ analysis requirements
- **Swap Space**: Not recommended - can cause performance issues

## Security Considerations
- **Data Privacy**: Ensure personal genomic data is properly secured
- **Access Control**: Limit environment access to authorized users
- **Backup**: Regular backup of analysis results and configurations

## Next Steps
1. **Reference Preparation**: Download and index GRCh38 reference genome
2. **Data Import**: Place raw FASTQ files in appropriate directories
3. **Quality Assessment**: Begin with FastQC analysis
4. **Pipeline Execution**: Follow analysis workflow documentation

## Maintenance
- **Updates**: Regularly update tools using `conda update`
- **Cleanup**: Remove temporary files and logs periodically
- **Monitoring**: Track disk usage and computational resources

## Additional Resources
- **Bioconda**: https://bioconda.github.io/
- **BWA Manual**: http://bio-bwa.sourceforge.net/bwa.shtml
- **samtools Documentation**: http://www.htslib.org/
- **FastQC Documentation**: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/