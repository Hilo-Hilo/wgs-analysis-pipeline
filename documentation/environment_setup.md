# WGS Analysis Environment Setup

## Overview
This document provides detailed instructions for setting up the conda environment for personal whole genome sequencing (WGS) analysis on macOS ARM64 (Apple Silicon) systems.

## System Requirements
- **Operating System**: macOS (ARM64/Apple Silicon)
- **Conda/Miniconda**: Version 4.0 or higher
- **Available Storage**: Minimum 100GB for reference genome and analysis
- **RAM**: Minimum 16GB recommended for WGS analysis

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

# Install bioinformatics tools
conda install -c bioconda fastp fastqc bwa samtools bcftools -y
```

## Installed Software Stack

### Core Bioinformatics Tools
| Tool | Version | Purpose | Channel |
|------|---------|---------|---------|
| **fastp** | 1.0.1 | High-performance FASTQ preprocessing | bioconda |
| **fastqc** | 0.12.1 | Quality control for sequencing data | bioconda |
| **BWA** | 0.7.19 | Burrows-Wheeler Aligner for mapping | bioconda |
| **samtools** | 1.22.1 | SAM/BAM file manipulation | bioconda |
| **bcftools** | 1.22 | VCF file manipulation and variant calling | bioconda |
| **htslib** | 1.22.1 | High-throughput sequencing data formats | bioconda |

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

### Performance Optimization
- **Memory**: Increase available RAM for large genome analysis
- **Storage**: Use SSD for faster I/O operations
- **CPU**: Utilize multi-threading options in tools

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