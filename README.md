# Whole Genome Sequencing Analysis Pipeline

A local-optimized, comprehensive pipeline for analyzing whole genome sequencing (WGS) data from raw FASTQ files to annotated variants. **Designed for 16GB RAM systems** and users with minimal bioinformatics experience but comfortable with command-line tools.

## What Goes In, What Comes Out

**INPUT**: Raw FASTQ files (paired-end, 50-100GB)  
**OUTPUT**: Annotated variants ready for analysis (see [INPUT_OUTPUT_SPECIFICATION.md](INPUT_OUTPUT_SPECIFICATION.md))

## Quick Navigation

| I want to... | Go to... |
|---------------|----------|
| **See exactly what files are created** | [INPUT_OUTPUT_SPECIFICATION.md](INPUT_OUTPUT_SPECIFICATION.md) |
| **Complete 16GB system setup guide** | [16GB_SYSTEM_GUIDE.md](16GB_SYSTEM_GUIDE.md) |
| **Get started quickly** | [GETTING_STARTED.md](GETTING_STARTED.md) |
| **Check if my system is ready** | `./scripts/check_requirements.sh --min-ram 16 --min-disk 400` |
| **Use 16GB RAM optimized settings** | Built-in (config/default.conf) |
| **Fix problems** | [TROUBLESHOOTING.md](TROUBLESHOOTING.md) |

**First time?** Run this test:
```bash
./scripts/check_requirements.sh --min-ram 16 --min-disk 400
```

## Overview

This pipeline processes paired-end Illumina WGS data through quality control, alignment, variant calling, and annotation stages. **Optimized for local machines with 16GB RAM**. Supports GRCh38 reference genome with memory-efficient processing.

### Key Features
- **16GB RAM optimized** - Runs efficiently on local workstations
- **Memory-aware processing** - Sequential steps prevent system crashes
- **Complete workflow** - Raw FASTQ to annotated variants
- **GRCh38 reference** - Standard human genome analysis
- **Automatic cleanup** - Manages disk space during analysis
- **Clear outputs** - Exactly specified file structure
- **Progress tracking** - Real-time status updates

## Pipeline Stages

1. **Quality Control** (10-30 minutes, 2GB peak RAM)
   - FastQC analysis of raw reads
   - Quality metrics and contamination detection
   - HTML reports for visual inspection

2. **Read Cleaning** (30-90 minutes, 4GB peak RAM)
   - Adapter trimming with fastp
   - Quality filtering and read length filtering
   - Statistics on removed sequences

3. **Alignment** (4-8 hours, 12GB peak RAM)
   - BWA-MEM mapping to GRCh38 reference
   - SAMtools sorting and indexing
   - Alignment statistics and coverage metrics

4. **Variant Calling** (2-4 hours, 6GB peak RAM)
   - BCFtools variant detection
   - Quality filtering and normalization
   - VCF indexing and statistics

5. **Annotation** (1-2 hours, 8GB peak RAM)
   - VEP variant annotation with gnomAD and ClinVar
   - High-impact variant identification
   - Clinical significance assessment

## Quick Start for 16GB Systems

### 1. Check System Requirements
```bash
# Verify your system meets 16GB requirements
./scripts/check_requirements.sh --min-ram 16 --min-disk 400
```

### 2. Set Up Environment
```bash
# Create conda environment
conda create -n wgs_analysis -c bioconda -c conda-forge \
    python=3.9 fastqc fastp bwa samtools bcftools vep

# Activate environment
conda activate wgs_analysis
```

### 3. Use 16GB Optimized Configuration
```bash
# Copy optimized configuration
cp config/local_16gb.conf config/my_analysis.conf

# Or load directly in scripts:
source config/local_16gb.conf
```

### 4. Run with Your Data
```bash
# 1. Place your FASTQ files
mkdir -p data/raw
cp your_sample_R1.fastq.gz data/raw/
cp your_sample_R2.fastq.gz data/raw/

# 2. Run complete pipeline (sequential, memory-safe)
./scripts/quality_control.sh --threads 4
./scripts/data_cleaning.sh --threads 4
./scripts/alignment.sh --threads 4
./scripts/variant_calling.sh --threads 2
./scripts/vep_annotation.sh --threads 2

# Total time: 8-15 hours
# Total storage needed: 300-400GB
```

## 16GB RAM Optimization Notes

### Memory Management Strategy
**Key Finding**: Original BWA works efficiently with 16GB RAM, but BWA-MEM2 requires 128GB+ for indexing.

**Our Solution**: Use standard BWA with optimized parameters
```bash
# Memory-efficient alignment (used automatically)
bwa mem -t 4 -M reference.fna reads_R1.fq reads_R2.fq
```

### Expected Performance (16GB system)
- **Input**: 50-100GB paired-end FASTQ files
- **Processing Time**: 8-15 hours on 4-8 cores
- **Peak Memory**: 12GB during alignment
- **Output**: 4-5 million variants (filtered)
- **Storage Required**: 300-400GB total

### Critical Success Factors
- **Sequential processing** prevents memory conflicts
- **Automatic cleanup** manages disk space
- **Progress monitoring** prevents system hangs
- **Optimized parameters** balance speed and memory

## Features for 16GB Systems

### Memory-Optimized Scripts
- **`config/local_16gb.conf`** - Pre-configured settings for 16GB RAM
- **`check_requirements.sh`** - Validates memory and disk requirements
- **Sequential processing** - Prevents memory conflicts between steps

### Smart Resource Management
- **Memory monitoring** - Stops processes if RAM usage gets dangerous
- **Automatic cleanup** - Removes intermediate files to save space
- **Progress tracking** - Shows current step and estimated completion
- **Error recovery** - Can resume from failed steps

### Input/Output Clarity
- **[INPUT_OUTPUT_SPECIFICATION.md](INPUT_OUTPUT_SPECIFICATION.md)** - Exact file list and sizes
- **Predictable structure** - All outputs in clearly organized directories
- **Size estimates** - Know storage requirements before starting

### Documentation for Local Use
- **[GETTING_STARTED.md](GETTING_STARTED.md)** - Complete setup guide
- **[TROUBLESHOOTING.md](TROUBLESHOOTING.md)** - Local system problems and solutions

## Repository Structure

```
wgs-analysis-pipeline/
├── INPUT_OUTPUT_SPECIFICATION.md  # WHAT FILES ARE CREATED
├── GETTING_STARTED.md              # Setup guide for beginners
├── TROUBLESHOOTING.md              # Local system problems
├── scripts/                        # Memory-optimized pipeline scripts
│   ├── check_requirements.sh       # 16GB RAM validation
│   ├── quality_control.sh          # FastQC analysis
│   ├── data_cleaning.sh            # Read trimming
│   ├── alignment.sh                # BWA alignment
│   └── load_config.sh              # Configuration management
├── config/                         # Memory configurations
│   ├── default.conf                # Standard settings
│   ├── local_16gb.conf             # 16GB RAM optimized
│   └── example.conf                # Template
├── analysis/                       # Post-processing tools
├── documentation/                  # Advanced guides
└── templates/                      # Configuration examples
```

## 🛠️ Script Usage Examples

All scripts now include comprehensive help documentation. Use `--help` with any script:

### Quality Control
```bash
# Basic usage
./scripts/quality_control.sh

# With custom options
./scripts/quality_control.sh --input-dir /path/to/fastq --threads 16 --verbose

# Dry run to see what will be processed
./scripts/quality_control.sh --dry-run

# Get help
./scripts/quality_control.sh --help
```

### System Validation
```bash
# Check if your system is ready for analysis
./scripts/check_requirements.sh

# Check with custom memory requirements
./scripts/check_requirements.sh --min-ram 32 --min-disk 1000

# Verbose output showing all checks
./scripts/check_requirements.sh --verbose
```

### Sample Data Management
```bash
# Download small test dataset
./scripts/download_sample_data.sh --type small

# Download medium dataset for comprehensive testing
./scripts/download_sample_data.sh --type medium

# List available datasets
./scripts/download_sample_data.sh --list
```

### Configuration Management
```bash
# Show current configuration
./scripts/load_config.sh show

# Create new configuration from template
./scripts/load_config.sh create my_analysis

# List available configurations
./scripts/load_config.sh list

# Validate configuration
./scripts/load_config.sh validate config/my_analysis.conf
```

## 🧬 Analysis Tools

### Pharmacogenomics Analysis
Analyze genetic variants affecting drug metabolism:

```python
python analysis/pharmacogenomics.py --vcf input.vcf.gz
```

### Comprehensive Analysis
Run multiple analysis modules:

```python
python analysis/comprehensive_analysis.py --vcf input.vcf.gz
```

## ☁️ Cloud Deployment

### Google Cloud Platform
Optimized configuration for WGS analysis:

```bash
# Create high-memory instance for alignment
gcloud compute instances create wgs-analysis \
    --machine-type=n2d-highmem-16 \
    --boot-disk-size=500GB \
    --zone=us-central1-a
```

Cost optimization strategies:
- Use preemptible instances for non-critical stages
- Downgrade to n2-standard-4 after alignment
- Store results in Cloud Storage

See [cloud_infrastructure_guide.md](documentation/cloud_infrastructure_guide.md) for detailed setup.

## 📚 Documentation

- [Environment Setup](documentation/environment_setup.md) - Tool installation and configuration
- [Complete Analysis Workflow](documentation/complete_analysis_workflow.md) - Step-by-step pipeline guide
- [Cloud Infrastructure Guide](documentation/cloud_infrastructure_guide.md) - GCP deployment
- [BWA-MEM2 Analysis](documentation/BWA_MEM2_ROOT_CAUSE_ANALYSIS.md) - Memory requirement investigation

## 🔧 Troubleshooting

**Having issues?** → Check [TROUBLESHOOTING.md](TROUBLESHOOTING.md) for comprehensive solutions

### Quick Fixes for Common Issues

1. **"Command not found" errors**
   ```bash
   conda activate wgs_analysis  # Activate environment first
   ./scripts/check_requirements.sh  # Verify installation
   ```

2. **Out of memory errors**
   ```bash
   ./scripts/quality_control.sh --threads 4  # Reduce threads
   ./scripts/check_requirements.sh --min-ram 16  # Check requirements
   ```

3. **No FASTQ files found**
   ```bash
   # Check file extensions (.fq.gz or .fastq.gz)
   ls -la data/raw/
   # Or specify custom directory:
   ./scripts/quality_control.sh --input-dir /path/to/your/files
   ```

4. **Low mapping rate (<80%)**
   ```bash
   ./scripts/quality_control.sh  # Check read quality first
   # Verify reference genome matches your sample species
   ```

**For detailed solutions**, see [TROUBLESHOOTING.md](TROUBLESHOOTING.md)

## 📊 Expected Results

From 30x WGS data, expect:
- **Total variants**: 4-5 million
- **SNVs**: ~4 million
- **Indels**: ~0.8 million
- **Novel variants**: 5-10%
- **High-impact variants**: 300-500

## 🤝 Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Areas for Contribution
- Additional variant callers (GATK, Strelka2)
- Structural variant detection
- Copy number variation analysis
- Improved visualization tools
- Additional population databases

## 📖 Citations

If you use this pipeline, please cite:

```bibtex
@software{wgs_analysis_pipeline,
  title = {Whole Genome Sequencing Analysis Pipeline},
  year = {2025},
  url = {https://github.com/yourusername/wgs-analysis-pipeline}
}
```

### Key Tools Used
- BWA: [Li and Durbin, 2009](https://doi.org/10.1093/bioinformatics/btp324)
- SAMtools: [Danecek et al., 2021](https://doi.org/10.1093/gigascience/giab008)
- BCFtools: [Danecek et al., 2021](https://doi.org/10.1093/gigascience/giab008)
- DeepVariant: [Poplin et al., 2018](https://doi.org/10.1038/nbt.4235)
- VEP: [McLaren et al., 2016](https://doi.org/10.1186/s13059-016-0974-4)

## 📝 License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

## ⚠️ Disclaimer

This pipeline is for research and educational purposes. For clinical applications, please use validated clinical pipelines and consult with qualified geneticists and healthcare providers.

## 🙏 Acknowledgments

- Genome Reference Consortium for the GRCh38 human reference
- gnomAD team for population frequency data
- Ensembl VEP team for annotation tools
- The open-source bioinformatics community

## 📧 Contact

For questions or issues, please open a GitHub issue or see the documentation.

---

**Note**: This repository contains only pipeline code and documentation. No personal genomic data or identifying information is included.