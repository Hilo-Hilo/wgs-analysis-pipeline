# Getting Started with WGS Analysis Pipeline

This guide will walk you through running your first whole genome sequencing (WGS) analysis using this pipeline, designed for users with minimal bioinformatics experience but comfort with command-line tools.

## 🎯 Quick Start (5 minutes)

If you want to jump right in and test the pipeline:

```bash
# 1. Check if your system is ready
./scripts/check_requirements.sh

# 2. Download sample data for testing
./scripts/download_sample_data.sh --type small

# 3. Test quality control (takes ~2 minutes)
conda activate wgs_analysis
./scripts/quality_control.sh --input-dir data/samples/fastq --dry-run
```

## 📋 Prerequisites

Before starting, ensure you have:

- **Operating System**: Linux or macOS
- **RAM**: Minimum 16GB (32GB+ recommended)
- **Disk Space**: Minimum 500GB free
- **Internet Connection**: For downloading references and databases
- **Command Line**: Basic familiarity with terminal/command prompt

## 🔧 Step 1: System Setup

### Install Conda (if not already installed)

```bash
# Download Miniconda (lightweight conda installer)
# For Linux:
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# For macOS:
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh

# Restart your terminal or run:
source ~/.bashrc  # Linux
source ~/.zshrc   # macOS with zsh
```

### Create WGS Analysis Environment

```bash
# Create conda environment with all required tools
conda create -n wgs_analysis -c bioconda -c conda-forge \
    python=3.9 fastqc fastp bwa samtools bcftools multiqc

# Activate the environment
conda activate wgs_analysis

# Verify installation
conda list | grep -E "(fastqc|fastp|bwa|samtools|bcftools)"
```

## ✅ Step 2: Verify System Requirements

Run the comprehensive system check:

```bash
# Check all requirements
./scripts/check_requirements.sh

# If you see any failures, follow the provided recommendations
# The script will tell you exactly what to install or fix
```

**What this checks:**
- System RAM and disk space
- Required software tools
- Conda environment setup
- File permissions
- Resource estimates for analysis

## 📁 Step 3: Prepare Your Data

### Option A: Test with Sample Data (Recommended for First Run)

```bash
# Download small test dataset (100MB, ~5 minutes to analyze)
./scripts/download_sample_data.sh --type small

# This creates:
# data/samples/fastq/          - Test FASTQ files
# data/samples/expected_results/ - What results should look like
```

### Option B: Use Your Own Data

```bash
# Create data directory structure
mkdir -p data/raw data/processed results logs

# Copy your FASTQ files to data/raw/
# Files should be named: *_R1.fastq.gz and *_R2.fastq.gz
# Example: MySample_R1.fastq.gz, MySample_R2.fastq.gz

ls data/raw/  # Verify your files are there
```

## 🚀 Step 4: Run Your First Analysis

### Configure the Analysis

```bash
# Option 1: Use default settings (recommended for beginners)
# No configuration needed - defaults work for most cases

# Option 2: Create custom configuration (for advanced users)
./scripts/load_config.sh create my_first_analysis
# Then edit config/my_first_analysis.conf
```

### Quality Control Analysis

This is always the first step:

```bash
# Activate conda environment
conda activate wgs_analysis

# Run quality control (10-30 minutes depending on data size)
./scripts/quality_control.sh

# Or with custom input directory:
./scripts/quality_control.sh --input-dir data/samples/fastq --output-dir results/sample_qc

# View help for all options:
./scripts/quality_control.sh --help
```

**What happens:**
- Analyzes read quality, length distribution, GC content
- Identifies adapter contamination and quality issues
- Generates HTML reports you can view in a web browser
- Creates summary statistics

**Expected output:**
- `results/fastqc_raw/` - Quality control reports
- `results/fastqc_raw/quality_control_summary.txt` - Summary statistics

### Review Quality Control Results

```bash
# View the summary report
cat results/fastqc_raw/quality_control_summary.txt

# Open HTML reports in your browser (macOS):
open results/fastqc_raw/*.html

# Open HTML reports in your browser (Linux):
firefox results/fastqc_raw/*.html &
```

**What to look for:**
- ✅ **Good**: Mean quality >30, <5% adapter contamination
- ⚠️ **Warning**: Quality 20-30, some adapter contamination
- ❌ **Poor**: Quality <20, >10% adapters (consider re-sequencing)

## 🧬 Step 5: Complete Pipeline (Optional)

If quality control looks good, continue with the full pipeline:

### Download Reference Genome

```bash
# For beginners, we recommend GRCh38 (standard human reference)
mkdir -p data/reference/GRCh38

# Download GRCh38 reference (this may take 30-60 minutes)
cd data/reference/GRCh38
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
mv GCF_000001405.40_GRCh38.p14_genomic.fna GRCh38_latest_genomic.fna

# Index the reference (15-30 minutes)
samtools faidx GRCh38_latest_genomic.fna
bwa index GRCh38_latest_genomic.fna

cd ../../..  # Return to main directory
```

### Read Cleaning and Alignment

```bash
# Clean reads (removes adapters, low-quality bases)
./scripts/data_cleaning.sh

# Align to reference genome (2-6 hours for WGS)
./scripts/alignment.sh --reference data/reference/GRCh38/GRCh38_latest_genomic.fna
```

### Variant Calling and Annotation

```bash
# Call variants (1-4 hours)
./scripts/variant_calling.sh

# Annotate variants (30-90 minutes)
./scripts/vep_annotation.sh
```

## 📊 Understanding Your Results

### Directory Structure After Analysis

```
your_project/
├── data/
│   ├── raw/              # Your original FASTQ files
│   ├── processed/        # Cleaned FASTQ files
│   └── reference/        # Reference genome files
├── results/
│   ├── fastqc_raw/       # Quality control reports
│   ├── alignment/        # BAM files and mapping stats
│   └── variants/         # VCF files and annotations
└── logs/                 # Analysis log files
```

### Key Result Files

| File | Description | What to Check |
|------|-------------|---------------|
| `results/fastqc_raw/*.html` | Quality control reports | Read quality, contamination |
| `results/alignment/*_stats.txt` | Mapping statistics | Mapping rate (should be >85%) |
| `results/variants/*.vcf.gz` | Variant calls | Number of variants found |
| `results/variants/*_annotated.txt` | Annotated variants | High-impact variants |

## 🚨 Troubleshooting Common Issues

### Problem: "Command not found" errors
**Solution:**
```bash
# Make sure conda environment is activated
conda activate wgs_analysis

# Verify tools are installed
which fastqc bwa samtools
```

### Problem: "Out of memory" errors
**Solution:**
```bash
# Check available RAM
./scripts/check_requirements.sh

# Reduce thread count if needed
./scripts/quality_control.sh --threads 4
```

### Problem: "No space left on device"
**Solution:**
```bash
# Check disk space
df -h

# Clean up intermediate files
rm -rf temp/ results/*/temp/
```

### Problem: Low mapping rates (<80%)
**Possible causes:**
- Wrong reference genome
- Poor quality reads
- Contamination

**Solution:**
```bash
# Check read quality first
./scripts/quality_control.sh

# Try read cleaning
./scripts/data_cleaning.sh
```

## 📚 Next Steps

### For Research Users
1. **Explore variant annotations** - Look for variants in genes of interest
2. **Compare with databases** - Check against ClinVar, gnomAD
3. **Pathway analysis** - Analyze variants by biological pathways

### For Clinical Users
1. **Focus on high-impact variants** - Filter for likely pathogenic variants
2. **Check known disease genes** - Look at ACMG recommended genes
3. **Consult genetics professionals** - Interpretation requires expertise

### For Method Developers
1. **Benchmark against known variants** - Use reference samples
2. **Optimize parameters** - Tune variant calling thresholds
3. **Add new tools** - Integrate additional analysis methods

## 🔗 Additional Resources

### Documentation
- [Complete Analysis Workflow](documentation/complete_analysis_workflow.md) - Detailed pipeline steps
- [Cloud Infrastructure Guide](documentation/cloud_infrastructure_guide.md) - Running on cloud platforms
- [Environment Setup](documentation/environment_setup.md) - Advanced installation options
- [TROUBLESHOOTING.md](TROUBLESHOOTING.md) - Common problems and solutions

### External Resources
- [FastQC Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [BWA Manual](http://bio-bwa.sourceforge.net/bwa.shtml)
- [samtools Documentation](http://www.htslib.org/)
- [VEP Documentation](https://ensembl.org/info/docs/tools/vep/index.html)

### Getting Help
1. **Check logs** - Look in `logs/` directory for error messages
2. **Run diagnostics** - Use `./scripts/check_requirements.sh`
3. **Search documentation** - Check TROUBLESHOOTING.md
4. **Community forums** - Bioinformatics Stack Exchange, Reddit r/bioinformatics
5. **File issues** - Create GitHub issues for bugs or feature requests

## ⚠️ Important Notes

### Data Privacy
- This pipeline is for research/educational use
- **Never** upload personal genomic data to public repositories
- Follow your institution's data handling policies
- Consider encryption for sensitive data

### Resource Planning
- **Small dataset (1M reads)**: ~15 minutes, ~1GB storage
- **Medium dataset (10M reads)**: ~2 hours, ~10GB storage  
- **Full WGS (3B reads)**: ~24 hours, ~500GB storage

### Quality Expectations
- **Mapping rate**: 85-95% (human samples)
- **Variant count**: 4-5 million for whole genome
- **Ti/Tv ratio**: 2.0-2.1 (transition/transversion)
- **Het/Hom ratio**: 1.5-2.0 (heterozygous/homozygous)

---

**🎉 Congratulations!** You're now ready to start analyzing genomic data. Remember to start small with test data and gradually work up to full-scale analyses as you become more comfortable with the tools and concepts.

**Questions?** Check out [TROUBLESHOOTING.md](TROUBLESHOOTING.md) or the documentation in the `documentation/` folder.