# Troubleshooting Guide

This guide addresses common issues you may encounter when running the WGS analysis pipeline. Issues are organized by category with clear solutions and prevention strategies.

## 🔍 Quick Diagnostics

Before diving into specific issues, run these diagnostic commands:

```bash
# Check system requirements
./scripts/check_requirements.sh

# Verify conda environment
conda activate wgs_analysis
conda list | grep -E "(fastqc|fastp|bwa|samtools|bcftools)"

# Check disk space
df -h

# Check memory usage
free -h  # Linux
vm_stat  # macOS

# Check recent logs
ls -la logs/
tail logs/*.log
```

## 📊 Issue Categories

- [Environment & Installation](#environment--installation)
- [File & Directory Issues](#file--directory-issues)
- [Resource Problems](#resource-problems)
- [Quality Control Issues](#quality-control-issues)
- [Alignment Problems](#alignment-problems)
- [Variant Calling Issues](#variant-calling-issues)
- [Performance & Speed](#performance--speed)
- [Data Quality Concerns](#data-quality-concerns)

---

## 🔧 Environment & Installation

### Error: "Command not found" (fastqc, bwa, samtools, etc.)

**Symptoms:**
```bash
./scripts/quality_control.sh
bash: fastqc: command not found
```

**Cause:** Tools not installed or conda environment not activated

**Solution:**
```bash
# 1. Activate conda environment
conda activate wgs_analysis

# 2. If environment doesn't exist, create it
conda create -n wgs_analysis -c bioconda -c conda-forge \
    python=3.9 fastqc fastp bwa samtools bcftools

# 3. Verify installation
which fastqc bwa samtools
```

**Prevention:** Always activate environment before running scripts

### Error: "Please activate conda environment"

**Symptoms:**
```bash
ERROR: Please activate conda environment: conda activate wgs_analysis
```

**Solution:**
```bash
# Check current environment
echo $CONDA_DEFAULT_ENV

# Activate correct environment
conda activate wgs_analysis

# Make activation permanent (add to ~/.bashrc or ~/.zshrc)
echo "conda activate wgs_analysis" >> ~/.bashrc
```

### Error: Conda environment activation fails

**Symptoms:**
```bash
conda activate wgs_analysis
CommandNotFoundError: Your shell has not been configured to use 'conda activate'
```

**Solution:**
```bash
# Initialize conda for your shell
conda init bash  # or zsh, fish, etc.

# Restart terminal or run:
source ~/.bashrc

# Alternative: use conda directly
/path/to/conda/envs/wgs_analysis/bin/fastqc --version
```

---

## 📁 File & Directory Issues

### Error: "No FASTQ files found"

**Symptoms:**
```bash
ERROR: No FASTQ files found in: data/raw
Expected files with extensions: .fq.gz or .fastq.gz
```

**Solutions:**

**Option 1: Fix file extensions**
```bash
# Check what files you have
ls -la data/raw/

# Rename files if needed
cd data/raw/
rename 's/\.fastq$/.fastq.gz/' *.fastq  # If uncompressed
rename 's/\.fq$/.fq.gz/' *.fq          # Fix extension
gzip *.fastq  # Compress if needed
```

**Option 2: Specify correct directory**
```bash
./scripts/quality_control.sh --input-dir /path/to/your/fastq/files
```

**Option 3: Generate synthetic sample data for testing**
```bash
mkdir -p data/raw
python3 tests/generate_sample_data.py --output-dir data/raw --sample-name troubleshoot --num-reads 5000
./scripts/quality_control.sh --input-dir data/raw
```

### Error: "Permission denied"

**Symptoms:**
```bash
mkdir: cannot create directory 'results': Permission denied
```

**Solutions:**
```bash
# Check permissions
ls -la

# Fix permissions
chmod 755 .
sudo chown $USER:$USER . -R  # If needed

# Or run in a directory you own
cd $HOME
mkdir my_wgs_analysis
cd my_wgs_analysis
# Copy pipeline files here
```

### Error: "No space left on device"

**Symptoms:**
```bash
ERROR: No space left on device
```

**Immediate Solution:**
```bash
# Check disk usage
df -h
du -sh * | sort -hr  # Find large directories

# Clean up space
rm -rf logs/*.log        # Old logs
rm -rf temp/            # Temporary files
rm -rf results/*/temp/  # Intermediate files

# Compress old results
gzip results/*.bam results/*.sam
```

**Long-term Solution:**
- Move to partition with more space
- Use external storage
- Enable automatic cleanup in config
- Use cloud storage for large files

---

## 💾 Resource Problems

### Error: "Out of memory" / Process killed

**Symptoms:**
```bash
BWA alignment failed
Killed
# or
java.lang.OutOfMemoryError
```

**Immediate Solution:**
```bash
# Check available memory
free -h

# Reduce thread count (uses less memory)
./scripts/quality_control.sh --threads 4
./scripts/alignment.sh --threads 4

# Close other applications
# Kill memory-intensive processes
```

**Permanent Solution:**
```bash
# Create memory-optimized config
./scripts/load_config.sh create low_memory

# Edit config/low_memory.conf:
THREADS=4
MAX_MEMORY_BWA=8
MAX_MEMORY_VEP=6
KEEP_INTERMEDIATE_FILES=false
AUTO_CLEANUP=true

# Use the config
source config/low_memory.conf
./scripts/quality_control.sh
```

**Hardware Recommendations:**
- **Minimum**: 16GB RAM
- **Recommended**: 32GB RAM
- **High-throughput**: 64GB+ RAM

### Error: BWA-MEM2 memory requirements

**Symptoms:**
```bash
BWA-MEM2 index failed
Segmentation fault (core dumped)
```

**Cause:** BWA-MEM2 requires 128GB RAM for GRCh38 indexing

**Solution:**
```bash
# Use original BWA instead (requires less memory)
bwa index reference.fasta  # Instead of bwa-mem2 index

# Or use pre-built indexes
# Download from: https://bwa-mem2.github.io/
```

### Error: Process runs too slowly

**Symptoms:**
- Analysis taking much longer than expected
- High CPU usage but low progress

**Solutions:**
```bash
# Check system load
top
htop  # If available

# Optimize thread usage
# Rule: threads = CPU cores - 1
nproc  # Check number of CPU cores
./scripts/quality_control.sh --threads 7  # If you have 8 cores

# Use faster storage (SSD vs HDD)
# Move working directory to SSD if available

# Check I/O usage
iotop  # Linux
iostat 1  # Monitor I/O

# Enable compression to reduce I/O
COMPRESSION_LEVEL=9  # In config file
```

---

## 🔬 Quality Control Issues

### Warning: Low quality scores (mean < 30)

**Symptoms:**
```bash
WARNING: Mean quality score: 25 (below recommended 30)
```

**Assessment:**
```bash
# Check detailed FastQC report
open results/quality_control/*.html

# Look for:
# - Per base sequence quality
# - Per sequence quality scores
# - Overrepresented sequences (adapters)
```

**Solutions:**

**Option 1: Clean reads more aggressively**
```bash
# Edit quality thresholds in config
FASTP_QUALITY_THRESHOLD=25  # Increase from 20
FASTP_CUT_TAIL_QUALITY=25   # More aggressive tail trimming
FASTP_MIN_LENGTH=75         # Longer minimum length

# Re-run cleaning
./scripts/data_cleaning.sh --force
```

**Option 2: Filter during alignment**
```bash
# Use higher mapping quality filter
MIN_MAPPING_QUALITY=30  # Increase from 20
```

**Option 3: Contact sequencing provider**
- If quality is consistently poor across samples
- May indicate sequencing run issues

### Error: High adapter contamination (>10%)

**Symptoms:**
```bash
WARNING: High adapter contamination detected: 15%
```

**Solution:**
```bash
# More aggressive adapter trimming
FASTP_ADAPTER_TRIMMING=true
FASTP_AUTO_DETECT_ADAPTERS=true

# Manual adapter sequences (if known)
FASTP_ADAPTER_SEQUENCE_R1="AGATCGGAAGAGC"
FASTP_ADAPTER_SEQUENCE_R2="AGATCGGAAGAGC"

# Re-run with cleaning
./scripts/data_cleaning.sh --force
```

### Warning: Unusual GC content

**Symptoms:**
```bash
WARNING: GC content deviation from expected (~42% for human)
Observed: 65%
```

**Possible Causes:**
1. **Contamination** (bacterial, viral)
2. **Adapter dimers**
3. **PCR bias**
4. **Non-human sample**

**Solutions:**
```bash
# Check for contamination
# Run alignment to see what maps

# Screen against common contaminants
# Use tools like FastQ Screen

# If confirmed contamination, contact lab
```

---

## 🎯 Alignment Problems

### Error: Low mapping rate (<80%)

**Symptoms:**
```bash
WARNING: Low mapping rate: 75% (expected >85%)
```

**Diagnostic Steps:**
```bash
# 1. Check reference genome
ls -la data/reference/
head -5 data/reference/GRCh38/GRCh38_latest_genomic.fna

# 2. Check read quality
./scripts/quality_control.sh

# 3. Check for contamination
samtools view results/alignment/*.bam | head -100
```

**Solutions:**

**Option 1: Wrong reference genome**
```bash
# Human samples should use GRCh38
# Check your sample species matches reference

# Download GRCh38 reference manually
mkdir -p data/reference/GRCh38
cd data/reference/GRCh38
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
mv GCF_000001405.40_GRCh38.p14_genomic.fna GRCh38_latest_genomic.fna
samtools faidx GRCh38_latest_genomic.fna
bwa index GRCh38_latest_genomic.fna
cd ../../..
```

**Option 2: Poor read quality**
```bash
# More aggressive read cleaning
FASTP_QUALITY_THRESHOLD=25
FASTP_MIN_LENGTH=100

./scripts/data_cleaning.sh --force
./scripts/alignment.sh --force
```

**Option 3: Contamination**
```bash
# Screen for contamination
# Map to multiple references
# Contact sequencing facility
```

### Error: BWA index files missing

**Symptoms:**
```bash
ERROR: BWA index not found for reference
Please run: bwa index reference.fasta
```

**Solution:**
```bash
# Index your reference genome
cd data/reference/GRCh38/
bwa index GRCh38_latest_genomic.fna

# Also create samtools index
samtools faidx GRCh38_latest_genomic.fna

# Verify indexes created
ls -la *.bwt *.pac *.ann *.amb *.sa *.fai
```

### Error: Reference genome not found

**Symptoms:**
```bash
ERROR: Reference file not found: data/reference/GRCh38/GRCh38_latest_genomic.fna
```

**Solutions:**
```bash
# Option 1: Download reference
mkdir -p data/reference/GRCh38
cd data/reference/GRCh38
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
mv GCF_000001405.40_GRCh38.p14_genomic.fna GRCh38_latest_genomic.fna
samtools faidx GRCh38_latest_genomic.fna
bwa index GRCh38_latest_genomic.fna
cd ../../..

# Option 2: Specify correct path
./scripts/alignment.sh --reference /path/to/your/reference.fa

# Option 3: Create symlink
ln -s /path/to/reference.fa data/reference/GRCh38/GRCh38_latest_genomic.fna
```

---

## 🧬 Variant Calling Issues

### Warning: Very few variants called

**Symptoms:**
```bash
WARNING: Only 1,000 variants called (expected 4-5 million for WGS)
```

**Possible Causes:**
1. **Low coverage data**
2. **Stringent filtering**
3. **Reference genome mismatch**
4. **Poor alignment quality**

**Diagnostic Steps:**
```bash
# Check coverage depth
samtools depth results/alignment/*.bam | awk '{sum+=$3} END {print "Average depth:", sum/NR}'

# Check alignment stats
cat results/alignment/*_stats.txt
```

**Solutions:**

**Option 1: Reduce stringent filters**
```bash
# Lower quality thresholds
VARIANT_MIN_DEPTH=5     # Reduce from 10
VARIANT_MIN_QUALITY=20  # Reduce from 30

# Re-run variant calling
./scripts/variant_calling.sh --force
```

**Option 2: Check data type**
```bash
# Confirm this is WGS (not exome or targeted)
# Exome typically yields 20,000-50,000 variants
# WGS typically yields 4-5 million variants
```

### Error: Too many variants called

**Symptoms:**
```bash
WARNING: 20 million variants called (expected 4-5 million for WGS)
```

**Possible Causes:**
1. **Contamination**
2. **Poor quality data**
3. **Wrong reference genome**
4. **Too lenient filtering**

**Solutions:**
```bash
# More stringent filtering
VARIANT_MIN_DEPTH=15
VARIANT_MIN_QUALITY=50
BCFTOOLS_MIN_BASE_QUALITY=30

# Check for contamination
# Review alignment statistics
```

---

## ⚡ Performance & Speed

### Problem: Analysis taking much longer than expected

**Expected Times (30x WGS, 8 cores):**
- Quality Control: 10-30 minutes
- Read Cleaning: 30-60 minutes
- Alignment: 2-6 hours
- Variant Calling: 1-4 hours
- Annotation: 30-90 minutes

**Solutions:**

**Optimize thread usage:**
```bash
# Check CPU cores
nproc

# Use optimal thread count (cores - 1)
./scripts/quality_control.sh --threads 15  # If 16 cores
```

**Optimize I/O:**
```bash
# Use SSD storage if available
# Move temp files to ramdisk (if enough RAM)
export TMPDIR=/dev/shm  # Linux ramdisk

# Reduce compression during processing
COMPRESSION_LEVEL=1  # Fast compression
```

**Monitor system resources:**
```bash
# Real-time monitoring
htop  # CPU and memory usage
iotop  # Disk I/O usage
iftop  # Network usage (if using network storage)
```

### Problem: Running out of temporary disk space

**Symptoms:**
```bash
ERROR: No space left on device
```

**Solutions:**
```bash
# Change temp directory location
export TMPDIR=/path/to/larger/partition

# Clean temp files more frequently
AUTO_CLEANUP=true

# Use less disk-intensive options
KEEP_INTERMEDIATE_FILES=false
```

---

## 📈 Data Quality Concerns

### Question: How do I know if my results are good?

**Quality Metrics to Check:**

**Sequencing Quality:**
- Mean quality score: >30
- Adapter contamination: <5%
- Duplication rate: <20%

**Alignment Quality:**
- Mapping rate: >85%
- Properly paired: >95%
- Mean insert size: 300-500bp

**Variant Quality:**
- Total variants: 4-5M (WGS), 20-50K (exome)
- Ti/Tv ratio: 2.0-2.1
- Het/Hom ratio: 1.5-2.0

**Check these with:**
```bash
# Quality summary
cat results/quality_control/quality_control_summary.txt

# Alignment summary
cat results/alignment/*_alignment_summary.txt

# Variant summary
bcftools stats results/variants/*_filtered.vcf.gz
```

### Question: Should I be concerned about these warnings?

**Common warnings and their significance:**

**🟢 Usually OK:**
- Slight GC bias (38-45% for human)
- Some adapter content (<5%)
- Mapping rate 80-85% (still usable)

**🟡 Investigate further:**
- High duplication (>20%)
- Low mapping rate (70-80%)
- Unusual insert sizes

**🔴 Serious concerns:**
- Very low quality (<20)
- High adapter content (>10%)
- Very low mapping rate (<70%)
- Extremely high/low GC content

---

## 🆘 Emergency Fixes

### System completely out of disk space

```bash
# Quick cleanup (be careful!)
rm -rf logs/*.log
rm -rf */temp/
rm -rf */*.tmp
find . -name "*.sam" -size +1G -delete  # Remove large SAM files if BAM exists

# Emergency disk space
sudo apt-get clean  # Linux
brew cleanup        # macOS
```

### Process consuming all memory

```bash
# Find memory hogs
ps aux --sort=-%mem | head

# Kill specific process (replace PID)
kill -9 PID

# Kill all BWA processes
killall bwa

# Prevent future issues
ulimit -v 16777216  # Limit to 16GB
```

### Complete analysis restart

```bash
# Clean slate restart
rm -rf results/ logs/ temp/
mkdir -p results logs temp

# Start over with fresh settings
./scripts/check_requirements.sh
./scripts/quality_control.sh --verbose
```

---

## 📞 Getting Help

### Before asking for help, collect this information:

```bash
# System information
uname -a
conda --version
conda list -n wgs_analysis

# Error details
tail -50 logs/*.log
./scripts/check_requirements.sh

# File information
ls -la data/
df -h
```

### Where to get help:

1. **Documentation**
   - [GETTING_STARTED.md](GETTING_STARTED.md)
   - [16GB_SYSTEM_GUIDE.md](16GB_SYSTEM_GUIDE.md)
   - [INPUT_OUTPUT_SPECIFICATION.md](INPUT_OUTPUT_SPECIFICATION.md)

2. **Online Communities**
   - [Bioinformatics Stack Exchange](https://bioinformatics.stackexchange.com/)
   - [Reddit r/bioinformatics](https://reddit.com/r/bioinformatics)
   - [SeqAnswers Forum](http://seqanswers.com/)

3. **Tool-Specific Help**
   - FastQC: `fastqc --help`
   - BWA: `bwa mem`
   - samtools: `samtools help`
   - bcftools: `bcftools --help`

4. **GitHub Issues**
   - File bug reports
   - Feature requests
   - Documentation improvements

### When reporting issues, include:

- **What you were trying to do**
- **What command you ran**
- **Complete error message**
- **Your system information**
- **Log files (relevant portions)**
- **Steps to reproduce**

---

## 🔄 Prevention Strategies

### Before starting any analysis:

```bash
# 1. Always check requirements
./scripts/check_requirements.sh

# 2. Start with small test data
python3 tests/generate_sample_data.py --output-dir data/raw --sample-name smoke --num-reads 2000

# 3. Use dry-run mode first
./scripts/quality_control.sh --dry-run

# 4. Monitor resources during runs
htop &  # Keep running in background
```

### Regular maintenance:

```bash
# Weekly cleanup
rm -rf logs/*.log
rm -rf temp/*

# Monthly updates
conda update -n wgs_analysis --all

# Backup important configs
cp -r config/ backup/
```

### Best practices:

1. **Start small** - Test with sample data first
2. **Check resources** - Always verify before starting
3. **Monitor progress** - Watch logs and system usage  
4. **Keep backups** - Save important results
5. **Document changes** - Note what parameters you modify
6. **Version control** - Track configuration changes

---

**Remember:** Most issues have simple solutions. Take time to read error messages carefully, and don't hesitate to ask for help when needed. The bioinformatics community is generally very helpful!