#!/bin/bash

# VEP Installation and Configuration Script for CHM13 T2T Annotation
# Optimized for n2-standard-32 VM (32 cores, 128GB RAM)
# Author: WGS Analysis Pipeline
# Date: July 22, 2025

set -e  # Exit on any error
set -u  # Exit on undefined variables

# Configuration
GENOMICS_DIR="/mnt/genomics"
VEP_DIR="${GENOMICS_DIR}/vep"
CACHE_DIR="${VEP_DIR}/cache"
DB_DIR="${VEP_DIR}/databases"
INPUT_VCF="${GENOMICS_DIR}/results/bcftools_variants_parallel.vcf.gz"
OUTPUT_DIR="${GENOMICS_DIR}/results/annotation"
TEMP_DIR="${VEP_DIR}/temp"
LOG_FILE="${GENOMICS_DIR}/vep_setup.log"

# VM Resources
CORES=$(nproc)
MEMORY_GB=$(($(free -g | awk '/^Mem:/{print $2}') - 8))  # Leave 8GB for system
FORK_PROCESSES=$((CORES / 4))  # VEP forks, optimize based on I/O

echo "=== VEP Setup for CHM13 T2T Annotation ===" | tee -a $LOG_FILE
echo "VM Resources: ${CORES} cores, ${MEMORY_GB}GB available memory" | tee -a $LOG_FILE
echo "VEP Configuration: ${FORK_PROCESSES} parallel processes" | tee -a $LOG_FILE
echo "Started: $(date)" | tee -a $LOG_FILE

# Create directory structure
echo "Creating directory structure..." | tee -a $LOG_FILE
mkdir -p $VEP_DIR $CACHE_DIR $DB_DIR $OUTPUT_DIR $TEMP_DIR

# Update system and install dependencies
echo "Installing dependencies..." | tee -a $LOG_FILE
sudo apt-get update -y
sudo apt-get install -y \
    build-essential \
    git \
    perl \
    cpanminus \
    libdbi-perl \
    libdbd-mysql-perl \
    libbz2-dev \
    zlib1g-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libperl-dev \
    libmodule-build-perl \
    wget \
    curl \
    unzip \
    tabix

# Install VEP from Ensembl
echo "Installing VEP..." | tee -a $LOG_FILE
cd $VEP_DIR
if [ ! -d "ensembl-vep" ]; then
    git clone https://github.com/Ensembl/ensembl-vep.git
    cd ensembl-vep
    
    # Install VEP with all dependencies
    perl INSTALL.pl \
        --AUTO a \
        --SPECIES homo_sapiens \
        --ASSEMBLY CHM13 \
        --NO_UPDATE \
        --CACHEDIR $CACHE_DIR \
        --CONVERT
        
    echo "VEP installation completed" | tee -a $LOG_FILE
else
    echo "VEP already installed, updating..." | tee -a $LOG_FILE
    cd ensembl-vep
    git pull
fi

# Set up environment
export PATH="${VEP_DIR}/ensembl-vep:$PATH"
export PERL5LIB="${VEP_DIR}/ensembl-vep:$PERL5LIB"

# Create VEP configuration file for optimal performance
cat > ${VEP_DIR}/vep.ini << EOF
# VEP Configuration for CHM13 T2T Annotation
# Optimized for 32-core VM with 128GB RAM

# Performance settings
fork                = ${FORK_PROCESSES}
buffer_size         = 50000
check_existing      = 1

# Assembly and cache
assembly            = CHM13
cache               = 1
dir_cache           = ${CACHE_DIR}
offline             = 1

# Output format
format              = vcf
compress_output     = bgzip
vcf                 = 1

# Annotation options
everything          = 1
gene_phenotype      = 1
regulatory          = 1
protein             = 1
symbol              = 1
ccds                = 1
uniprot             = 1
biotype             = 1
tsl                 = 1
canonical           = 1
variant_class       = 1

# Population frequencies
af_gnomad           = 1
af_1kg              = 1

# Prediction scores
sift                = both
polyphen            = both
domains             = 1
numbers             = 1

# Clinical significance
clin_sig_allele     = 1

# Additional annotations
pubmed              = 1
var_synonyms        = 1

# Temporary directory
tmpdir              = ${TEMP_DIR}

# Statistics
stats_file          = ${OUTPUT_DIR}/vep_summary.html
EOF

echo "VEP configuration file created" | tee -a $LOG_FILE

# Create database download script
cat > ${VEP_DIR}/download_databases.sh << 'EOF'
#!/bin/bash

echo "=== Downloading Annotation Databases ==="
DB_DIR="/mnt/genomics/vep/databases"
CACHE_DIR="/mnt/genomics/vep/cache"

# Function to download with retry
download_with_retry() {
    local url=$1
    local output=$2
    local max_attempts=3
    local attempt=1
    
    while [ $attempt -le $max_attempts ]; do
        echo "Downloading $url (attempt $attempt/$max_attempts)"
        if wget -c -O "$output" "$url"; then
            echo "Download successful"
            return 0
        else
            echo "Download failed, retrying..."
            ((attempt++))
            sleep 10
        fi
    done
    
    echo "Failed to download $url after $max_attempts attempts"
    return 1
}

# Download CHM13 VEP cache
echo "Downloading VEP cache for CHM13..."
cd $CACHE_DIR
VEP_CACHE_URL="ftp://ftp.ensembl.org/pub/current_variation/indexed_vep_cache/homo_sapiens_vep_110_CHM13.tar.gz"
if [ ! -f "homo_sapiens_vep_110_CHM13.tar.gz" ]; then
    download_with_retry $VEP_CACHE_URL "homo_sapiens_vep_110_CHM13.tar.gz"
    echo "Extracting VEP cache..."
    tar -xzf homo_sapiens_vep_110_CHM13.tar.gz
fi

# Download ClinVar
echo "Downloading ClinVar..."
cd $DB_DIR
CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
if [ ! -f "clinvar.vcf.gz" ]; then
    download_with_retry $CLINVAR_URL "clinvar.vcf.gz"
    tabix -p vcf clinvar.vcf.gz
fi

# Download gnomAD v4 (subset for CHM13 if available)
echo "Downloading gnomAD v4..."
GNOMAD_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.vcf.bgz"
if [ ! -f "gnomad.genomes.v4.0.sites.vcf.bgz" ]; then
    download_with_retry $GNOMAD_URL "gnomad.genomes.v4.0.sites.vcf.bgz"
    tabix -p vcf gnomad.genomes.v4.0.sites.vcf.bgz
fi

# Download dbNSFP
echo "Downloading dbNSFP..."
DBNSFP_URL="https://dbnsfp.softgenetics.com/dbNSFP4.5a.zip"
if [ ! -f "dbNSFP4.5a.zip" ]; then
    download_with_retry $DBNSFP_URL "dbNSFP4.5a.zip"
    echo "Extracting dbNSFP..."
    unzip -o dbNSFP4.5a.zip
    # Prepare for VEP
    zcat dbNSFP*_variant.chr* | head -1 > dbNSFP_header
    zcat dbNSFP*_variant.chr* | grep -v "^#" | sort -k1,1 -k2,2n | bgzip > dbNSFP.txt.gz
    tabix -s 1 -b 2 -e 2 -c '#' dbNSFP.txt.gz
fi

echo "All databases downloaded successfully"
EOF

chmod +x ${VEP_DIR}/download_databases.sh

# Create optimized annotation script
cat > ${VEP_DIR}/run_annotation.sh << EOF
#!/bin/bash

# Optimized VEP Annotation Script for 32-core VM
# Input: bcftools_variants_parallel.vcf.gz (4.75M variants)

set -e

GENOMICS_DIR="/mnt/genomics"
VEP_DIR="\${GENOMICS_DIR}/vep"
INPUT_VCF="\${GENOMICS_DIR}/results/bcftools_variants_parallel.vcf.gz"
OUTPUT_DIR="\${GENOMICS_DIR}/results/annotation"
CONFIG_FILE="\${VEP_DIR}/vep.ini"
LOG_FILE="\${OUTPUT_DIR}/annotation.log"

# Verify input exists
if [ ! -f "\$INPUT_VCF" ]; then
    echo "Error: Input VCF not found: \$INPUT_VCF"
    exit 1
fi

echo "=== VEP Annotation Started ===" | tee -a \$LOG_FILE
echo "Input: \$INPUT_VCF" | tee -a \$LOG_FILE
echo "Cores: $(nproc)" | tee -a \$LOG_FILE
echo "Memory: \$(free -h | grep ^Mem: | awk '{print \$2}')" | tee -a \$LOG_FILE
echo "Started: \$(date)" | tee -a \$LOG_FILE

# Set environment
export PATH="\${VEP_DIR}/ensembl-vep:\$PATH"

# Run VEP annotation with optimal settings
\${VEP_DIR}/ensembl-vep/vep \\
    --config \$CONFIG_FILE \\
    --input_file \$INPUT_VCF \\
    --output_file \${OUTPUT_DIR}/annotated_variants.vcf.gz \\
    --stats_file \${OUTPUT_DIR}/annotation_stats.html \\
    --warning_file \${OUTPUT_DIR}/annotation_warnings.txt \\
    --custom \${VEP_DIR}/databases/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \\
    --custom \${VEP_DIR}/databases/gnomad.genomes.v4.0.sites.vcf.bgz,gnomAD,vcf,exact,0,AF,AF_popmax \\
    --plugin dbNSFP,\${VEP_DIR}/databases/dbNSFP.txt.gz,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,CADD_raw,CADD_phred,REVEL_score \\
    --verbose | tee -a \$LOG_FILE

echo "=== VEP Annotation Completed ===" | tee -a \$LOG_FILE
echo "Finished: \$(date)" | tee -a \$LOG_FILE

# Generate summary statistics
echo "Generating annotation summary..." | tee -a \$LOG_FILE
bcftools stats \${OUTPUT_DIR}/annotated_variants.vcf.gz > \${OUTPUT_DIR}/annotated_variants_stats.txt

echo "Annotation pipeline completed successfully" | tee -a \$LOG_FILE
EOF

chmod +x ${VEP_DIR}/run_annotation.sh

echo "=== VEP Setup Script Creation Completed ===" | tee -a $LOG_FILE
echo "Created files:" | tee -a $LOG_FILE
echo "  - VEP configuration: ${VEP_DIR}/vep.ini" | tee -a $LOG_FILE
echo "  - Database downloader: ${VEP_DIR}/download_databases.sh" | tee -a $LOG_FILE
echo "  - Annotation runner: ${VEP_DIR}/run_annotation.sh" | tee -a $LOG_FILE
echo "" | tee -a $LOG_FILE
echo "Next steps:" | tee -a $LOG_FILE
echo "1. Run this setup script on VM: bash vm_vep_setup.sh" | tee -a $LOG_FILE
echo "2. Download databases: bash ${VEP_DIR}/download_databases.sh" | tee -a $LOG_FILE
echo "3. Run annotation: bash ${VEP_DIR}/run_annotation.sh" | tee -a $LOG_FILE
echo "Completed: $(date)" | tee -a $LOG_FILE