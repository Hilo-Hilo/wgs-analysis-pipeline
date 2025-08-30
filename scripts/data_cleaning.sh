#!/bin/bash

# Data Cleaning Pipeline with fastp
# This script performs minimal cleaning on high-quality WGS data
# Usage: ./scripts/data_cleaning.sh

set -e  # Exit on any error

# Configuration
CONDA_ENV="wgs_analysis"
RAW_R1="data/raw/SAMPLE001_L01_UDB-406_1.fq.gz"
RAW_R2="data/raw/SAMPLE001_L01_UDB-406_2.fq.gz"
CLEANED_R1="data/processed/SAMPLE001_clean_R1.fq.gz"
CLEANED_R2="data/processed/SAMPLE001_clean_R2.fq.gz"
LOG_DIR="logs"
THREADS=8

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Logging function
log() {
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1" | tee -a "$LOG_DIR/data_cleaning.log"
}

error() {
    echo -e "${RED}[$(date '+%Y-%m-%d %H:%M:%S')] ERROR:${NC} $1" | tee -a "$LOG_DIR/data_cleaning.log"
}

warning() {
    echo -e "${YELLOW}[$(date '+%Y-%m-%d %H:%M:%S')] WARNING:${NC} $1" | tee -a "$LOG_DIR/data_cleaning.log"
}

# Check prerequisites
check_prerequisites() {
    log "Checking prerequisites for data cleaning..."
    
    # Check conda environment
    if [[ "$CONDA_DEFAULT_ENV" != "$CONDA_ENV" ]]; then
        error "Please activate conda environment: conda activate $CONDA_ENV"
        exit 1
    fi
    
    # Check raw reads
    if [[ ! -f "$RAW_R1" ]]; then
        error "Raw forward reads not found: $RAW_R1"
        exit 1
    fi
    
    if [[ ! -f "$RAW_R2" ]]; then
        error "Raw reverse reads not found: $RAW_R2"
        exit 1
    fi
    
    # Check tools
    if ! command -v fastp &> /dev/null; then
        error "fastp not found. Please install fastp."
        exit 1
    fi
    
    log "✓ Prerequisites check passed"
}

# Setup directories
setup_directories() {
    log "Setting up directories..."
    
    mkdir -p data/processed
    mkdir -p "$LOG_DIR"
    
    log "✓ Directories created"
}

# Get file information
get_file_info() {
    log "Getting file information..."
    
    if [[ -f "$RAW_R1" ]]; then
        r1_size=$(du -h "$RAW_R1" | cut -f1)
        log "Raw R1 reads: $r1_size"
    fi
    
    if [[ -f "$RAW_R2" ]]; then
        r2_size=$(du -h "$RAW_R2" | cut -f1)
        log "Raw R2 reads: $r2_size"
    fi
}

# Run fastp cleaning
run_fastp_cleaning() {
    log "Starting fastp data cleaning..."
    
    local fastp_log="$LOG_DIR/fastp_cleaning.log"
    local fastp_report="results/fastp_report.html"
    local fastp_json="results/fastp_report.json"
    
    # Check if output already exists
    if [[ -f "$CLEANED_R1" && -f "$CLEANED_R2" ]]; then
        warning "Cleaned reads already exist"
        log "Skipping fastp cleaning step"
        return 0
    fi
    
    # Run fastp with minimal cleaning (data quality is excellent)
    log "Running fastp with minimal cleaning parameters..."
    
    if fastp \
        -i "$RAW_R1" \
        -I "$RAW_R2" \
        -o "$CLEANED_R1" \
        -O "$CLEANED_R2" \
        --html "$fastp_report" \
        --json "$fastp_json" \
        --thread "$THREADS" \
        --detect_adapter_for_pe \
        --cut_tail \
        --cut_tail_window_size 4 \
        --cut_tail_mean_quality 20 \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 10 \
        --length_required 50 \
        --correction \
        --overrepresentation_analysis \
        2> "$fastp_log"; then
        
        log "✓ fastp cleaning completed successfully"
        
        # Get cleaned file sizes
        cleaned_r1_size=$(du -h "$CLEANED_R1" | cut -f1)
        cleaned_r2_size=$(du -h "$CLEANED_R2" | cut -f1)
        log "Cleaned R1 reads: $cleaned_r1_size"
        log "Cleaned R2 reads: $cleaned_r2_size"
        
        return 0
    else
        error "fastp cleaning failed"
        return 1
    fi
}

# Generate summary report
generate_summary() {
    log "Generating cleaning summary report..."
    
    local summary_file="data/processed/cleaning_summary.txt"
    
    cat > "$summary_file" << EOF
# Data Cleaning Summary - SAMPLE001
# Cleaning Date: $(date)
# Tool: fastp with minimal cleaning parameters

## Input Files:
- Forward reads: $RAW_R1
- Reverse reads: $RAW_R2

## Cleaning Parameters:
- Tool: fastp
- Threads: $THREADS
- Quality threshold: Q20
- Minimum length: 50bp
- Adapter detection: Automatic paired-end
- Tail cutting: Window size 4, mean quality 20
- Correction: Enabled
- Overrepresentation analysis: Enabled

## Output Files:
- Cleaned forward reads: $CLEANED_R1
- Cleaned reverse reads: $CLEANED_R2
- HTML report: results/fastp_report.html
- JSON report: results/fastp_report.json

## File Sizes:
EOF
    
    # Add file sizes
    if [[ -f "$RAW_R1" ]]; then
        raw_r1_size=$(du -h "$RAW_R1" | cut -f1)
        echo "- Raw R1: $raw_r1_size" >> "$summary_file"
    fi
    
    if [[ -f "$RAW_R2" ]]; then
        raw_r2_size=$(du -h "$RAW_R2" | cut -f1)
        echo "- Raw R2: $raw_r2_size" >> "$summary_file"
    fi
    
    if [[ -f "$CLEANED_R1" ]]; then
        cleaned_r1_size=$(du -h "$CLEANED_R1" | cut -f1)
        echo "- Cleaned R1: $cleaned_r1_size" >> "$summary_file"
    fi
    
    if [[ -f "$CLEANED_R2" ]]; then
        cleaned_r2_size=$(du -h "$CLEANED_R2" | cut -f1)
        echo "- Cleaned R2: $cleaned_r2_size" >> "$summary_file"
    fi
    
    cat >> "$summary_file" << EOF

## Quality Assessment:
- Original data quality: EXCELLENT (0 poor quality sequences from FastQC)
- Cleaning approach: Minimal cleaning to preserve high-quality data
- Expected impact: Minimal read loss, adapter removal, quality trimming

## Next Steps:
1. Review fastp HTML report
2. Proceed with CHM13 mapping pipeline
3. Run: ./scripts/chm13_mapping.sh

## Cleaning completed: $(date)
EOF
    
    log "✓ Summary report generated: $summary_file"
}

# Main execution
main() {
    log "Starting Data Cleaning Pipeline for WGS Data"
    log "==========================================="
    
    # Run all steps
    check_prerequisites
    setup_directories
    get_file_info
    run_fastp_cleaning || exit 1
    generate_summary
    
    log "Data Cleaning Pipeline Completed Successfully!"
    log "============================================="
    log "Review the fastp report: results/fastp_report.html"
    log "Summary file: data/processed/cleaning_summary.txt"
    log "Next: Run CHM13 mapping pipeline with ./scripts/chm13_mapping.sh"
}

# Handle interrupts gracefully
trap 'error "Script interrupted by user"; exit 1' INT TERM

# Run main function
main "$@"