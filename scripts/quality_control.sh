#!/bin/bash

# WGS Quality Control Analysis Script
# This script performs comprehensive quality control analysis on raw FASTQ files
# Usage: ./scripts/quality_control.sh

set -e  # Exit on any error

# Configuration
CONDA_ENV="wgs_analysis"
RAW_DATA_DIR="data/raw"
OUTPUT_DIR="results/fastqc_raw"
LOG_DIR="logs"
THREADS=8

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Logging function
log() {
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1" | tee -a "$LOG_DIR/quality_control.log"
}

error() {
    echo -e "${RED}[$(date '+%Y-%m-%d %H:%M:%S')] ERROR:${NC} $1" | tee -a "$LOG_DIR/quality_control.log"
}

warning() {
    echo -e "${YELLOW}[$(date '+%Y-%m-%d %H:%M:%S')] WARNING:${NC} $1" | tee -a "$LOG_DIR/quality_control.log"
}

# Check prerequisites
check_prerequisites() {
    log "Checking prerequisites..."
    
    # Check conda environment
    if [[ "$CONDA_DEFAULT_ENV" != "$CONDA_ENV" ]]; then
        error "Please activate conda environment: conda activate $CONDA_ENV"
        exit 1
    fi
    
    # Check FastQC installation
    if ! command -v fastqc &> /dev/null; then
        error "FastQC not found. Please install FastQC."
        exit 1
    fi
    
    # Check input files
    if [[ ! -f "$RAW_DATA_DIR/SAMPLE001_L01_UDB-406_1.fq.gz" ]]; then
        error "Forward reads file not found: $RAW_DATA_DIR/SAMPLE001_L01_UDB-406_1.fq.gz"
        exit 1
    fi
    
    if [[ ! -f "$RAW_DATA_DIR/SAMPLE001_L01_UDB-406_2.fq.gz" ]]; then
        error "Reverse reads file not found: $RAW_DATA_DIR/SAMPLE001_L01_UDB-406_2.fq.gz"
        exit 1
    fi
    
    log "✓ Prerequisites check passed"
}

# Setup directories
setup_directories() {
    log "Setting up directories..."
    
    mkdir -p "$OUTPUT_DIR"
    mkdir -p "$LOG_DIR"
    
    log "✓ Directories created"
}

# Get file information
get_file_info() {
    log "Getting file information..."
    
    for file in "$RAW_DATA_DIR"/*.fq.gz; do
        if [[ -f "$file" ]]; then
            filename=$(basename "$file")
            filesize=$(du -h "$file" | cut -f1)
            log "File: $filename, Size: $filesize"
        fi
    done
}

# Run FastQC analysis
run_fastqc() {
    local input_file=$1
    local filename=$(basename "$input_file" .fq.gz)
    
    log "Starting FastQC analysis for $filename..."
    
    # Check if output already exists
    if [[ -f "$OUTPUT_DIR/${filename}_fastqc.html" ]]; then
        warning "FastQC output already exists for $filename. Skipping..."
        return 0
    fi
    
    # Run FastQC with error handling
    if fastqc \
        --threads "$THREADS" \
        --quiet \
        --extract \
        --outdir "$OUTPUT_DIR" \
        "$input_file" 2>&1 | tee -a "$LOG_DIR/fastqc_${filename}.log"; then
        log "✓ FastQC analysis completed for $filename"
        return 0
    else
        error "FastQC analysis failed for $filename"
        return 1
    fi
}

# Analyze FastQC results
analyze_results() {
    local fastqc_data_file=$1
    local filename=$(basename "$fastqc_data_file" _fastqc)
    
    log "Analyzing results for $filename..."
    
    if [[ ! -f "$fastqc_data_file/fastqc_data.txt" ]]; then
        error "FastQC data file not found: $fastqc_data_file/fastqc_data.txt"
        return 1
    fi
    
    # Extract key metrics
    local total_sequences=$(grep "Total Sequences" "$fastqc_data_file/fastqc_data.txt" | cut -f2)
    local poor_quality=$(grep "Sequences flagged as poor quality" "$fastqc_data_file/fastqc_data.txt" | cut -f2)
    local sequence_length=$(grep "Sequence length" "$fastqc_data_file/fastqc_data.txt" | cut -f2)
    local gc_content=$(grep "%GC" "$fastqc_data_file/fastqc_data.txt" | cut -f2)
    
    log "Results for $filename:"
    log "  Total sequences: $total_sequences"
    log "  Poor quality sequences: $poor_quality"
    log "  Sequence length: $sequence_length"
    log "  GC content: $gc_content%"
    
    # Check for warnings and failures
    local warnings=$(grep -c "WARN" "$fastqc_data_file/fastqc_data.txt" || echo "0")
    local failures=$(grep -c "FAIL" "$fastqc_data_file/fastqc_data.txt" || echo "0")
    
    if [[ $failures -gt 0 ]]; then
        warning "  $failures quality checks FAILED"
    fi
    
    if [[ $warnings -gt 0 ]]; then
        warning "  $warnings quality checks generated WARNINGS"
    fi
    
    if [[ $failures -eq 0 && $warnings -eq 0 ]]; then
        log "  ✓ All quality checks PASSED"
    fi
}

# Generate summary report
generate_summary() {
    local summary_file="$OUTPUT_DIR/quality_control_summary.txt"
    
    log "Generating summary report..."
    
    cat > "$summary_file" << EOF
# Quality Control Summary - SAMPLE001
# Analysis Date: $(date)
# Analyst: Automated Script

## Input Files:
- Forward reads: SAMPLE001_L01_UDB-406_1.fq.gz
- Reverse reads: SAMPLE001_L01_UDB-406_2.fq.gz

## Analysis Parameters:
- FastQC version: $(fastqc --version | head -1)
- Threads used: $THREADS
- Output directory: $OUTPUT_DIR

## File Information:
EOF
    
    for file in "$RAW_DATA_DIR"/*.fq.gz; do
        if [[ -f "$file" ]]; then
            filename=$(basename "$file")
            filesize=$(du -h "$file" | cut -f1)
            echo "- $filename: $filesize" >> "$summary_file"
        fi
    done
    
    cat >> "$summary_file" << EOF

## Quality Metrics:
EOF
    
    # Add quality metrics for each file
    for fastqc_dir in "$OUTPUT_DIR"/*_fastqc; do
        if [[ -d "$fastqc_dir" ]]; then
            filename=$(basename "$fastqc_dir" _fastqc)
            echo "### $filename:" >> "$summary_file"
            
            if [[ -f "$fastqc_dir/fastqc_data.txt" ]]; then
                local total_sequences=$(grep "Total Sequences" "$fastqc_dir/fastqc_data.txt" | cut -f2)
                local poor_quality=$(grep "Sequences flagged as poor quality" "$fastqc_dir/fastqc_data.txt" | cut -f2)
                local sequence_length=$(grep "Sequence length" "$fastqc_dir/fastqc_data.txt" | cut -f2)
                local gc_content=$(grep "%GC" "$fastqc_dir/fastqc_data.txt" | cut -f2)
                
                echo "- Total sequences: $total_sequences" >> "$summary_file"
                echo "- Poor quality sequences: $poor_quality" >> "$summary_file"
                echo "- Sequence length: $sequence_length" >> "$summary_file"
                echo "- GC content: $gc_content%" >> "$summary_file"
                
                local warnings=$(grep -c "WARN" "$fastqc_dir/fastqc_data.txt" || echo "0")
                local failures=$(grep -c "FAIL" "$fastqc_dir/fastqc_data.txt" || echo "0")
                
                echo "- Warnings: $warnings" >> "$summary_file"
                echo "- Failures: $failures" >> "$summary_file"
            fi
            echo "" >> "$summary_file"
        fi
    done
    
    cat >> "$summary_file" << EOF
## Recommendations:
- Review HTML reports for detailed quality assessment
- Check for adapter contamination and quality score distribution
- Proceed with data cleaning if quality is acceptable
- Consider resequencing if major quality issues detected

## Next Steps:
1. Review FastQC HTML reports
2. Make quality-based decisions for data cleaning
3. Proceed with fastp trimming if quality is acceptable
4. Document any quality issues for downstream analysis

## Files Generated:
EOF
    
    ls -la "$OUTPUT_DIR"/*.html >> "$summary_file" 2>/dev/null || echo "No HTML files found" >> "$summary_file"
    
    log "✓ Summary report generated: $summary_file"
}

# Open results for review
open_results() {
    log "Opening results for review..."
    
    # Open HTML reports
    for html_file in "$OUTPUT_DIR"/*.html; do
        if [[ -f "$html_file" ]]; then
            log "Opening: $html_file"
            open "$html_file" 2>/dev/null || {
                warning "Could not open $html_file automatically"
                log "Please manually open: $html_file"
            }
        fi
    done
}

# Cleanup function
cleanup() {
    log "Cleaning up temporary files..."
    
    # Remove any temporary files if needed
    # Currently no cleanup needed
    
    log "✓ Cleanup completed"
}

# Main execution
main() {
    log "Starting Quality Control Analysis for WGS Data"
    log "============================================="
    
    # Run all steps
    check_prerequisites
    setup_directories
    get_file_info
    
    # Run FastQC on both files
    run_fastqc "$RAW_DATA_DIR/SAMPLE001_L01_UDB-406_1.fq.gz" || exit 1
    run_fastqc "$RAW_DATA_DIR/SAMPLE001_L01_UDB-406_2.fq.gz" || exit 1
    
    # Analyze results
    for fastqc_dir in "$OUTPUT_DIR"/*_fastqc; do
        if [[ -d "$fastqc_dir" ]]; then
            analyze_results "$fastqc_dir"
        fi
    done
    
    # Generate summary
    generate_summary
    
    # Open results
    open_results
    
    # Cleanup
    cleanup
    
    log "Quality Control Analysis Completed Successfully!"
    log "=============================================="
    log "Review the HTML reports and summary file for detailed results"
    log "Summary file: $OUTPUT_DIR/quality_control_summary.txt"
}

# Handle interrupts gracefully
trap 'error "Script interrupted by user"; exit 1' INT TERM

# Run main function
main "$@"