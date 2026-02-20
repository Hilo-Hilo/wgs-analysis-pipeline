#!/bin/bash

# WGS Quality Control Analysis Script
# This script performs comprehensive quality control analysis on raw FASTQ files
# Usage: ./scripts/quality_control.sh [OPTIONS]

set -e  # Exit on any error

# Script version and info
SCRIPT_VERSION="2.0"
SCRIPT_NAME="quality_control.sh"

# Default Configuration (16GB RAM optimized)
CONDA_ENV="wgs_analysis"
RAW_DATA_DIR="data/raw"
OUTPUT_DIR="results/quality_control"
LOG_DIR="logs"
THREADS=4
DRY_RUN=false
VERBOSE=false
FORCE_OVERWRITE=false

# Help function
show_help() {
    cat << EOF
$SCRIPT_NAME v$SCRIPT_VERSION - WGS Quality Control Analysis

DESCRIPTION:
    Performs quality control analysis optimized for 16GB RAM systems.
    Uses FastQC with memory-efficient settings for local analysis.

USAGE:
    $0 [OPTIONS]

OPTIONS:
    -h, --help              Show this help message and exit
    -i, --input-dir DIR     Input directory containing FASTQ files (default: $RAW_DATA_DIR)
    -o, --output-dir DIR    Output directory for results (default: $OUTPUT_DIR)
    -t, --threads NUM       Number of threads to use (default: $THREADS)
    -e, --env ENV           Environment label for logs (default: $CONDA_ENV)
    -l, --log-dir DIR       Log directory (default: $LOG_DIR)
    --dry-run               Show what would be done without executing
    --force                 Overwrite existing output files
    --verbose               Enable verbose output
    --version               Show version information

EXAMPLES:
    # Basic usage (uses default directories)
    $0

    # Specify custom input directory
    $0 --input-dir /path/to/fastq/files

    # Run with more threads
    $0 --threads 16

    # Dry run to see what will be processed
    $0 --dry-run

    # Force overwrite of existing results
    $0 --force --verbose

REQUIREMENTS:
    - FastQC available in PATH
    - FASTQ files (.fq.gz or .fastq.gz) in input directory
    - Sufficient disk space for output files

OUTPUT:
    - FastQC HTML reports for each input file
    - Quality control summary report
    - Detailed analysis logs

AUTHOR:
    WGS Analysis Pipeline
    https://github.com/Hilo-Hilo/wgs-analysis-pipeline

EOF
}

# Version function
show_version() {
    echo "$SCRIPT_NAME version $SCRIPT_VERSION"
    exit 0
}

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Parse command line arguments
parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                show_help
                exit 0
                ;;
            --version)
                show_version
                ;;
            -i|--input-dir)
                RAW_DATA_DIR="$2"
                shift 2
                ;;
            -o|--output-dir)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            -t|--threads)
                THREADS="$2"
                shift 2
                ;;
            -e|--env)
                CONDA_ENV="$2"
                shift 2
                ;;
            -l|--log-dir)
                LOG_DIR="$2"
                shift 2
                ;;
            --dry-run)
                DRY_RUN=true
                shift
                ;;
            --force)
                FORCE_OVERWRITE=true
                shift
                ;;
            --verbose)
                VERBOSE=true
                shift
                ;;
            *)
                echo "Unknown option: $1"
                echo "Use --help for usage information"
                exit 1
                ;;
        esac
    done
}

# Logging functions
log() {
    local message="$1"
    local timestamp="[$(date '+%Y-%m-%d %H:%M:%S')]"
    echo -e "${GREEN}${timestamp}${NC} $message" | tee -a "$LOG_DIR/quality_control.log"
}

error() {
    local message="$1"
    local timestamp="[$(date '+%Y-%m-%d %H:%M:%S')]"
    echo -e "${RED}${timestamp} ERROR:${NC} $message" | tee -a "$LOG_DIR/quality_control.log" >&2
}

warning() {
    local message="$1"
    local timestamp="[$(date '+%Y-%m-%d %H:%M:%S')]"
    echo -e "${YELLOW}${timestamp} WARNING:${NC} $message" | tee -a "$LOG_DIR/quality_control.log"
}

info() {
    local message="$1"
    local timestamp="[$(date '+%Y-%m-%d %H:%M:%S')]"
    if [[ "$VERBOSE" == "true" ]]; then
        echo -e "${BLUE}${timestamp} INFO:${NC} $message" | tee -a "$LOG_DIR/quality_control.log"
    fi
}

# Progress indicator
show_progress() {
    local current=$1
    local total=$2
    local desc="$3"
    local percent=$((current * 100 / total))
    printf "\rProgress: [%-20s] %d%% %s" $(printf "#%.0s" $(seq 1 $((percent / 5)))) "$percent" "$desc"
    if [[ $current -eq $total ]]; then
        echo ""
    fi
}

# Validate input parameters
validate_parameters() {
    info "Validating parameters..."
    
    # Validate threads
    if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -le 0 ]]; then
        error "Invalid threads value: $THREADS. Must be a positive integer."
        exit 1
    fi
    
    # Convert relative paths to absolute paths
    RAW_DATA_DIR=$(realpath "$RAW_DATA_DIR" 2>/dev/null || echo "$RAW_DATA_DIR")
    OUTPUT_DIR=$(realpath "$OUTPUT_DIR" 2>/dev/null || echo "$OUTPUT_DIR")
    LOG_DIR=$(realpath "$LOG_DIR" 2>/dev/null || echo "$LOG_DIR")
    
    info "Using input directory: $RAW_DATA_DIR"
    info "Using output directory: $OUTPUT_DIR"
    info "Using threads: $THREADS"
}

# Check prerequisites
check_prerequisites() {
    log "Checking prerequisites..."
    
    # Check FastQC installation directly (portable across conda/docker/system installs)
    if ! command -v fastqc &> /dev/null; then
        error "FastQC not found in PATH."
        echo "Install FastQC (for example):"
        echo "  conda install -c bioconda fastqc"
        echo "Or ensure your environment/container PATH includes fastqc."
        exit 1
    fi
    
    # Check input directory
    if [[ ! -d "$RAW_DATA_DIR" ]]; then
        error "Input directory not found: $RAW_DATA_DIR"
        echo "Create the directory or specify a different path with --input-dir"
        exit 1
    fi
    
    # Find FASTQ files
    local fastq_files
    fastq_files=($(find "$RAW_DATA_DIR" \( -name "*.fq.gz" -o -name "*.fastq.gz" \) 2>/dev/null))
    
    if [[ ${#fastq_files[@]} -eq 0 ]]; then
        error "No FASTQ files found in: $RAW_DATA_DIR"
        echo "Expected files with extensions: .fq.gz or .fastq.gz"
        echo "Available files in directory:"
        ls -la "$RAW_DATA_DIR" 2>/dev/null || echo "  Directory is empty or unreadable"
        exit 1
    fi
    
    log "✓ Found ${#fastq_files[@]} FASTQ file(s)"
    for file in "${fastq_files[@]}"; do
        info "  $(basename "$file")"
    done
    
    # Check disk space
    local available_space
    available_space=$(df "$OUTPUT_DIR" 2>/dev/null | awk 'NR==2 {print $4}' || echo "0")
    if [[ "$available_space" -lt 1048576 ]]; then # Less than 1GB
        warning "Low disk space available (less than 1GB). Quality control may fail."
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
    local filename=$(basename "$input_file")
    local basename_clean="${filename%.fq.gz}"
    basename_clean="${basename_clean%.fastq.gz}"
    
    log "Starting FastQC analysis for $filename..."
    
    # Check if output already exists
    if [[ -f "$OUTPUT_DIR/${basename_clean}_fastqc.html" ]] && [[ "$FORCE_OVERWRITE" != "true" ]]; then
        warning "FastQC output already exists for $filename. Use --force to overwrite."
        return 0
    fi
    
    # Dry run check
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "Would run: fastqc --threads $THREADS --quiet --extract --outdir $OUTPUT_DIR $input_file"
        return 0
    fi
    
    # Validate input file
    if [[ ! -f "$input_file" ]]; then
        error "Input file not found: $input_file"
        return 1
    fi
    
    # Check file format (basic validation)
    if ! file "$input_file" | grep -q "gzip"; then
        warning "File may not be gzip compressed: $input_file"
    fi
    
    # Run FastQC with error handling
    info "Running FastQC with $THREADS threads on $filename"
    local fastqc_cmd="fastqc --threads $THREADS --quiet --extract --outdir $OUTPUT_DIR $input_file"
    
    if eval "$fastqc_cmd" 2>&1 | tee -a "$LOG_DIR/fastqc_${basename_clean}.log"; then
        log "✓ FastQC analysis completed for $filename"
        return 0
    else
        error "FastQC analysis failed for $filename"
        error "Check log file: $LOG_DIR/fastqc_${basename_clean}.log"
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
    # Parse command line arguments first
    parse_arguments "$@"
    
    log "Starting Quality Control Analysis for WGS Data"
    log "============================================="
    
    if [[ "$DRY_RUN" == "true" ]]; then
        log "DRY RUN MODE - No files will be processed"
    fi
    
    # Run all steps
    validate_parameters
    check_prerequisites
    setup_directories
    get_file_info
    
    # Find all FASTQ files in input directory
    local fastq_files
    fastq_files=($(find "$RAW_DATA_DIR" \( -name "*.fq.gz" -o -name "*.fastq.gz" \) 2>/dev/null))
    
    if [[ ${#fastq_files[@]} -eq 0 ]]; then
        error "No FASTQ files found to process"
        exit 1
    fi
    
    log "Processing ${#fastq_files[@]} FASTQ file(s)..."
    
    # Process each FASTQ file with progress indicator
    local failed_files=0
    for i in "${!fastq_files[@]}"; do
        local file="${fastq_files[$i]}"
        local current=$((i + 1))
        show_progress "$current" "${#fastq_files[@]}" "Processing $(basename "$file")"
        
        if ! run_fastqc "$file"; then
            failed_files=$((failed_files + 1))
            warning "Failed to process: $file"
        fi
    done
    
    echo "" # New line after progress indicator
    
    if [[ $failed_files -gt 0 ]]; then
        warning "$failed_files file(s) failed to process"
    fi
    
    # Skip analysis and summary for dry run
    if [[ "$DRY_RUN" == "true" ]]; then
        log "Dry run completed. No analysis performed."
        exit 0
    fi
    
    # Analyze results
    log "Analyzing FastQC results..."
    for fastqc_dir in "$OUTPUT_DIR"/*_fastqc; do
        if [[ -d "$fastqc_dir" ]]; then
            analyze_results "$fastqc_dir"
        fi
    done
    
    # Generate summary
    generate_summary
    
    # Open results (only if not running in batch mode)
    if [[ -t 1 ]]; then # Check if stdout is a terminal
        open_results
    else
        log "Batch mode detected - skipping automatic file opening"
    fi
    
    # Cleanup
    cleanup
    
    if [[ $failed_files -eq 0 ]]; then
        log "✓ Quality Control Analysis Completed Successfully!"
    else
        warning "Quality Control Analysis completed with $failed_files error(s)"
    fi
    log "=============================================="
    log "Review the HTML reports and summary file for detailed results"
    log "Summary file: $OUTPUT_DIR/quality_control_summary.txt"
}

# Handle interrupts gracefully
trap 'error "Script interrupted by user"; cleanup; exit 1' INT TERM

# Run main function
main "$@"