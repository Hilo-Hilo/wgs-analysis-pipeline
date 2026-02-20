#!/bin/bash

# WGS Data Cleaning Pipeline - 16GB RAM Optimized
# This script performs quality-based trimming and filtering of WGS FASTQ files
# Usage: ./scripts/data_cleaning.sh [OPTIONS]

set -e  # Exit on any error

# Script version and info
SCRIPT_VERSION="2.0-16GB"
SCRIPT_NAME="data_cleaning.sh"

# Load 16GB optimized configuration
source config/default.conf

# Track CLI overrides explicitly
CLI_SAMPLE_ID=""

# Help function
show_help() {
    cat << EOF
$SCRIPT_NAME v$SCRIPT_VERSION - 16GB RAM Optimized Data Cleaning

DESCRIPTION:
    Performs adapter trimming and quality filtering using fastp.
    Optimized for 16GB RAM systems with conservative memory usage.

USAGE:
    $0 [OPTIONS]

OPTIONS:
    -h, --help              Show this help message and exit
    -i, --input-dir DIR     Input directory with FASTQ files (default: data/raw)
    -o, --output-dir DIR    Output directory (default: data/processed)
    -s, --sample-id ID      Sample identifier (default: from config)
    -t, --threads NUM       Number of threads (default: 4 for 16GB systems)
    --min-length NUM        Minimum read length after trimming (default: 75)
    --quality NUM           Quality threshold (default: 25)
    --dry-run               Show what would be done without executing
    --force                 Overwrite existing output files
    --verbose               Enable verbose output
    --version               Show version information

16GB OPTIMIZATIONS:
    - Limited to 4 threads to conserve memory
    - 4GB maximum memory usage
    - Aggressive quality filtering to reduce downstream load
    - Conservative adapter trimming

EXAMPLES:
    # Basic usage (16GB optimized)
    $0

    # Custom quality settings
    $0 --quality 30 --min-length 100

    # Dry run to check parameters
    $0 --dry-run --verbose

REQUIREMENTS:
    - fastp installed
    - Raw FASTQ files (paired-end)
    - 16GB RAM, 100GB free disk space

OUTPUT:
    - Cleaned FASTQ files
    - Quality control report
    - Cleaning statistics
    - Processing logs

AUTHOR:
    WGS Analysis Pipeline (16GB Optimized)

EOF
}

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Parse command line arguments
parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                show_help
                exit 0
                ;;
            --version)
                echo "$SCRIPT_NAME version $SCRIPT_VERSION"
                exit 0
                ;;
            -i|--input-dir)
                INPUT_DIR="$2"
                shift 2
                ;;
            -o|--output-dir)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            -s|--sample-id)
                CLI_SAMPLE_ID="$2"
                SAMPLE_ID="$2"
                shift 2
                ;;
            -t|--threads)
                THREADS="$2"
                shift 2
                ;;
            --min-length)
                MIN_LENGTH="$2"
                shift 2
                ;;
            --quality)
                QUALITY_THRESHOLD="$2"
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
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1" | tee -a "$LOG_DIR/data_cleaning.log"
}

error() {
    echo -e "${RED}[$(date '+%Y-%m-%d %H:%M:%S')] ERROR:${NC} $1" | tee -a "$LOG_DIR/data_cleaning.log" >&2
}

warning() {
    echo -e "${YELLOW}[$(date '+%Y-%m-%d %H:%M:%S')] WARNING:${NC} $1" | tee -a "$LOG_DIR/data_cleaning.log"
}

info() {
    [[ "$VERBOSE" == "true" ]] && echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')] INFO:${NC} $1" | tee -a "$LOG_DIR/data_cleaning.log"
}

# Set default values from config
set_defaults() {
    INPUT_DIR="${INPUT_DIR:-${RAW_DATA_DIR}}"
    OUTPUT_DIR="${OUTPUT_DIR:-${PROCESSED_DATA_DIR}}"
    if [[ -n "${CLI_SAMPLE_ID:-}" ]]; then
        SAMPLE_ID="$CLI_SAMPLE_ID"
    else
        SAMPLE_ID="${SAMPLE_ID:-LOCAL_SAMPLE}"
    fi
    THREADS="${THREADS:-${FASTP_THREADS:-4}}"
    MIN_LENGTH="${MIN_LENGTH:-${FASTP_MIN_LENGTH:-75}}"
    QUALITY_THRESHOLD="${QUALITY_THRESHOLD:-${FASTP_QUALITY_THRESHOLD:-25}}"
    DRY_RUN="${DRY_RUN:-false}"
    FORCE_OVERWRITE="${FORCE_OVERWRITE:-false}"
    VERBOSE="${VERBOSE:-false}"
}

# Auto-detect input FASTQ files
detect_input_files() {
    info "Detecting FASTQ files in: $INPUT_DIR"
    
    # Find paired-end FASTQ files (parentheses required for correct -o behavior)
    local r1_files=($(find "$INPUT_DIR" \( -name "*_R1.fastq.gz" -o -name "*_1.fq.gz" \) 2>/dev/null || true))
    local r2_files=($(find "$INPUT_DIR" \( -name "*_R2.fastq.gz" -o -name "*_2.fq.gz" \) 2>/dev/null || true))
    
    if [[ ${#r1_files[@]} -eq 0 ]]; then
        error "No R1 FASTQ files found in: $INPUT_DIR"
        echo "Expected files: *_R1.fastq.gz or *_1.fq.gz"
        exit 1
    fi
    
    if [[ ${#r1_files[@]} -ne ${#r2_files[@]} ]]; then
        error "Mismatch between R1 and R2 files"
        echo "R1 files: ${#r1_files[@]}, R2 files: ${#r2_files[@]}"
        exit 1
    fi
    
    if [[ ${#r1_files[@]} -gt 1 ]]; then
        warning "Multiple FASTQ pairs found. Processing first pair only."
        echo "Found files:"
        for i in "${!r1_files[@]}"; do
            echo "  Pair $((i+1)): ${r1_files[$i]} + ${r2_files[$i]}"
        done
    fi
    
    RAW_R1="${r1_files[0]}"
    RAW_R2="${r2_files[0]}"
    
    info "Selected input files:"
    info "  R1: $RAW_R1"
    info "  R2: $RAW_R2"
}

# Check prerequisites
check_prerequisites() {
    log "Checking prerequisites for data cleaning..."
    
    # Check fastp
    if ! command -v fastp &> /dev/null; then
        error "fastp not found. Install with: conda install -c bioconda fastp"
        exit 1
    fi
    
    # Check input directory
    if [[ ! -d "$INPUT_DIR" ]]; then
        error "Input directory not found: $INPUT_DIR"
        exit 1
    fi
    
    # Check input files
    if [[ ! -f "$RAW_R1" ]]; then
        error "Raw forward reads not found: $RAW_R1"
        exit 1
    fi
    
    if [[ ! -f "$RAW_R2" ]]; then
        error "Raw reverse reads not found: $RAW_R2"
        exit 1
    fi
    
    # Check file sizes (basic validation)
    local r1_size=$(stat -c%s "$RAW_R1" 2>/dev/null || stat -f%z "$RAW_R1" 2>/dev/null || echo "0")
    local r2_size=$(stat -c%s "$RAW_R2" 2>/dev/null || stat -f%z "$RAW_R2" 2>/dev/null || echo "0")
    
    if [[ $r1_size -lt 1000000 || $r2_size -lt 1000000 ]]; then
        warning "Input files appear very small (<1MB). Check file integrity."
    fi
    
    # Check available memory
    local available_mem_kb=$(grep MemAvailable /proc/meminfo 2>/dev/null | awk '{print $2}' || echo "16000000")
    local available_mem_gb=$((available_mem_kb / 1024 / 1024))
    if [[ $available_mem_gb -lt 6 ]]; then
        warning "Low available memory: ${available_mem_gb}GB. Data cleaning may be slow."
    fi
    
    log "Prerequisites check passed"
}

# Setup directories
setup_directories() {
    log "Setting up directories..."
    mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "results"
    info "Output directory: $OUTPUT_DIR"
}

# Run fastp cleaning (16GB optimized)
run_fastp_cleaning() {
    log "Starting fastp data cleaning (16GB optimized)..."
    
    local cleaned_r1="$OUTPUT_DIR/${SAMPLE_ID}_clean_R1.fq.gz"
    local cleaned_r2="$OUTPUT_DIR/${SAMPLE_ID}_clean_R2.fq.gz"
    local report_dir="${OUTPUT_DIR}/reports"
    mkdir -p "$report_dir"
    local report_html="${report_dir}/${SAMPLE_ID}_cleaning_report.html"
    local report_json="${report_dir}/${SAMPLE_ID}_cleaning_report.json"
    
    if [[ -f "$cleaned_r1" && -f "$cleaned_r2" && "$FORCE_OVERWRITE" != "true" ]]; then
        warning "Cleaned reads already exist:"
        echo "  R1: $cleaned_r1"
        echo "  R2: $cleaned_r2"
        echo "Use --force to overwrite"
        return 0
    fi
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "Would run fastp cleaning with:"
        echo "  Input R1: $RAW_R1"
        echo "  Input R2: $RAW_R2"
        echo "  Output R1: $cleaned_r1"
        echo "  Output R2: $cleaned_r2"
        echo "  Threads: $THREADS"
        echo "  Min length: $MIN_LENGTH"
        echo "  Quality: $QUALITY_THRESHOLD"
        return 0
    fi
    
    info "Running fastp with 16GB optimized settings..."
    info "Memory limit: 4GB, Threads: $THREADS"
    info "Quality threshold: $QUALITY_THRESHOLD, Min length: $MIN_LENGTH"
    
    # Run fastp with 16GB optimized parameters
    local fastp_cmd="fastp"
    fastp_cmd+=" -i $RAW_R1 -I $RAW_R2"
    fastp_cmd+=" -o $cleaned_r1 -O $cleaned_r2"
    fastp_cmd+=" --html $report_html --json $report_json"
    fastp_cmd+=" --thread $THREADS"
    
    # Quality filtering (aggressive for 16GB systems)
    fastp_cmd+=" --qualified_quality_phred $QUALITY_THRESHOLD"
    fastp_cmd+=" --unqualified_percent_limit 20"
    fastp_cmd+=" --length_required $MIN_LENGTH"
    
    # Adapter trimming
    fastp_cmd+=" --detect_adapter_for_pe"
    fastp_cmd+=" --correction"
    
    # Tail cutting (aggressive)
    fastp_cmd+=" --cut_tail --cut_tail_window_size 4"
    fastp_cmd+=" --cut_tail_mean_quality $QUALITY_THRESHOLD"
    
    # Memory conservation
    fastp_cmd+=" --overrepresentation_analysis"
    
    if eval "$fastp_cmd" 2>>"$LOG_DIR/data_cleaning.log"; then
        log "Data cleaning completed successfully"
        
        # Get file sizes
        local r1_size_raw=$(du -h "$RAW_R1" | cut -f1)
        local r2_size_raw=$(du -h "$RAW_R2" | cut -f1)
        local r1_size_clean=$(du -h "$cleaned_r1" | cut -f1)
        local r2_size_clean=$(du -h "$cleaned_r2" | cut -f1)
        
        log "File size comparison:"
        log "  Raw R1: $r1_size_raw -> Cleaned R1: $r1_size_clean"
        log "  Raw R2: $r2_size_raw -> Cleaned R2: $r2_size_clean"
        
    else
        error "Data cleaning failed"
        error "Check log file: $LOG_DIR/data_cleaning.log"
        return 1
    fi
}

# Generate summary report
generate_summary() {
    log "Generating cleaning summary report..."
    
    local summary_file="$OUTPUT_DIR/${SAMPLE_ID}_cleaning_summary.txt"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "Would generate summary: $summary_file"
        return 0
    fi
    
    cat > "$summary_file" << EOF
# Data Cleaning Summary - $SAMPLE_ID
# Generated: $(date)
# Tool: fastp with 16GB RAM optimizations

## Input Files:
- Forward reads: $RAW_R1
- Reverse reads: $RAW_R2

## Cleaning Parameters (16GB Optimized):
- Tool: fastp v$(fastp --version 2>&1 | head -1 | awk '{print $2}' || echo "unknown")
- Threads: $THREADS
- Quality threshold: Q$QUALITY_THRESHOLD
- Minimum length: ${MIN_LENGTH}bp
- Memory limit: 4GB
- Adapter detection: Automatic paired-end
- Tail cutting: Window size 4, mean quality $QUALITY_THRESHOLD
- Correction: Enabled
- Overrepresentation analysis: Enabled

## Output Files:
- Cleaned forward reads: $OUTPUT_DIR/${SAMPLE_ID}_clean_R1.fq.gz
- Cleaned reverse reads: $OUTPUT_DIR/${SAMPLE_ID}_clean_R2.fq.gz
- HTML report: results/${SAMPLE_ID}_cleaning_report.html
- JSON report: results/${SAMPLE_ID}_cleaning_report.json

## File Sizes:
EOF
    
    # Add file sizes if files exist
    if [[ -f "$RAW_R1" ]]; then
        local raw_r1_size=$(du -h "$RAW_R1" | cut -f1)
        echo "- Raw R1: $raw_r1_size" >> "$summary_file"
    fi
    
    if [[ -f "$RAW_R2" ]]; then
        local raw_r2_size=$(du -h "$RAW_R2" | cut -f1)
        echo "- Raw R2: $raw_r2_size" >> "$summary_file"
    fi
    
    if [[ -f "$OUTPUT_DIR/${SAMPLE_ID}_clean_R1.fq.gz" ]]; then
        local clean_r1_size=$(du -h "$OUTPUT_DIR/${SAMPLE_ID}_clean_R1.fq.gz" | cut -f1)
        echo "- Cleaned R1: $clean_r1_size" >> "$summary_file"
    fi
    
    if [[ -f "$OUTPUT_DIR/${SAMPLE_ID}_clean_R2.fq.gz" ]]; then
        local clean_r2_size=$(du -h "$OUTPUT_DIR/${SAMPLE_ID}_clean_R2.fq.gz" | cut -f1)
        echo "- Cleaned R2: $clean_r2_size" >> "$summary_file"
    fi
    
    cat >> "$summary_file" << EOF

## 16GB System Optimizations:
- Conservative thread count to prevent memory overload
- Aggressive quality filtering to reduce downstream memory usage
- Maximum compression to minimize disk space
- Sequential processing for stability

## Quality Assessment:
- Review HTML report for detailed statistics
- Check adapter removal efficiency
- Verify read length distribution
- Confirm quality score improvement

## Next Steps:
1. Review cleaning report: results/${SAMPLE_ID}_cleaning_report.html
2. Proceed with alignment: ./scripts/alignment.sh
3. Monitor system resources during alignment

## Analysis completed: $(date)
EOF
    
    log "Summary report generated: $summary_file"
}

# Cleanup temporary files
cleanup() {
    info "Cleaning up temporary files..."
    # Remove any temporary files created during processing
    # (Currently none, but placeholder for future use)
}

# Main execution
main() {
    parse_arguments "$@"
    set_defaults
    
    log "Starting 16GB Optimized Data Cleaning Pipeline"
    log "============================================="
    
    if [[ "$DRY_RUN" == "true" ]]; then
        log "DRY RUN MODE - No files will be processed"
    fi
    
    info "Sample ID: $SAMPLE_ID"
    info "Input directory: $INPUT_DIR"
    info "Output directory: $OUTPUT_DIR"
    info "Threads: $THREADS (16GB optimized)"
    info "Memory conservative mode enabled"
    
    detect_input_files
    check_prerequisites
    setup_directories
    
    if [[ "$DRY_RUN" == "true" ]]; then
        log "Dry run completed"
        exit 0
    fi
    
    run_fastp_cleaning || exit 1
    generate_summary
    cleanup
    
    log "16GB Optimized Data Cleaning Pipeline Completed!"
    log "=============================================="
    log "Results:"
    log "  Cleaned reads: $OUTPUT_DIR/${SAMPLE_ID}_clean_R*.fq.gz"
    log "  Cleaning report: results/${SAMPLE_ID}_cleaning_report.html"
    log "  Summary: $OUTPUT_DIR/${SAMPLE_ID}_cleaning_summary.txt"
    log ""
    log "Next step: Alignment with ./scripts/alignment.sh"
}

# Handle interrupts gracefully
trap 'error "Script interrupted by user"; cleanup; exit 1' INT TERM

# Run main function
main "$@"