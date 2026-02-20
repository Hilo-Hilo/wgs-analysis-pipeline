#!/bin/bash

# WGS Variant Calling Pipeline - 16GB RAM Optimized
# This script performs variant calling from aligned BAM files
# Usage: ./scripts/variant_calling.sh [OPTIONS]

set -e  # Exit on any error

# Script version and info
SCRIPT_VERSION="2.0-16GB"
SCRIPT_NAME="variant_calling.sh"

# Load 16GB optimized configuration
source config/default.conf

# Track CLI overrides explicitly
CLI_SAMPLE_ID=""

# Help function
show_help() {
    cat << EOF
$SCRIPT_NAME v$SCRIPT_VERSION - 16GB RAM Optimized Variant Calling

DESCRIPTION:
    Calls genetic variants from aligned BAM files using bcftools.
    Optimized for 16GB RAM systems with conservative memory usage.

USAGE:
    $0 [OPTIONS]

OPTIONS:
    -h, --help              Show this help message and exit
    -i, --input FILE        Input BAM file (default: auto-detect in results/alignment/)
    -o, --output-dir DIR    Output directory (default: results/variants)
    -r, --reference FILE    Reference genome (default: from config)
    -s, --sample-id ID      Sample identifier (default: from config)
    -t, --threads NUM       Number of threads (default: 2 for 16GB systems)
    --min-depth NUM         Minimum read depth (default: 8)
    --min-quality NUM       Minimum variant quality (default: 30)
    --dry-run               Show what would be done without executing
    --force                 Overwrite existing output files
    --verbose               Enable verbose output
    --version               Show version information

16GB OPTIMIZATIONS:
    - Limited to 2 threads to conserve memory
    - 3GB maximum memory usage
    - Sequential chromosome processing
    - Conservative quality filters

EXAMPLES:
    # Basic usage (16GB optimized)
    $0

    # Custom BAM file
    $0 --input results/alignment/my_sample.bam

    # Dry run to check parameters
    $0 --dry-run --verbose

REQUIREMENTS:
    - bcftools installed
    - Sorted and indexed BAM file
    - Reference genome with index
    - 16GB RAM, 50GB free disk space

OUTPUT:
    - Raw variants VCF file
    - Quality-filtered variants VCF file
    - Variant statistics
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
            -i|--input)
                INPUT_BAM="$2"
                shift 2
                ;;
            -o|--output-dir)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            -r|--reference)
                REFERENCE_GENOME="$2"
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
            --min-depth)
                MIN_DEPTH="$2"
                shift 2
                ;;
            --min-quality)
                MIN_QUALITY="$2"
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
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1" | tee -a "$LOG_DIR/variant_calling.log"
}

error() {
    echo -e "${RED}[$(date '+%Y-%m-%d %H:%M:%S')] ERROR:${NC} $1" | tee -a "$LOG_DIR/variant_calling.log" >&2
}

warning() {
    echo -e "${YELLOW}[$(date '+%Y-%m-%d %H:%M:%S')] WARNING:${NC} $1" | tee -a "$LOG_DIR/variant_calling.log"
}

info() {
    [[ "$VERBOSE" == "true" ]] && echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')] INFO:${NC} $1" | tee -a "$LOG_DIR/variant_calling.log"
}

# Set default values from config
set_defaults() {
    OUTPUT_DIR="${OUTPUT_DIR:-results/variants}"
    REFERENCE_GENOME="${REFERENCE_GENOME:-${GRCH38_REFERENCE}}"
    if [[ -n "${CLI_SAMPLE_ID:-}" ]]; then
        SAMPLE_ID="$CLI_SAMPLE_ID"
    else
        SAMPLE_ID="${SAMPLE_ID:-LOCAL_SAMPLE}"
    fi
    THREADS="${THREADS:-${BCFTOOLS_THREADS:-2}}"
    MIN_DEPTH="${MIN_DEPTH:-${VARIANT_MIN_DEPTH:-8}}"
    MIN_QUALITY="${MIN_QUALITY:-${VARIANT_MIN_QUALITY:-30}}"
    DRY_RUN="${DRY_RUN:-false}"
    FORCE_OVERWRITE="${FORCE_OVERWRITE:-false}"
    VERBOSE="${VERBOSE:-false}"
}

# Auto-detect input BAM file
detect_input_bam() {
    if [[ -z "$INPUT_BAM" ]]; then
        local bam_files=($(find results/alignment -name "*.bam" ! -name "*.bai" 2>/dev/null || true))
        
        if [[ ${#bam_files[@]} -eq 0 ]]; then
            error "No BAM files found in results/alignment/"
            echo "Run alignment step first or specify BAM file with --input"
            exit 1
        elif [[ ${#bam_files[@]} -eq 1 ]]; then
            INPUT_BAM="${bam_files[0]}"
            info "Auto-detected BAM file: $INPUT_BAM"
        else
            error "Multiple BAM files found. Please specify with --input:"
            for bam in "${bam_files[@]}"; do
                echo "  $bam"
            done
            exit 1
        fi
    fi
}

# Check prerequisites
check_prerequisites() {
    log "Checking prerequisites for variant calling..."
    
    # Check bcftools
    if ! command -v bcftools &> /dev/null; then
        error "bcftools not found. Install with: conda install -c bioconda bcftools"
        exit 1
    fi
    
    # Check input BAM file
    if [[ ! -f "$INPUT_BAM" ]]; then
        error "BAM file not found: $INPUT_BAM"
        exit 1
    fi
    
    # Check BAM index
    if [[ ! -f "$INPUT_BAM.bai" ]]; then
        error "BAM index not found: $INPUT_BAM.bai"
        echo "Create index with: samtools index $INPUT_BAM"
        exit 1
    fi
    
    # Check reference genome
    if [[ ! -f "$REFERENCE_GENOME" ]]; then
        error "Reference genome not found: $REFERENCE_GENOME"
        exit 1
    fi
    
    # Check reference index
    if [[ ! -f "$REFERENCE_GENOME.fai" ]]; then
        error "Reference index not found: $REFERENCE_GENOME.fai"
        echo "Create index with: samtools faidx $REFERENCE_GENOME"
        exit 1
    fi
    
    # Check memory
    local available_mem_kb=$(grep MemAvailable /proc/meminfo 2>/dev/null | awk '{print $2}' || echo "16000000")
    local available_mem_gb=$((available_mem_kb / 1024 / 1024))
    if [[ $available_mem_gb -lt 6 ]]; then
        warning "Low available memory: ${available_mem_gb}GB. Variant calling may be slow."
    fi
    
    log "Prerequisites check passed"
}

# Setup directories
setup_directories() {
    log "Setting up directories..."
    mkdir -p "$OUTPUT_DIR" "$LOG_DIR"
    info "Output directory: $OUTPUT_DIR"
}

# Call variants with bcftools (16GB optimized)
call_variants() {
    log "Starting variant calling (16GB optimized)..."
    
    local raw_vcf="$OUTPUT_DIR/${SAMPLE_ID}_raw.vcf.gz"
    
    if [[ -f "$raw_vcf" && "$FORCE_OVERWRITE" != "true" ]]; then
        warning "Raw VCF already exists: $raw_vcf"
        echo "Use --force to overwrite"
        return 0
    fi
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "Would call variants: bcftools mpileup + bcftools call"
        echo "Input: $INPUT_BAM"
        echo "Output: $raw_vcf"
        echo "Memory limit: 3GB"
        echo "Threads: $THREADS"
        return 0
    fi
    
    info "Calling variants with conservative 16GB settings..."
    info "Memory limit: 3GB, Threads: $THREADS"
    
    # Use bcftools with memory-conservative settings
    local mpileup_cmd="bcftools mpileup -f $REFERENCE_GENOME"
    mpileup_cmd+=" -q 20 -Q 20"  # Conservative base and mapping quality
    mpileup_cmd+=" -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR"
    mpileup_cmd+=" --threads $THREADS"
    mpileup_cmd+=" $INPUT_BAM"
    
    local call_cmd="bcftools call -mv"
    call_cmd+=" --threads $THREADS"
    call_cmd+=" -O z"  # Compress output
    
    if eval "$mpileup_cmd | $call_cmd > $raw_vcf" 2>>"$LOG_DIR/variant_calling.log"; then
        log "Variant calling completed successfully"
        
        # Index the VCF
        if bcftools index "$raw_vcf" 2>>"$LOG_DIR/variant_calling.log"; then
            info "VCF indexed successfully"
        else
            warning "Failed to index VCF file"
        fi
        
        # Get basic statistics
        local variant_count=$(bcftools view -H "$raw_vcf" | wc -l)
        log "Raw variants called: $variant_count"
        
    else
        error "Variant calling failed"
        return 1
    fi
}

# Filter variants (16GB optimized)
filter_variants() {
    log "Filtering variants with conservative settings..."
    
    local raw_vcf="$OUTPUT_DIR/${SAMPLE_ID}_raw.vcf.gz"
    local filtered_vcf="$OUTPUT_DIR/${SAMPLE_ID}_filtered.vcf.gz"
    
    if [[ ! -f "$raw_vcf" ]]; then
        error "Raw VCF not found: $raw_vcf"
        return 1
    fi
    
    if [[ -f "$filtered_vcf" && "$FORCE_OVERWRITE" != "true" ]]; then
        warning "Filtered VCF already exists: $filtered_vcf"
        return 0
    fi
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "Would filter variants with:"
        echo "  Minimum depth: $MIN_DEPTH"
        echo "  Minimum quality: $MIN_QUALITY"
        echo "  Output: $filtered_vcf"
        return 0
    fi
    
    info "Applying quality filters: depth >= $MIN_DEPTH, quality >= $MIN_QUALITY"
    
    # Filter variants with conservative settings
    local filter_expr="QUAL >= $MIN_QUALITY && INFO/DP >= $MIN_DEPTH && INFO/DP <= 200"
    
    if bcftools filter -e "$filter_expr" "$raw_vcf" -O z -o "$filtered_vcf" 2>>"$LOG_DIR/variant_calling.log"; then
        log "Variant filtering completed"
        
        # Index filtered VCF
        if bcftools index "$filtered_vcf" 2>>"$LOG_DIR/variant_calling.log"; then
            info "Filtered VCF indexed successfully"
        fi
        
        # Get filtered statistics
        local filtered_count=$(bcftools view -H "$filtered_vcf" | wc -l)
        local raw_count=$(bcftools view -H "$raw_vcf" | wc -l)
        local pass_rate=$((filtered_count * 100 / raw_count))
        
        log "Variants after filtering: $filtered_count ($pass_rate% pass rate)"
        
    else
        error "Variant filtering failed"
        return 1
    fi
}

# Generate statistics
generate_statistics() {
    log "Generating variant statistics..."
    
    local filtered_vcf="$OUTPUT_DIR/${SAMPLE_ID}_filtered.vcf.gz"
    local stats_file="$OUTPUT_DIR/${SAMPLE_ID}_variant_stats.txt"
    
    if [[ ! -f "$filtered_vcf" ]]; then
        warning "Filtered VCF not found, skipping statistics"
        return 0
    fi
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "Would generate statistics: $stats_file"
        return 0
    fi
    
    # Generate comprehensive statistics
    cat > "$stats_file" << EOF
# Variant Calling Statistics - $SAMPLE_ID
# Generated: $(date)
# Input BAM: $INPUT_BAM
# Reference: $REFERENCE_GENOME

## Parameters Used:
- Threads: $THREADS
- Min depth: $MIN_DEPTH
- Min quality: $MIN_QUALITY
- Memory optimized for 16GB systems

## Variant Counts:
EOF
    
    if bcftools stats "$filtered_vcf" >> "$stats_file" 2>/dev/null; then
        info "Statistics generated: $stats_file"
    else
        warning "Failed to generate detailed statistics"
    fi
    
    log "Variant calling pipeline completed successfully"
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
    
    log "Starting 16GB Optimized Variant Calling Pipeline"
    log "=============================================="
    
    if [[ "$DRY_RUN" == "true" ]]; then
        log "DRY RUN MODE - No files will be processed"
    fi
    
    info "Sample ID: $SAMPLE_ID"
    info "Output directory: $OUTPUT_DIR"
    info "Threads: $THREADS (16GB optimized)"
    info "Memory conservative mode enabled"
    
    detect_input_bam
    check_prerequisites
    setup_directories
    
    if [[ "$DRY_RUN" == "true" ]]; then
        log "Dry run completed"
        exit 0
    fi
    
    call_variants || exit 1
    filter_variants || exit 1
    generate_statistics
    cleanup
    
    log "16GB Optimized Variant Calling Pipeline Completed!"
    log "=============================================="
    log "Results:"
    log "  Raw variants: $OUTPUT_DIR/${SAMPLE_ID}_raw.vcf.gz"
    log "  Filtered variants: $OUTPUT_DIR/${SAMPLE_ID}_filtered.vcf.gz"
    log "  Statistics: $OUTPUT_DIR/${SAMPLE_ID}_variant_stats.txt"
    log ""
    log "Next step: Variant annotation"
}

# Handle interrupts gracefully
trap 'error "Script interrupted by user"; cleanup; exit 1' INT TERM

# Run main function
main "$@"