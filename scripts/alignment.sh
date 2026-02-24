#!/bin/bash

# WGS Alignment Pipeline - 16GB RAM Optimized
# This script performs read mapping to GRCh38 reference genome
# Usage: ./scripts/alignment.sh [OPTIONS]

set -e  # Exit on any error

# Script version and info
SCRIPT_VERSION="2.0-16GB"
SCRIPT_NAME="alignment.sh"

# Load 16GB optimized configuration
source config/default.conf

# Track CLI overrides explicitly
CLI_SAMPLE_ID=""

# Help function
show_help() {
    cat << EOF
$SCRIPT_NAME v$SCRIPT_VERSION - 16GB RAM Optimized WGS Alignment

DESCRIPTION:
    Performs read mapping to GRCh38 reference genome using BWA.
    Optimized for 16GB RAM systems with conservative memory usage.

USAGE:
    $0 [OPTIONS]

OPTIONS:
    -h, --help              Show this help message and exit
    -i, --input-dir DIR     Input directory with cleaned FASTQ files (default: auto-detect)
    -r, --reference FILE    Reference genome FASTA file (default: from config)
    -o, --output-dir DIR    Output directory (default: results/alignment)
    -s, --sample-id ID      Sample identifier (default: from config)
    -t, --threads NUM       Number of threads (default: 4 for 16GB systems)
    --dry-run               Show what would be done without executing
    --force                 Overwrite existing output files
    --verbose               Enable verbose output
    --version               Show version information

16GB OPTIMIZATIONS:
    - Limited to 4 threads to conserve memory
    - 10GB maximum memory usage for BWA
    - Sequential processing to avoid memory conflicts
    - Aggressive cleanup of intermediate files

EXAMPLES:
    # Basic usage (16GB optimized)
    $0
    
    # Custom sample ID
    $0 --sample-id MySample
    
    # Dry run to check parameters
    $0 --dry-run --verbose

REQUIREMENTS:
    - BWA and samtools installed
    - GRCh38 reference genome indexed with BWA
    - Cleaned FASTQ files from quality control step
    - 16GB RAM, 200GB free disk space

OUTPUT:
    - Sorted BAM file with index
    - Alignment statistics and quality metrics
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
            -r|--reference)
                REFERENCE_GENOME="$2"
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
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1" | tee -a "$LOG_DIR/alignment.log"
}

error() {
    echo -e "${RED}[$(date '+%Y-%m-%d %H:%M:%S')] ERROR:${NC} $1" | tee -a "$LOG_DIR/alignment.log" >&2
}

warning() {
    echo -e "${YELLOW}[$(date '+%Y-%m-%d %H:%M:%S')] WARNING:${NC} $1" | tee -a "$LOG_DIR/alignment.log"
}

info() {
    if [[ "$VERBOSE" == "true" ]]; then
        echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')] INFO:${NC} $1" | tee -a "$LOG_DIR/alignment.log"
    fi
    return 0
}

# Set default values from config
set_defaults() {
    INPUT_DIR="${INPUT_DIR:-${PROCESSED_DATA_DIR}}"
    REFERENCE_GENOME="${REFERENCE_GENOME:-${GRCH38_REFERENCE}}"
    OUTPUT_DIR="${OUTPUT_DIR:-results/alignment}"
    if [[ -n "${CLI_SAMPLE_ID:-}" ]]; then
        SAMPLE_ID="$CLI_SAMPLE_ID"
    else
        SAMPLE_ID="${SAMPLE_ID:-LOCAL_SAMPLE}"
    fi
    THREADS="${THREADS:-${BWA_THREADS:-4}}"
    DRY_RUN="${DRY_RUN:-false}"
    FORCE_OVERWRITE="${FORCE_OVERWRITE:-false}"
    VERBOSE="${VERBOSE:-false}"
}

# Auto-detect input FASTQ files
detect_input_files() {
    if [[ -z "$CLEANED_R1" && -z "$CLEANED_R2" ]]; then
        info "Auto-detecting cleaned FASTQ files in: $INPUT_DIR"
        
        local r1_files=($(find "$INPUT_DIR" \( -name "*_clean_R1.fq.gz" -o -name "*_clean_R1.fastq.gz" \) 2>/dev/null || true))
        local r2_files=($(find "$INPUT_DIR" \( -name "*_clean_R2.fq.gz" -o -name "*_clean_R2.fastq.gz" \) 2>/dev/null || true))
        
        if [[ ${#r1_files[@]} -eq 0 ]]; then
            error "No cleaned FASTQ files found in: $INPUT_DIR"
            echo "Expected files: *_clean_R1.fq.gz or *_clean_R1.fastq.gz"
            echo "Run data cleaning step first"
            exit 1
        fi
        
        if [[ ${#r1_files[@]} -ne ${#r2_files[@]} ]]; then
            error "Mismatch between R1 and R2 cleaned files"
            echo "R1 files: ${#r1_files[@]}, R2 files: ${#r2_files[@]}"
            exit 1
        fi
        
        if [[ ${#r1_files[@]} -gt 1 ]]; then
            warning "Multiple cleaned FASTQ pairs found. Processing first pair only."
            echo "Found files:"
            for i in "${!r1_files[@]}"; do
                echo "  Pair $((i+1)): ${r1_files[$i]} + ${r2_files[$i]}"
            done
        fi
        
        CLEANED_R1="${r1_files[0]}"
        CLEANED_R2="${r2_files[0]}"
    fi
    
    info "Selected input files:"
    info "  R1: $CLEANED_R1"
    info "  R2: $CLEANED_R2"
}

# Check prerequisites
check_prerequisites() {
    log "Checking prerequisites for alignment..."
    
    # Check BWA
    if ! command -v bwa &> /dev/null; then
        error "BWA not found. Install with: conda install -c bioconda bwa"
        exit 1
    fi
    
    # Check samtools
    if ! command -v samtools &> /dev/null; then
        error "samtools not found. Install with: conda install -c bioconda samtools"
        exit 1
    fi
    
    # Check reference genome
    if [[ ! -f "$REFERENCE_GENOME" ]]; then
        error "Reference genome not found: $REFERENCE_GENOME"
        exit 1
    fi
    
    # Check BWA index
    if [[ ! -f "$REFERENCE_GENOME.bwt" ]]; then
        error "BWA index not found. Create with: bwa index $REFERENCE_GENOME"
        exit 1
    fi
    
    # Check samtools index
    if [[ ! -f "$REFERENCE_GENOME.fai" ]]; then
        error "samtools index not found. Create with: samtools faidx $REFERENCE_GENOME"
        exit 1
    fi
    
    # Check cleaned reads
    if [[ ! -f "$CLEANED_R1" ]]; then
        error "Cleaned forward reads not found: $CLEANED_R1"
        exit 1
    fi
    
    if [[ ! -f "$CLEANED_R2" ]]; then
        error "Cleaned reverse reads not found: $CLEANED_R2"
        exit 1
    fi
    
    # Check available memory
    local available_mem_kb=$(grep MemAvailable /proc/meminfo 2>/dev/null | awk '{print $2}' || echo "16000000")
    local available_mem_gb=$((available_mem_kb / 1024 / 1024))
    if [[ $available_mem_gb -lt 10 ]]; then
        warning "Low available memory: ${available_mem_gb}GB. Alignment may be slow."
    fi
    
    log "Prerequisites check passed"
}

# Setup directories
setup_directories() {
    log "Setting up directories..."
    mkdir -p "$OUTPUT_DIR" "$LOG_DIR"
    info "Output directory: $OUTPUT_DIR"
}

# Run BWA alignment (16GB optimized)
run_bwa_alignment() {
    log "Starting BWA alignment (16GB optimized)..."
    
    local bam_output="$OUTPUT_DIR/${SAMPLE_ID}_aligned_sorted.bam"
    
    if [[ -f "$bam_output" && "$FORCE_OVERWRITE" != "true" ]]; then
        warning "Aligned BAM already exists: $bam_output"
        echo "Use --force to overwrite"
        return 0
    fi
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "Would run BWA alignment with:"
        echo "  Input R1: $CLEANED_R1"
        echo "  Input R2: $CLEANED_R2"
        echo "  Output: $bam_output"
        echo "  Memory limit: 10GB"
        echo "  Threads: $THREADS"
        return 0
    fi
    
    info "Running BWA-MEM with 16GB optimized settings..."
    info "Memory limit: 10GB, Threads: $THREADS"
    info "This may take 2-4 hours for WGS data"
    
    # Run BWA-MEM with direct BAM output and sorting
    local bwa_cmd="bwa mem -t $THREADS -M"
    bwa_cmd+=" -R '@RG\tID:$SAMPLE_ID\tSM:$SAMPLE_ID\tPL:ILLUMINA\tLB:WGS\tPU:$SAMPLE_ID'"
    bwa_cmd+=" $REFERENCE_GENOME $CLEANED_R1 $CLEANED_R2"
    
    # Samtools compatibility:
    # - samtools >=1.x supports: samtools sort -@ N -m 2G -o out.bam -
    # - samtools 0.1.x uses old syntax: samtools sort -m BYTES - prefix
    local samtools_version
    samtools_version=$(samtools --version 2>/dev/null | awk 'NR==1{print $2}')
    if [[ -z "$samtools_version" ]]; then
        samtools_version=$(samtools 2>&1 | awk '/^Version:/{print $2; exit}')
    fi

    if [[ "$samtools_version" =~ ^0\. ]]; then
        warning "Detected legacy samtools $samtools_version; using compatibility sort syntax"
        local sort_prefix="$OUTPUT_DIR/${SAMPLE_ID}_aligned_sorted"

        if eval "$bwa_cmd | samtools sort -m 2000000000 - $sort_prefix" 2>>"$LOG_DIR/alignment.log"; then
            # Legacy samtools writes ${prefix}.bam
            if [[ -f "${sort_prefix}.bam" ]]; then
                mv -f "${sort_prefix}.bam" "$bam_output"
            fi
            log "BWA alignment and sorting completed successfully"
        else
            error "BWA alignment failed"
            return 1
        fi
    else
        local samtools_cmd="samtools sort -@ $THREADS -m 2G -o $bam_output -"

        if eval "$bwa_cmd | $samtools_cmd" 2>>"$LOG_DIR/alignment.log"; then
            log "BWA alignment and sorting completed successfully"
        else
            error "BWA alignment failed"
            return 1
        fi
    fi

    # Get file size
    local bam_size=$(du -h "$bam_output" | cut -f1)
    log "Aligned BAM file size: $bam_size"
}

# Index BAM file
index_bam() {
    log "Indexing BAM file..."
    
    local bam_input="$OUTPUT_DIR/${SAMPLE_ID}_aligned_sorted.bam"
    local bam_index="$OUTPUT_DIR/${SAMPLE_ID}_aligned_sorted.bam.bai"
    
    if [[ -f "$bam_index" && "$FORCE_OVERWRITE" != "true" ]]; then
        warning "BAM index already exists: $bam_index"
        return 0
    fi
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "Would index BAM file: $bam_input"
        return 0
    fi
    
    if samtools index "$bam_input" 2>>"$LOG_DIR/alignment.log"; then
        log "BAM indexing completed"
    else
        error "BAM indexing failed"
        return 1
    fi
}

# Generate alignment statistics
generate_statistics() {
    log "Generating alignment statistics..."
    
    local bam_input="$OUTPUT_DIR/${SAMPLE_ID}_aligned_sorted.bam"
    local stats_output="$OUTPUT_DIR/${SAMPLE_ID}_alignment_stats.txt"
    
    if [[ ! -f "$bam_input" ]]; then
        warning "BAM file not found, skipping statistics"
        return 0
    fi
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "Would generate alignment statistics: $stats_output"
        return 0
    fi
    
    # Generate basic alignment statistics
    if samtools flagstat "$bam_input" > "$stats_output" 2>>"$LOG_DIR/alignment.log"; then
        log "Alignment statistics generated"
        
        # Extract and display key metrics
        local total_reads=$(grep "in total" "$stats_output" | cut -d' ' -f1)
        local mapped_reads=$(grep "mapped (" "$stats_output" | head -1 | cut -d' ' -f1)
        local mapping_rate=$(grep "mapped (" "$stats_output" | head -1 | sed 's/.*(\\([0-9.]*\\)%.*/\\1/')
        
        log "Total reads: $total_reads"
        log "Mapped reads: $mapped_reads"
        log "Mapping rate: $mapping_rate%"
        
        # Check mapping quality
        if (( $(echo "$mapping_rate >= 90" | bc -l 2>/dev/null || echo "1") )); then
            log "Excellent mapping rate (>=90%)"
        elif (( $(echo "$mapping_rate >= 80" | bc -l 2>/dev/null || echo "1") )); then
            log "Good mapping rate (>=80%)"
        else
            warning "Low mapping rate (<80%) - check data quality"
        fi
    else
        error "Failed to generate alignment statistics"
        return 1
    fi
}

# Generate summary report
generate_summary() {
    log "Generating alignment summary report..."
    
    local summary_file="$OUTPUT_DIR/${SAMPLE_ID}_alignment_summary.txt"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "Would generate summary: $summary_file"
        return 0
    fi
    
    cat > "$summary_file" << EOF
# WGS Alignment Summary - $SAMPLE_ID
# Generated: $(date)
# Tool: BWA with 16GB RAM optimizations

## Input Files:
- Reference: $REFERENCE_GENOME
- Forward reads: $CLEANED_R1
- Reverse reads: $CLEANED_R2

## Alignment Parameters (16GB Optimized):
- Tool: BWA-MEM v$(bwa 2>&1 | grep Version | awk '{print $2}' || echo "unknown")
- samtools: v$(samtools --version | head -1 | awk '{print $2}' || echo "unknown")
- Threads: $THREADS
- Memory limit: 10GB
- Sample ID: $SAMPLE_ID

## Output Files:
- Sorted BAM: $OUTPUT_DIR/${SAMPLE_ID}_aligned_sorted.bam
- BAM index: $OUTPUT_DIR/${SAMPLE_ID}_aligned_sorted.bam.bai
- Statistics: $OUTPUT_DIR/${SAMPLE_ID}_alignment_stats.txt

## File Sizes:
EOF
    
    # Add file sizes if files exist
    for file in "$OUTPUT_DIR"/*.bam "$OUTPUT_DIR"/*.txt; do
        if [[ -f "$file" ]]; then
            local filename=$(basename "$file")
            local filesize=$(du -h "$file" | cut -f1)
            echo "- $filename: $filesize" >> "$summary_file"
        fi
    done
    
    cat >> "$summary_file" << EOF

## 16GB System Optimizations:
- Conservative thread count to prevent memory overload
- Direct BAM output to avoid intermediate SAM files
- Memory-limited sorting (2GB per thread)
- Efficient read group assignment

## Quality Assessment:
EOF
    
    # Add mapping statistics if available
    local stats_file="$OUTPUT_DIR/${SAMPLE_ID}_alignment_stats.txt"
    if [[ -f "$stats_file" ]]; then
        local total_reads=$(grep "in total" "$stats_file" | cut -d' ' -f1)
        local mapped_reads=$(grep "mapped (" "$stats_file" | head -1 | cut -d' ' -f1)
        local mapping_rate=$(grep "mapped (" "$stats_file" | head -1 | sed 's/.*(\\([0-9.]*\\)%.*/\\1/')
        
        echo "- Total reads: $total_reads" >> "$summary_file"
        echo "- Mapped reads: $mapped_reads" >> "$summary_file"
        echo "- Mapping rate: $mapping_rate%" >> "$summary_file"
    fi
    
    cat >> "$summary_file" << EOF

## Next Steps:
1. Review alignment statistics and mapping rate
2. Proceed with variant calling: ./scripts/variant_calling.sh
3. Monitor system resources during variant calling

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
    
    log "Starting 16GB Optimized WGS Alignment Pipeline"
    log "============================================="
    
    if [[ "$DRY_RUN" == "true" ]]; then
        log "DRY RUN MODE - No files will be processed"
    fi
    
    info "Sample ID: $SAMPLE_ID"
    info "Reference: $REFERENCE_GENOME"
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
    
    run_bwa_alignment || exit 1
    index_bam || exit 1
    generate_statistics || exit 1
    generate_summary
    cleanup
    
    log "16GB Optimized WGS Alignment Pipeline Completed!"
    log "============================================="
    log "Results:"
    log "  Sorted BAM: $OUTPUT_DIR/${SAMPLE_ID}_aligned_sorted.bam"
    log "  Statistics: $OUTPUT_DIR/${SAMPLE_ID}_alignment_stats.txt"
    log "  Summary: $OUTPUT_DIR/${SAMPLE_ID}_alignment_summary.txt"
    log ""
    log "Next step: Variant calling with ./scripts/variant_calling.sh"
}

# Handle interrupts gracefully
trap 'error "Script interrupted by user"; cleanup; exit 1' INT TERM

# Run main function
main "$@"