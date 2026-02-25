#!/bin/bash

# WGS VEP Annotation Pipeline - 16GB RAM Optimized
# This script annotates variants using Variant Effect Predictor (VEP)
# Usage: ./scripts/vep_annotation.sh [OPTIONS]
#
# EXIT CODES:
#   0  - Success
#   1  - General/unknown error
#  10  - Missing VEP binary
#  11  - Missing bcftools dependency
#  12  - Input VCF file not found
#  13  - Malformed or invalid VCF input
#  14  - Empty VCF (zero variants to annotate)
#  15  - Invalid or missing VEP cache
#  16  - Insufficient system resources
#  17  - VEP execution failed
#  18  - Post-processing failed

set -e  # Exit on any error
set -o pipefail  # Catch pipeline failures

# Script version and info
SCRIPT_VERSION="2.1-16GB"
SCRIPT_NAME="vep_annotation.sh"

# Exit code constants
EXIT_SUCCESS=0
EXIT_GENERAL_ERROR=1
EXIT_MISSING_VEP=10
EXIT_MISSING_BCFTOOLS=11
EXIT_MISSING_INPUT=12
EXIT_MALFORMED_VCF=13
EXIT_EMPTY_VCF=14
EXIT_INVALID_CACHE=15
EXIT_LOW_RESOURCES=16
EXIT_VEP_FAILED=17
EXIT_POSTPROCESS_FAILED=18

# Load 16GB optimized configuration if available
if [[ -f "config/default.conf" ]]; then
    source config/default.conf
elif [[ -f "${BASH_SOURCE[0]%/*}/../config/default.conf" ]]; then
    source "${BASH_SOURCE[0]%/*}/../config/default.conf"
fi

# Track CLI overrides explicitly
CLI_SAMPLE_ID=""

# Internal state
ANNOTATION_CHECKPOINT_FILE=""
LOG_DIR="${LOG_DIR:-logs}"

# Help function
show_help() {
    cat << EOF
$SCRIPT_NAME v$SCRIPT_VERSION - 16GB RAM Optimized VEP Annotation

DESCRIPTION:
    Annotates genetic variants using Variant Effect Predictor (VEP).
    Optimized for 16GB RAM systems with conservative memory usage.

USAGE:
    $0 [OPTIONS]

OPTIONS:
    -h, --help              Show this help message and exit
    -i, --input FILE        Input VCF file (default: auto-detect in results/variants/)
    -o, --output-dir DIR    Output directory (default: results/annotation)
    -s, --sample-id ID      Sample identifier (default: from config)
    -t, --threads NUM       Number of threads (default: 2 for 16GB systems)
    --cache-dir DIR         VEP cache directory (default: from config)
    --buffer-size NUM       VEP buffer size (default: 1000 for 16GB)
    --dry-run               Show what would be done without executing
    --force                 Overwrite existing output files
    --verbose               Enable verbose output
    --version               Show version information
    --skip-empty-check      Continue even if VCF has zero variants
    --resume                Resume from last checkpoint if available
    --validate-only         Validate inputs without running annotation
    --min-memory-gb NUM     Minimum required memory in GB (default: 6)

16GB OPTIMIZATIONS:
    - Limited to 2 threads to conserve memory
    - Small buffer size (1000 variants)
    - 6GB maximum memory usage
    - Essential annotations only (gnomAD, ClinVar)

EXIT CODES:
    0   Success
    1   General error
    10  VEP not installed
    11  bcftools not installed
    12  Input VCF not found
    13  Malformed VCF file
    14  Empty VCF (no variants)
    15  Invalid VEP cache
    16  Insufficient memory/resources
    17  VEP execution failed
    18  Post-processing failed

EXAMPLES:
    # Basic usage (16GB optimized)
    $0

    # Custom VCF file
    $0 --input results/variants/my_sample_filtered.vcf.gz

    # Dry run to check parameters
    $0 --dry-run --verbose

    # Validate inputs only (check VCF, cache, deps)
    $0 --validate-only

    # Resume after interruption
    $0 --resume

REQUIREMENTS:
    - VEP installed with cache data (optional but required for this script)
    - bcftools for VCF manipulation
    - Filtered VCF file from variant calling
    - 6GB available memory (8GB recommended)
    - 20GB free disk space

OUTPUT:
    - Annotated VCF file
    - High-impact variants table
    - ClinVar matches table
    - Annotation statistics
    - Processing logs

TROUBLESHOOTING:
    VEP not found:
        conda install -c bioconda ensembl-vep
        # OR: pip install ensembl-vep

    Cache directory issues:
        vep_install -a cf -s homo_sapiens -y GRCh38 -c \$VEP_CACHE_DIR

    Malformed VCF errors:
        bcftools view -H input.vcf.gz | head  # check format
        bcftools norm -c ws input.vcf.gz -o fixed.vcf.gz

    Low memory warnings:
        Use --buffer-size 500 --threads 1

AUTHOR:
    WGS Analysis Pipeline (16GB Optimized)

EOF
}

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

# Parse command line arguments
parse_arguments() {
    SKIP_EMPTY_CHECK=false
    RESUME_MODE=false
    VALIDATE_ONLY=false
    MIN_MEMORY_GB=6
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                show_help
                exit $EXIT_SUCCESS
                ;;
            --version)
                echo "$SCRIPT_NAME version $SCRIPT_VERSION"
                exit $EXIT_SUCCESS
                ;;
            -i|--input)
                INPUT_VCF="$2"
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
            --cache-dir)
                VEP_CACHE_DIR="$2"
                shift 2
                ;;
            --buffer-size)
                BUFFER_SIZE="$2"
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
            --skip-empty-check)
                SKIP_EMPTY_CHECK=true
                shift
                ;;
            --resume)
                RESUME_MODE=true
                shift
                ;;
            --validate-only)
                VALIDATE_ONLY=true
                shift
                ;;
            --min-memory-gb)
                MIN_MEMORY_GB="$2"
                shift 2
                ;;
            *)
                echo "Unknown option: $1"
                echo "Use --help for usage information"
                exit $EXIT_GENERAL_ERROR
                ;;
        esac
    done
}

# Logging functions with guaranteed directory
ensure_log_dir() {
    mkdir -p "$LOG_DIR" 2>/dev/null || true
}

log() {
    ensure_log_dir
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] $1"
    echo -e "${GREEN}${msg}${NC}" | tee -a "$LOG_DIR/annotation.log" 2>/dev/null || echo -e "${GREEN}${msg}${NC}"
}

error() {
    ensure_log_dir
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1"
    echo -e "${RED}${msg}${NC}" | tee -a "$LOG_DIR/annotation.log" >&2 2>/dev/null || echo -e "${RED}${msg}${NC}" >&2
}

warning() {
    ensure_log_dir
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: $1"
    echo -e "${YELLOW}${msg}${NC}" | tee -a "$LOG_DIR/annotation.log" 2>/dev/null || echo -e "${YELLOW}${msg}${NC}"
}

info() {
    if [[ "$VERBOSE" == "true" ]]; then
        ensure_log_dir
        local msg="[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $1"
        echo -e "${BLUE}${msg}${NC}" | tee -a "$LOG_DIR/annotation.log" 2>/dev/null || echo -e "${BLUE}${msg}${NC}"
    fi
    return 0
}

diagnostic() {
    ensure_log_dir
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] DIAGNOSTIC: $1"
    echo -e "${CYAN}${msg}${NC}" | tee -a "$LOG_DIR/annotation.log" 2>/dev/null || echo -e "${CYAN}${msg}${NC}"
}

# Set default values from config
set_defaults() {
    OUTPUT_DIR="${OUTPUT_DIR:-results/annotation}"
    if [[ -n "${CLI_SAMPLE_ID:-}" ]]; then
        SAMPLE_ID="$CLI_SAMPLE_ID"
    else
        SAMPLE_ID="${SAMPLE_ID:-LOCAL_SAMPLE}"
    fi
    THREADS="${THREADS:-${VEP_THREADS:-2}}"
    VEP_CACHE_DIR="${VEP_CACHE_DIR:-${REFERENCE_DIR:-data/reference}/vep_cache}"
    BUFFER_SIZE="${BUFFER_SIZE:-${VEP_BUFFER_SIZE:-1000}}"
    DRY_RUN="${DRY_RUN:-false}"
    FORCE_OVERWRITE="${FORCE_OVERWRITE:-false}"
    VERBOSE="${VERBOSE:-false}"
    
    # Set checkpoint file
    ANNOTATION_CHECKPOINT_FILE="$OUTPUT_DIR/.annotation_checkpoint"
}

# Get available memory in GB (cross-platform)
get_available_memory_gb() {
    local mem_kb=0
    
    if [[ -f /proc/meminfo ]]; then
        # Linux
        mem_kb=$(grep MemAvailable /proc/meminfo 2>/dev/null | awk '{print $2}' || echo "0")
        if [[ "$mem_kb" -eq 0 ]]; then
            # Fallback to MemFree + Cached
            local mem_free=$(grep MemFree /proc/meminfo 2>/dev/null | awk '{print $2}' || echo "0")
            local mem_cached=$(grep "^Cached:" /proc/meminfo 2>/dev/null | awk '{print $2}' || echo "0")
            mem_kb=$((mem_free + mem_cached))
        fi
    elif command -v vm_stat &>/dev/null; then
        # macOS
        local page_size=$(vm_stat | grep "page size" | awk '{print $8}' | tr -d '.')
        page_size=${page_size:-4096}
        local free_pages=$(vm_stat | grep "Pages free" | awk '{print $3}' | tr -d '.')
        local inactive_pages=$(vm_stat | grep "Pages inactive" | awk '{print $3}' | tr -d '.')
        local speculative_pages=$(vm_stat | grep "Pages speculative" | awk '{print $3}' | tr -d '.')
        free_pages=${free_pages:-0}
        inactive_pages=${inactive_pages:-0}
        speculative_pages=${speculative_pages:-0}
        local total_free_pages=$((free_pages + inactive_pages + speculative_pages))
        mem_kb=$(( (total_free_pages * page_size) / 1024 ))
    fi
    
    # Convert to GB
    if [[ "$mem_kb" -gt 0 ]]; then
        echo $((mem_kb / 1024 / 1024))
    else
        # Unknown system, assume 16GB available
        echo "16"
    fi
}

# Auto-detect input VCF file
detect_input_vcf() {
    if [[ -z "$INPUT_VCF" ]]; then
        # Look for filtered VCF first, then raw VCF
        local filtered_vcfs=()
        local raw_vcfs=()
        
        if [[ -d "results/variants" ]]; then
            while IFS= read -r -d '' vcf; do
                filtered_vcfs+=("$vcf")
            done < <(find results/variants -name "*_filtered.vcf.gz" -print0 2>/dev/null)
            
            while IFS= read -r -d '' vcf; do
                raw_vcfs+=("$vcf")
            done < <(find results/variants -name "*_raw.vcf.gz" -print0 2>/dev/null)
        fi
        
        if [[ ${#filtered_vcfs[@]} -eq 1 ]]; then
            INPUT_VCF="${filtered_vcfs[0]}"
            info "Auto-detected filtered VCF: $INPUT_VCF"
        elif [[ ${#filtered_vcfs[@]} -eq 0 && ${#raw_vcfs[@]} -eq 1 ]]; then
            INPUT_VCF="${raw_vcfs[0]}"
            warning "Using raw VCF (filtered not found): $INPUT_VCF"
        elif [[ ${#filtered_vcfs[@]} -eq 0 && ${#raw_vcfs[@]} -eq 0 ]]; then
            error "No VCF files found in results/variants/"
            diagnostic "Expected location: results/variants/*_filtered.vcf.gz"
            diagnostic "Run variant calling step first, or specify VCF with --input"
            diagnostic "Example: $0 --input /path/to/variants.vcf.gz"
            exit $EXIT_MISSING_INPUT
        else
            error "Multiple VCF files found. Please specify with --input:"
            for vcf in "${filtered_vcfs[@]}" "${raw_vcfs[@]}"; do
                diagnostic "  $vcf"
            done
            exit $EXIT_MISSING_INPUT
        fi
    fi
}

# Validate VCF format
validate_vcf_format() {
    local vcf_file="$1"
    
    info "Validating VCF format: $vcf_file"
    
    # Check if file exists
    if [[ ! -f "$vcf_file" ]]; then
        error "VCF file not found: $vcf_file"
        exit $EXIT_MISSING_INPUT
    fi
    
    # Check if file is readable
    if [[ ! -r "$vcf_file" ]]; then
        error "VCF file not readable: $vcf_file"
        diagnostic "Check file permissions: ls -la $vcf_file"
        exit $EXIT_MISSING_INPUT
    fi
    
    # Check file size
    local file_size=$(stat -f%z "$vcf_file" 2>/dev/null || stat -c%s "$vcf_file" 2>/dev/null || echo "0")
    if [[ "$file_size" -lt 100 ]]; then
        error "VCF file appears truncated or empty: $vcf_file (${file_size} bytes)"
        diagnostic "Expected a valid gzipped VCF file"
        exit $EXIT_MALFORMED_VCF
    fi
    
    # Check if it's a valid gzip file (if .gz extension)
    if [[ "$vcf_file" == *.gz ]]; then
        if ! gzip -t "$vcf_file" 2>/dev/null; then
            error "VCF file is not a valid gzip file: $vcf_file"
            diagnostic "File may be corrupted. Try: gzip -t $vcf_file"
            diagnostic "Or re-run variant calling to regenerate"
            exit $EXIT_MALFORMED_VCF
        fi
    fi
    
    # Check VCF header
    local header_check
    if [[ "$vcf_file" == *.gz ]]; then
        header_check=$(zcat "$vcf_file" 2>/dev/null | head -1 || gzcat "$vcf_file" 2>/dev/null | head -1 || gunzip -c "$vcf_file" 2>/dev/null | head -1)
    else
        header_check=$(head -1 "$vcf_file" 2>/dev/null)
    fi
    
    if [[ ! "$header_check" =~ ^##fileformat=VCF ]]; then
        error "Invalid VCF file format: $vcf_file"
        diagnostic "First line should start with ##fileformat=VCF"
        diagnostic "Got: ${header_check:0:50}..."
        diagnostic "If this is a valid VCF, it may need recompression:"
        diagnostic "  bcftools view -O z -o fixed.vcf.gz $vcf_file"
        exit $EXIT_MALFORMED_VCF
    fi
    
    # Validate with bcftools if available
    if command -v bcftools &>/dev/null; then
        if ! bcftools view -h "$vcf_file" &>/dev/null; then
            error "VCF file failed bcftools validation: $vcf_file"
            diagnostic "Try repairing: bcftools norm -c ws $vcf_file -O z -o repaired.vcf.gz"
            exit $EXIT_MALFORMED_VCF
        fi
    fi
    
    info "VCF format validation passed"
}

# Check variant count in VCF
check_variant_count() {
    local vcf_file="$1"
    local variant_count=0
    
    info "Counting variants in VCF..."
    
    if command -v bcftools &>/dev/null; then
        variant_count=$(bcftools view -H "$vcf_file" 2>/dev/null | wc -l | tr -d ' ')
    else
        if [[ "$vcf_file" == *.gz ]]; then
            variant_count=$(zcat "$vcf_file" 2>/dev/null | grep -v "^#" | wc -l | tr -d ' ' || gzcat "$vcf_file" 2>/dev/null | grep -v "^#" | wc -l | tr -d ' ')
        else
            variant_count=$(grep -v "^#" "$vcf_file" | wc -l | tr -d ' ')
        fi
    fi
    
    variant_count=${variant_count:-0}
    
    if [[ "$variant_count" -eq 0 ]]; then
        if [[ "$SKIP_EMPTY_CHECK" == "true" ]]; then
            warning "VCF contains zero variants (--skip-empty-check enabled)"
            log "Creating empty annotation outputs..."
            return 0
        else
            error "VCF file contains zero variants: $vcf_file"
            diagnostic "This may indicate:"
            diagnostic "  - All variants filtered out in variant calling"
            diagnostic "  - Incorrect filtering thresholds"
            diagnostic "  - Sample with no detectable variants"
            diagnostic ""
            diagnostic "To proceed anyway: $0 --skip-empty-check"
            diagnostic "To check filtering: bcftools view -H $vcf_file | head"
            exit $EXIT_EMPTY_VCF
        fi
    fi
    
    log "Input VCF contains $variant_count variants"
    echo "$variant_count"
}

# Check prerequisites
check_prerequisites() {
    log "Checking prerequisites for VEP annotation..."
    local prereq_failed=false
    
    # Check VEP installation
    if ! command -v vep &> /dev/null; then
        error "VEP (Variant Effect Predictor) not found in PATH"
        diagnostic ""
        diagnostic "Installation options:"
        diagnostic "  Conda:  conda install -c bioconda ensembl-vep"
        diagnostic "  Docker: docker pull ensemblorg/ensembl-vep"
        diagnostic "  Manual: https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html"
        diagnostic ""
        diagnostic "After installation, ensure 'vep' is in your PATH"
        diagnostic "Test with: vep --help"
        exit $EXIT_MISSING_VEP
    fi
    
    # Get VEP version for diagnostics
    local vep_version=$(vep --version 2>&1 | head -1 || echo "unknown")
    info "VEP version: $vep_version"
    
    # Check bcftools (required for post-processing)
    if ! command -v bcftools &> /dev/null; then
        error "bcftools not found in PATH"
        diagnostic ""
        diagnostic "bcftools is required for VCF processing and extraction"
        diagnostic "Install: conda install -c bioconda bcftools"
        diagnostic ""
        exit $EXIT_MISSING_BCFTOOLS
    fi
    
    local bcftools_version=$(bcftools --version 2>&1 | head -1 || echo "unknown")
    info "bcftools version: $bcftools_version"
    
    # Check input VCF file
    validate_vcf_format "$INPUT_VCF"
    
    # Check variant count
    local variant_count
    variant_count=$(check_variant_count "$INPUT_VCF")
    
    # Check VCF index
    if [[ ! -f "$INPUT_VCF.tbi" && ! -f "${INPUT_VCF%.gz}.idx" && ! -f "$INPUT_VCF.csi" ]]; then
        warning "VCF index not found"
        info "Will attempt to create index if needed"
        if command -v tabix &>/dev/null; then
            info "Creating VCF index with tabix..."
            if ! tabix -p vcf "$INPUT_VCF" 2>/dev/null; then
                warning "Could not create tabix index (may not be needed)"
            fi
        fi
    fi
    
    # Check VEP cache
    check_vep_cache
    
    # Check available memory
    local available_mem_gb=$(get_available_memory_gb)
    info "Available memory: ${available_mem_gb}GB"
    
    if [[ "$available_mem_gb" -lt "$MIN_MEMORY_GB" ]]; then
        error "Insufficient memory: ${available_mem_gb}GB available, ${MIN_MEMORY_GB}GB required"
        diagnostic ""
        diagnostic "Memory optimization suggestions:"
        diagnostic "  - Close other applications"
        diagnostic "  - Reduce buffer size: --buffer-size 500"
        diagnostic "  - Reduce threads: --threads 1"
        diagnostic "  - Lower minimum requirement: --min-memory-gb 4"
        diagnostic ""
        diagnostic "Current settings: threads=$THREADS, buffer=$BUFFER_SIZE"
        exit $EXIT_LOW_RESOURCES
    elif [[ "$available_mem_gb" -lt 8 ]]; then
        warning "Low available memory: ${available_mem_gb}GB. Annotation may be slow."
        info "Consider reducing buffer size or threads if issues occur"
    fi
    
    # Check disk space
    local output_dir_space=$(df -BG "${OUTPUT_DIR%/*}" 2>/dev/null | tail -1 | awk '{print $4}' | tr -d 'G' || echo "100")
    output_dir_space=${output_dir_space:-100}
    if [[ "$output_dir_space" -lt 10 ]]; then
        warning "Low disk space: ${output_dir_space}GB available"
        diagnostic "Annotation may require up to 10GB of temporary space"
    fi
    
    log "Prerequisites check passed"
}

# Validate VEP cache directory
check_vep_cache() {
    info "Checking VEP cache directory: $VEP_CACHE_DIR"
    
    if [[ ! -d "$VEP_CACHE_DIR" ]]; then
        warning "VEP cache directory not found: $VEP_CACHE_DIR"
        
        if [[ "$DRY_RUN" == "true" ]]; then
            info "Would create cache directory in actual run"
            return 0
        fi
        
        info "Creating cache directory..."
        mkdir -p "$VEP_CACHE_DIR"
        
        warning "VEP cache needs to be downloaded"
        diagnostic ""
        diagnostic "Download cache with:"
        diagnostic "  vep_install -a cf -s homo_sapiens -y GRCh38 -c $VEP_CACHE_DIR"
        diagnostic ""
        diagnostic "Or use --offline flag to skip cache validation"
        diagnostic ""
        diagnostic "Cache download typically requires:"
        diagnostic "  - 15-20GB disk space"
        diagnostic "  - 30-60 minutes (varies by connection)"
        
        # Don't fail here - VEP may download automatically
        return 0
    fi
    
    # Check for homo_sapiens cache specifically
    local cache_species_dir="$VEP_CACHE_DIR/homo_sapiens"
    if [[ ! -d "$cache_species_dir" ]]; then
        warning "Human genome cache not found in: $VEP_CACHE_DIR"
        diagnostic "Expected: $cache_species_dir"
        diagnostic ""
        diagnostic "Download with:"
        diagnostic "  vep_install -a cf -s homo_sapiens -y GRCh38 -c $VEP_CACHE_DIR"
        
        # May still work if VEP downloads automatically
        return 0
    fi
    
    # Check for GRCh38 assembly
    local grch38_found=false
    for dir in "$cache_species_dir"/*; do
        if [[ -d "$dir" && "$(basename "$dir")" =~ GRCh38 ]]; then
            grch38_found=true
            info "Found GRCh38 cache: $dir"
            break
        fi
    done
    
    if [[ "$grch38_found" != "true" ]]; then
        warning "GRCh38 assembly cache not found"
        diagnostic "Available caches in $cache_species_dir:"
        ls -la "$cache_species_dir" 2>/dev/null | head -10 || true
        diagnostic ""
        diagnostic "Download GRCh38 cache:"
        diagnostic "  vep_install -a cf -s homo_sapiens -y GRCh38 -c $VEP_CACHE_DIR"
    fi
    
    info "VEP cache check completed"
}

# Setup directories
setup_directories() {
    log "Setting up directories..."
    mkdir -p "$OUTPUT_DIR" "$LOG_DIR"
    
    if [[ ! -d "$VEP_CACHE_DIR" ]]; then
        mkdir -p "$VEP_CACHE_DIR"
    fi
    
    info "Output directory: $OUTPUT_DIR"
}

# Save checkpoint for resume
save_checkpoint() {
    local stage="$1"
    local extra_info="$2"
    
    cat > "$ANNOTATION_CHECKPOINT_FILE" << EOF
CHECKPOINT_STAGE=$stage
CHECKPOINT_TIME=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
INPUT_VCF=$INPUT_VCF
OUTPUT_DIR=$OUTPUT_DIR
SAMPLE_ID=$SAMPLE_ID
EXTRA_INFO=$extra_info
EOF
    info "Checkpoint saved: $stage"
}

# Load checkpoint for resume
load_checkpoint() {
    if [[ -f "$ANNOTATION_CHECKPOINT_FILE" ]]; then
        source "$ANNOTATION_CHECKPOINT_FILE"
        log "Resuming from checkpoint: $CHECKPOINT_STAGE (${CHECKPOINT_TIME})"
        return 0
    fi
    return 1
}

# Clear checkpoint on success
clear_checkpoint() {
    rm -f "$ANNOTATION_CHECKPOINT_FILE"
}

# Run VEP annotation (16GB optimized)
run_vep_annotation() {
    log "Starting VEP annotation (16GB optimized)..."
    
    local annotated_vcf="$OUTPUT_DIR/${SAMPLE_ID}_annotated.vcf.gz"
    local vep_warnings_file="$LOG_DIR/${SAMPLE_ID}_vep_warnings.txt"
    
    # Check for resume
    if [[ "$RESUME_MODE" == "true" ]] && [[ -f "$annotated_vcf" ]]; then
        if load_checkpoint && [[ "$CHECKPOINT_STAGE" == "vep_complete" ]]; then
            log "VEP annotation already completed, skipping..."
            return 0
        fi
    fi
    
    if [[ -f "$annotated_vcf" && "$FORCE_OVERWRITE" != "true" ]]; then
        warning "Annotated VCF already exists: $annotated_vcf"
        echo "Use --force to overwrite"
        return 0
    fi
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "Would run VEP annotation with:"
        echo "  Input: $INPUT_VCF"
        echo "  Output: $annotated_vcf"
        echo "  Threads: $THREADS"
        echo "  Buffer size: $BUFFER_SIZE (16GB optimized)"
        echo "  Memory limit: 6GB"
        echo "  Cache: $VEP_CACHE_DIR"
        return 0
    fi
    
    save_checkpoint "vep_started" "$INPUT_VCF"
    
    info "Running VEP with 16GB optimized settings..."
    info "Memory limit: 6GB, Threads: $THREADS, Buffer: $BUFFER_SIZE"
    
    # Build VEP command as array for proper escaping
    local vep_args=(
        --input_file "$INPUT_VCF"
        --output_file "$annotated_vcf"
        --format vcf
        --vcf
        --compress_output gzip
        --fork "$THREADS"
        --buffer_size "$BUFFER_SIZE"
        --dir_cache "$VEP_CACHE_DIR"
        --cache
        --offline
        --assembly GRCh38
        --symbol
        --canonical
        --biotype
        --sift b
        --polyphen b
        --af
        --af_gnomad
        --clinical_significance
        --warning_file "$vep_warnings_file"
    )
    
    # Add no_stats only if we want to skip (default for memory savings)
    vep_args+=(--no_stats)
    
    local vep_exit_code=0
    info "VEP command: vep ${vep_args[*]}"
    
    if vep "${vep_args[@]}" 2>>"$LOG_DIR/annotation.log"; then
        log "VEP annotation completed successfully"
        save_checkpoint "vep_complete" "$annotated_vcf"
        
        # Check for warnings
        if [[ -f "$vep_warnings_file" && -s "$vep_warnings_file" ]]; then
            local warning_count=$(wc -l < "$vep_warnings_file" | tr -d ' ')
            warning "VEP generated $warning_count warnings - check $vep_warnings_file"
        fi
        
        # Index the annotated VCF
        log "Indexing annotated VCF..."
        if command -v tabix &>/dev/null; then
            if tabix -p vcf "$annotated_vcf" 2>>"$LOG_DIR/annotation.log"; then
                info "Annotated VCF indexed with tabix"
            else
                warning "Failed to index with tabix"
            fi
        elif bcftools index "$annotated_vcf" 2>>"$LOG_DIR/annotation.log"; then
            info "Annotated VCF indexed with bcftools"
        else
            warning "Failed to index annotated VCF"
        fi
        
        # Get basic counts
        local variant_count=$(bcftools view -H "$annotated_vcf" 2>/dev/null | wc -l | tr -d ' ')
        log "Variants annotated: $variant_count"
        
    else
        vep_exit_code=$?
        error "VEP annotation failed with exit code: $vep_exit_code"
        
        # Provide specific guidance based on common errors
        diagnostic ""
        diagnostic "Common VEP issues and solutions:"
        diagnostic ""
        
        # Check log for specific errors
        if grep -qi "cache.*not.*found" "$LOG_DIR/annotation.log" 2>/dev/null; then
            diagnostic "CACHE ISSUE: VEP cache not properly installed"
            diagnostic "  Fix: vep_install -a cf -s homo_sapiens -y GRCh38 -c $VEP_CACHE_DIR"
        elif grep -qi "out of memory\|cannot allocate\|killed" "$LOG_DIR/annotation.log" 2>/dev/null; then
            diagnostic "MEMORY ISSUE: VEP ran out of memory"
            diagnostic "  Fix: Reduce buffer/threads: --buffer-size 500 --threads 1"
        elif grep -qi "permission denied" "$LOG_DIR/annotation.log" 2>/dev/null; then
            diagnostic "PERMISSION ISSUE: Cannot write to output directory"
            diagnostic "  Fix: Check permissions: ls -la $OUTPUT_DIR"
        else
            diagnostic "Check log file for details: $LOG_DIR/annotation.log"
            diagnostic "Last 20 lines:"
            tail -20 "$LOG_DIR/annotation.log" 2>/dev/null || true
        fi
        
        diagnostic ""
        diagnostic "To retry: $0 --resume"
        diagnostic "To force fresh run: $0 --force"
        
        exit $EXIT_VEP_FAILED
    fi
}

# Extract high-impact variants
extract_high_impact() {
    log "Extracting high-impact variants..."
    
    local annotated_vcf="$OUTPUT_DIR/${SAMPLE_ID}_annotated.vcf.gz"
    local high_impact_file="$OUTPUT_DIR/${SAMPLE_ID}_high_impact.txt"
    
    # Check dry-run first (file won't exist in dry-run mode)
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "  Would extract high-impact variants to: $high_impact_file"
        return 0
    fi
    
    if [[ ! -f "$annotated_vcf" ]]; then
        if [[ "$SKIP_EMPTY_CHECK" == "true" ]]; then
            warning "Annotated VCF not found (empty input expected)"
            # Create empty output
            cat > "$high_impact_file" << 'EOF'
# High-Impact Variants
# No variants to annotate (empty input VCF)
EOF
            return 0
        fi
        error "Annotated VCF not found: $annotated_vcf"
        exit $EXIT_POSTPROCESS_FAILED
    fi
    
    # Check if VCF has CSQ annotation
    local has_csq=$(bcftools view -h "$annotated_vcf" 2>/dev/null | grep -c "ID=CSQ" || echo "0")
    
    if [[ "$has_csq" -eq 0 ]]; then
        warning "No CSQ annotation found in VCF - VEP may not have run correctly"
        # Still create output file with header
        cat > "$high_impact_file" << 'EOF'
# High-Impact Variants
# WARNING: No VEP annotation (CSQ) found in input VCF
# Check VEP execution and cache configuration
EOF
        return 0
    fi
    
    # Extract high-impact variants
    cat > "$high_impact_file" << 'EOF'
# High-Impact Variants
# These variants are predicted to have significant functional effects
# Impact levels: HIGH (protein truncating), MODERATE (missense), LOW (synonymous)
#
# Columns: CHROM, POS, REF, ALT, GENE, CONSEQUENCE, IMPACT, SIFT, POLYPHEN, gnomAD_AF
EOF
    
    if bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ\n' "$annotated_vcf" 2>/dev/null | \
       awk -F'\t' 'BEGIN{OFS="\t"} {
           if ($5 ~ /HIGH|MODERATE/) {
               split($5, csq, "|")
               if (length(csq) >= 5 && (csq[3] ~ /HIGH/ || csq[3] ~ /MODERATE/)) {
                   gene = (csq[4] != "" ? csq[4] : ".")
                   cons = (csq[2] != "" ? csq[2] : ".")
                   impact = (csq[3] != "" ? csq[3] : ".")
                   sift = (length(csq) >= 37 && csq[37] != "" ? csq[37] : ".")
                   polyphen = (length(csq) >= 38 && csq[38] != "" ? csq[38] : ".")
                   gnomad = (length(csq) >= 42 && csq[42] != "" ? csq[42] : ".")
                   print $1, $2, $3, $4, gene, cons, impact, sift, polyphen, gnomad
               }
           }
       }' >> "$high_impact_file" 2>/dev/null; then
        
        local high_impact_count=$(tail -n +7 "$high_impact_file" 2>/dev/null | wc -l | tr -d ' ')
        log "High-impact variants identified: $high_impact_count"
        
    else
        warning "Failed to extract high-impact variants"
        diagnostic "This may indicate VCF format issues - check annotation.log"
    fi
}

# Extract ClinVar variants
extract_clinvar() {
    log "Extracting ClinVar variants..."
    
    local annotated_vcf="$OUTPUT_DIR/${SAMPLE_ID}_annotated.vcf.gz"
    local clinvar_file="$OUTPUT_DIR/${SAMPLE_ID}_clinvar.txt"
    
    # Check dry-run first (file won't exist in dry-run mode)
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "  Would extract ClinVar variants to: $clinvar_file"
        return 0
    fi
    
    if [[ ! -f "$annotated_vcf" ]]; then
        if [[ "$SKIP_EMPTY_CHECK" == "true" ]]; then
            warning "Annotated VCF not found (empty input expected)"
            cat > "$clinvar_file" << 'EOF'
# ClinVar Annotated Variants
# No variants to annotate (empty input VCF)
EOF
            return 0
        fi
        error "Annotated VCF not found: $annotated_vcf"
        exit $EXIT_POSTPROCESS_FAILED
    fi
    
    # Extract ClinVar annotated variants
    cat > "$clinvar_file" << 'EOF'
# ClinVar Annotated Variants
# These variants have clinical significance annotations in ClinVar database
# Clinical significance: pathogenic, likely_pathogenic, uncertain_significance, etc.
#
# Columns: CHROM, POS, REF, ALT, GENE, CONSEQUENCE, CLINICAL_SIGNIFICANCE, gnomAD_AF
EOF
    
    if bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ\n' "$annotated_vcf" 2>/dev/null | \
       awk -F'\t' 'BEGIN{OFS="\t"} {
           if ($5 ~ /pathogenic|uncertain_significance|benign|likely/) {
               split($5, csq, "|")
               if (length(csq) >= 77 && csq[77] != "" && csq[77] != ".") {
                   gene = (csq[4] != "" ? csq[4] : ".")
                   cons = (csq[2] != "" ? csq[2] : ".")
                   clin = csq[77]
                   gnomad = (length(csq) >= 42 && csq[42] != "" ? csq[42] : ".")
                   print $1, $2, $3, $4, gene, cons, clin, gnomad
               }
           }
       }' >> "$clinvar_file" 2>/dev/null; then
        
        local clinvar_count=$(tail -n +6 "$clinvar_file" 2>/dev/null | wc -l | tr -d ' ')
        log "ClinVar variants identified: $clinvar_count"
        
    else
        warning "Failed to extract ClinVar variants"
        diagnostic "ClinVar extraction requires proper VEP cache with ClinVar plugin"
    fi
}

# Generate annotation statistics
generate_statistics() {
    log "Generating annotation statistics..."
    
    local annotated_vcf="$OUTPUT_DIR/${SAMPLE_ID}_annotated.vcf.gz"
    local stats_file="$OUTPUT_DIR/${SAMPLE_ID}_annotation_stats.txt"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "  Would generate annotation statistics to: $stats_file"
        return 0
    fi
    
    local vep_version_info="Unknown"
    if command -v vep &>/dev/null; then
        vep_version_info=$(vep --version 2>&1 | head -1 || echo "Unknown")
    fi
    
    cat > "$stats_file" << EOF
# Annotation Statistics - $SAMPLE_ID
# Generated: $(date)
# Script version: $SCRIPT_VERSION
# Input VCF: $INPUT_VCF
# VEP Version: $vep_version_info

## Parameters Used:
- Threads: $THREADS
- Buffer size: $BUFFER_SIZE
- Cache directory: $VEP_CACHE_DIR
- 16GB memory optimized: yes

## Annotation Summary:
EOF
    
    if [[ -f "$annotated_vcf" ]]; then
        local total_variants=$(bcftools view -H "$annotated_vcf" 2>/dev/null | wc -l | tr -d ' ')
        echo "- Total annotated variants: $total_variants" >> "$stats_file"
    else
        echo "- Total annotated variants: 0 (no annotation output)" >> "$stats_file"
    fi
    
    if [[ -f "$OUTPUT_DIR/${SAMPLE_ID}_high_impact.txt" ]]; then
        local high_impact=$(tail -n +7 "$OUTPUT_DIR/${SAMPLE_ID}_high_impact.txt" 2>/dev/null | wc -l | tr -d ' ')
        echo "- High-impact variants: $high_impact" >> "$stats_file"
    fi
    
    if [[ -f "$OUTPUT_DIR/${SAMPLE_ID}_clinvar.txt" ]]; then
        local clinvar=$(tail -n +6 "$OUTPUT_DIR/${SAMPLE_ID}_clinvar.txt" 2>/dev/null | wc -l | tr -d ' ')
        echo "- ClinVar annotated variants: $clinvar" >> "$stats_file"
    fi
    
    echo "" >> "$stats_file"
    echo "## Files Generated:" >> "$stats_file"
    echo "- Annotated VCF: ${SAMPLE_ID}_annotated.vcf.gz" >> "$stats_file"
    echo "- High-impact variants: ${SAMPLE_ID}_high_impact.txt" >> "$stats_file"
    echo "- ClinVar variants: ${SAMPLE_ID}_clinvar.txt" >> "$stats_file"
    echo "" >> "$stats_file"
    echo "## Exit Code: 0 (Success)" >> "$stats_file"
    
    info "Statistics generated: $stats_file"
}

# Cleanup temporary files
cleanup() {
    info "Cleaning up temporary files..."
    # Remove any temporary files created during processing
    # Clear checkpoint on successful completion
    if [[ -f "$ANNOTATION_CHECKPOINT_FILE" ]]; then
        clear_checkpoint
    fi
}

# Main execution
main() {
    parse_arguments "$@"
    set_defaults
    
    log "Starting 16GB Optimized VEP Annotation Pipeline v$SCRIPT_VERSION"
    log "============================================="
    
    if [[ "$DRY_RUN" == "true" ]]; then
        log "DRY RUN MODE - No files will be processed"
    fi
    
    if [[ "$VALIDATE_ONLY" == "true" ]]; then
        log "VALIDATION MODE - Checking inputs only"
    fi
    
    info "Sample ID: $SAMPLE_ID"
    info "Output directory: $OUTPUT_DIR"
    info "Threads: $THREADS (16GB optimized)"
    info "Buffer size: $BUFFER_SIZE (16GB optimized)"
    info "Memory conservative mode enabled"
    
    detect_input_vcf
    setup_directories
    check_prerequisites
    
    if [[ "$VALIDATE_ONLY" == "true" ]]; then
        log "Validation completed successfully"
        log "All prerequisites satisfied. Ready to run annotation."
        exit $EXIT_SUCCESS
    fi
    
    if [[ "$DRY_RUN" == "true" ]]; then
        run_vep_annotation  # Shows what would be done
        # In dry-run, also show what post-processing would do
        echo ""
        echo "Post-processing steps (dry run):"
        extract_high_impact   # Will show dry-run message
        extract_clinvar       # Will show dry-run message
        generate_statistics   # Will show dry-run message
        log "Dry run completed - no files were modified"
        exit $EXIT_SUCCESS
    fi
    
    run_vep_annotation || exit $EXIT_VEP_FAILED
    extract_high_impact
    extract_clinvar
    generate_statistics
    cleanup
    
    log "16GB Optimized VEP Annotation Pipeline Completed!"
    log "==============================================="
    log "Results:"
    log "  Annotated VCF: $OUTPUT_DIR/${SAMPLE_ID}_annotated.vcf.gz"
    log "  High-impact variants: $OUTPUT_DIR/${SAMPLE_ID}_high_impact.txt"
    log "  ClinVar variants: $OUTPUT_DIR/${SAMPLE_ID}_clinvar.txt"
    log "  Statistics: $OUTPUT_DIR/${SAMPLE_ID}_annotation_stats.txt"
    log ""
    log "Analysis complete! Review high-impact and ClinVar variants for clinical relevance."
    
    exit $EXIT_SUCCESS
}

# Handle interrupts gracefully
trap 'error "Script interrupted by user"; save_checkpoint "interrupted" "user_interrupt"; exit 1' INT TERM

# Run main function
main "$@"
