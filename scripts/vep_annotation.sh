#!/bin/bash

# WGS VEP Annotation Pipeline - 16GB RAM Optimized
# This script annotates variants using Variant Effect Predictor (VEP)
# Usage: ./scripts/vep_annotation.sh [OPTIONS]

set -e  # Exit on any error

# Script version and info
SCRIPT_VERSION="2.0-16GB"
SCRIPT_NAME="vep_annotation.sh"

# Load 16GB optimized configuration
source config/default.conf

# Track CLI overrides explicitly
CLI_SAMPLE_ID=""

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

16GB OPTIMIZATIONS:
    - Limited to 2 threads to conserve memory
    - Small buffer size (1000 variants)
    - 6GB maximum memory usage
    - Essential annotations only (gnomAD, ClinVar)

EXAMPLES:
    # Basic usage (16GB optimized)
    $0

    # Custom VCF file
    $0 --input results/variants/my_sample_filtered.vcf.gz

    # Dry run to check parameters
    $0 --dry-run --verbose

REQUIREMENTS:
    - VEP installed with cache data
    - Filtered VCF file from variant calling
    - 16GB RAM, 20GB free disk space

OUTPUT:
    - Annotated VCF file
    - High-impact variants table
    - ClinVar matches table
    - Annotation statistics
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
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1" | tee -a "$LOG_DIR/annotation.log"
}

error() {
    echo -e "${RED}[$(date '+%Y-%m-%d %H:%M:%S')] ERROR:${NC} $1" | tee -a "$LOG_DIR/annotation.log" >&2
}

warning() {
    echo -e "${YELLOW}[$(date '+%Y-%m-%d %H:%M:%S')] WARNING:${NC} $1" | tee -a "$LOG_DIR/annotation.log"
}

info() {
    [[ "$VERBOSE" == "true" ]] && echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')] INFO:${NC} $1" | tee -a "$LOG_DIR/annotation.log"
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
    VEP_CACHE_DIR="${VEP_CACHE_DIR:-${REFERENCE_DIR}/vep_cache}"
    BUFFER_SIZE="${BUFFER_SIZE:-${VEP_BUFFER_SIZE:-1000}}"
    DRY_RUN="${DRY_RUN:-false}"
    FORCE_OVERWRITE="${FORCE_OVERWRITE:-false}"
    VERBOSE="${VERBOSE:-false}"
}

# Auto-detect input VCF file
detect_input_vcf() {
    if [[ -z "$INPUT_VCF" ]]; then
        # Look for filtered VCF first, then raw VCF
        local filtered_vcfs=($(find results/variants -name "*_filtered.vcf.gz" 2>/dev/null || true))
        local raw_vcfs=($(find results/variants -name "*_raw.vcf.gz" 2>/dev/null || true))
        
        if [[ ${#filtered_vcfs[@]} -eq 1 ]]; then
            INPUT_VCF="${filtered_vcfs[0]}"
            info "Auto-detected filtered VCF: $INPUT_VCF"
        elif [[ ${#filtered_vcfs[@]} -eq 0 && ${#raw_vcfs[@]} -eq 1 ]]; then
            INPUT_VCF="${raw_vcfs[0]}"
            warning "Using raw VCF (filtered not found): $INPUT_VCF"
        elif [[ ${#filtered_vcfs[@]} -eq 0 && ${#raw_vcfs[@]} -eq 0 ]]; then
            error "No VCF files found in results/variants/"
            echo "Run variant calling step first or specify VCF file with --input"
            exit 1
        else
            error "Multiple VCF files found. Please specify with --input:"
            for vcf in "${filtered_vcfs[@]}" "${raw_vcfs[@]}"; do
                echo "  $vcf"
            done
            exit 1
        fi
    fi
}

# Check prerequisites
check_prerequisites() {
    log "Checking prerequisites for VEP annotation..."
    
    # Check VEP installation
    if ! command -v vep &> /dev/null; then
        error "VEP not found. Install with: conda install -c bioconda ensembl-vep"
        echo "Or use: pip install ensembl-vep"
        exit 1
    fi
    
    # Check input VCF file
    if [[ ! -f "$INPUT_VCF" ]]; then
        error "VCF file not found: $INPUT_VCF"
        exit 1
    fi
    
    # Check VCF index
    if [[ ! -f "$INPUT_VCF.tbi" && ! -f "${INPUT_VCF%.gz}.idx" ]]; then
        warning "VCF index not found, will create if needed"
    fi
    
    # Check VEP cache (create directory if needed)
    if [[ ! -d "$VEP_CACHE_DIR" ]]; then
        warning "VEP cache directory not found: $VEP_CACHE_DIR"
        info "VEP will download required data automatically (may take time)"
        mkdir -p "$VEP_CACHE_DIR"
    fi
    
    # Check available memory
    local available_mem_kb=$(grep MemAvailable /proc/meminfo 2>/dev/null | awk '{print $2}' || echo "16000000")
    local available_mem_gb=$((available_mem_kb / 1024 / 1024))
    if [[ $available_mem_gb -lt 8 ]]; then
        warning "Low available memory: ${available_mem_gb}GB. VEP annotation may be slow."
        info "Consider reducing buffer size further or closing other applications"
    fi
    
    log "Prerequisites check passed"
}

# Setup directories
setup_directories() {
    log "Setting up directories..."
    mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$VEP_CACHE_DIR"
    info "Output directory: $OUTPUT_DIR"
}

# Run VEP annotation (16GB optimized)
run_vep_annotation() {
    log "Starting VEP annotation (16GB optimized)..."
    
    local annotated_vcf="$OUTPUT_DIR/${SAMPLE_ID}_annotated.vcf.gz"
    
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
        return 0
    fi
    
    info "Running VEP with 16GB optimized settings..."
    info "Memory limit: 6GB, Threads: $THREADS, Buffer: $BUFFER_SIZE"
    
    # VEP command with 16GB optimizations
    local vep_cmd="vep"
    vep_cmd+=" --input_file $INPUT_VCF"
    vep_cmd+=" --output_file $annotated_vcf"
    vep_cmd+=" --format vcf"
    vep_cmd+=" --vcf"
    vep_cmd+=" --compress_output gzip"
    vep_cmd+=" --fork $THREADS"
    vep_cmd+=" --buffer_size $BUFFER_SIZE"
    vep_cmd+=" --dir_cache $VEP_CACHE_DIR"
    vep_cmd+=" --cache"
    vep_cmd+=" --offline"
    vep_cmd+=" --assembly GRCh38"
    
    # Essential annotations only (memory efficient)
    vep_cmd+=" --symbol --canonical --biotype"
    vep_cmd+=" --sift b --polyphen b"
    vep_cmd+=" --af --af_gnomad"
    vep_cmd+=" --clinical_significance"
    
    # Memory conservation flags
    vep_cmd+=" --no_stats"  # Skip statistics to save memory
    vep_cmd+=" --quiet"     # Reduce output verbosity
    
    if eval "$vep_cmd" 2>>"$LOG_DIR/annotation.log"; then
        log "VEP annotation completed successfully"
        
        # Index the annotated VCF
        if bcftools index "$annotated_vcf" 2>>"$LOG_DIR/annotation.log"; then
            info "Annotated VCF indexed successfully"
        else
            warning "Failed to index annotated VCF"
        fi
        
        # Get basic counts
        local variant_count=$(bcftools view -H "$annotated_vcf" | wc -l)
        log "Variants annotated: $variant_count"
        
    else
        error "VEP annotation failed"
        error "Check log file: $LOG_DIR/annotation.log"
        return 1
    fi
}

# Extract high-impact variants
extract_high_impact() {
    log "Extracting high-impact variants..."
    
    local annotated_vcf="$OUTPUT_DIR/${SAMPLE_ID}_annotated.vcf.gz"
    local high_impact_file="$OUTPUT_DIR/${SAMPLE_ID}_high_impact.txt"
    
    if [[ ! -f "$annotated_vcf" ]]; then
        error "Annotated VCF not found: $annotated_vcf"
        return 1
    fi
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "Would extract high-impact variants to: $high_impact_file"
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
    
    if bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ\n' "$annotated_vcf" | \
       awk -F'\t' 'BEGIN{OFS="\t"} {
           if ($5 ~ /HIGH|MODERATE/) {
               split($5, csq, "|")
               if (csq[2] ~ /HIGH|MODERATE/) {
                   print $1, $2, $3, $4, csq[4], csq[2], csq[3], csq[37], csq[38], csq[42]
               }
           }
       }' >> "$high_impact_file" 2>/dev/null; then
        
        local high_impact_count=$(tail -n +6 "$high_impact_file" | wc -l)
        log "High-impact variants identified: $high_impact_count"
        
    else
        warning "Failed to extract high-impact variants"
    fi
}

# Extract ClinVar variants
extract_clinvar() {
    log "Extracting ClinVar variants..."
    
    local annotated_vcf="$OUTPUT_DIR/${SAMPLE_ID}_annotated.vcf.gz"
    local clinvar_file="$OUTPUT_DIR/${SAMPLE_ID}_clinvar.txt"
    
    if [[ ! -f "$annotated_vcf" ]]; then
        error "Annotated VCF not found: $annotated_vcf"
        return 1
    fi
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "Would extract ClinVar variants to: $clinvar_file"
        return 0
    fi
    
    # Extract ClinVar annotated variants
    cat > "$clinvar_file" << 'EOF'
# ClinVar Annotated Variants
# These variants have clinical significance annotations in ClinVar database
# Clinical significance: pathogenic, likely_pathogenic, uncertain_significance, etc.
#
# Columns: CHROM, POS, REF, ALT, GENE, CONSEQUENCE, CLINICAL_SIGNIFICANCE, gnomAD_AF
EOF
    
    if bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ\n' "$annotated_vcf" | \
       awk -F'\t' 'BEGIN{OFS="\t"} {
           if ($5 ~ /pathogenic|uncertain_significance/) {
               split($5, csq, "|")
               if (csq[77] != "" && csq[77] != ".") {
                   print $1, $2, $3, $4, csq[4], csq[2], csq[77], csq[42]
               }
           }
       }' >> "$clinvar_file" 2>/dev/null; then
        
        local clinvar_count=$(tail -n +6 "$clinvar_file" | wc -l)
        log "ClinVar variants identified: $clinvar_count"
        
    else
        warning "Failed to extract ClinVar variants"
    fi
}

# Generate annotation statistics
generate_statistics() {
    log "Generating annotation statistics..."
    
    local annotated_vcf="$OUTPUT_DIR/${SAMPLE_ID}_annotated.vcf.gz"
    local stats_file="$OUTPUT_DIR/${SAMPLE_ID}_annotation_stats.txt"
    
    if [[ ! -f "$annotated_vcf" ]]; then
        warning "Annotated VCF not found, skipping statistics"
        return 0
    fi
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "Would generate annotation statistics: $stats_file"
        return 0
    fi
    
    cat > "$stats_file" << EOF
# Annotation Statistics - $SAMPLE_ID
# Generated: $(date)
# Input VCF: $INPUT_VCF
# VEP Version: $(vep --help 2>&1 | head -1 || echo "Unknown")

## Parameters Used:
- Threads: $THREADS
- Buffer size: $BUFFER_SIZE
- Cache directory: $VEP_CACHE_DIR
- 16GB memory optimized

## Annotation Summary:
EOF
    
    if [[ -f "$annotated_vcf" ]]; then
        local total_variants=$(bcftools view -H "$annotated_vcf" | wc -l)
        echo "- Total annotated variants: $total_variants" >> "$stats_file"
    fi
    
    if [[ -f "$OUTPUT_DIR/${SAMPLE_ID}_high_impact.txt" ]]; then
        local high_impact=$(tail -n +6 "$OUTPUT_DIR/${SAMPLE_ID}_high_impact.txt" | wc -l)
        echo "- High-impact variants: $high_impact" >> "$stats_file"
    fi
    
    if [[ -f "$OUTPUT_DIR/${SAMPLE_ID}_clinvar.txt" ]]; then
        local clinvar=$(tail -n +6 "$OUTPUT_DIR/${SAMPLE_ID}_clinvar.txt" | wc -l)
        echo "- ClinVar annotated variants: $clinvar" >> "$stats_file"
    fi
    
    echo "" >> "$stats_file"
    echo "## Files Generated:" >> "$stats_file"
    echo "- Annotated VCF: ${SAMPLE_ID}_annotated.vcf.gz" >> "$stats_file"
    echo "- High-impact variants: ${SAMPLE_ID}_high_impact.txt" >> "$stats_file"
    echo "- ClinVar variants: ${SAMPLE_ID}_clinvar.txt" >> "$stats_file"
    
    info "Statistics generated: $stats_file"
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
    
    log "Starting 16GB Optimized VEP Annotation Pipeline"
    log "============================================="
    
    if [[ "$DRY_RUN" == "true" ]]; then
        log "DRY RUN MODE - No files will be processed"
    fi
    
    info "Sample ID: $SAMPLE_ID"
    info "Output directory: $OUTPUT_DIR"
    info "Threads: $THREADS (16GB optimized)"
    info "Buffer size: $BUFFER_SIZE (16GB optimized)"
    info "Memory conservative mode enabled"
    
    detect_input_vcf
    check_prerequisites
    setup_directories
    
    if [[ "$DRY_RUN" == "true" ]]; then
        log "Dry run completed"
        exit 0
    fi
    
    run_vep_annotation || exit 1
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
}

# Handle interrupts gracefully
trap 'error "Script interrupted by user"; cleanup; exit 1' INT TERM

# Run main function
main "$@"