#!/bin/bash

# CHM13 T2T WGS Mapping Pipeline
# This script performs read mapping to CHM13 T2T v2.0 reference genome
# Usage: ./scripts/chm13_mapping.sh

set -e  # Exit on any error

# Configuration
CONDA_ENV="wgs_analysis"
REFERENCE="data/reference/CHM13/chm13v2.0.fa"
CLEANED_R1="data/processed/SAMPLE001_clean_R1.fq.gz"
CLEANED_R2="data/processed/SAMPLE001_clean_R2.fq.gz"
OUTPUT_DIR="results/alignment"
LOG_DIR="logs"
SAMPLE_ID="SAMPLE001"
THREADS=8

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Logging function
log() {
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1" | tee -a "$LOG_DIR/chm13_mapping.log"
}

error() {
    echo -e "${RED}[$(date '+%Y-%m-%d %H:%M:%S')] ERROR:${NC} $1" | tee -a "$LOG_DIR/chm13_mapping.log"
}

warning() {
    echo -e "${YELLOW}[$(date '+%Y-%m-%d %H:%M:%S')] WARNING:${NC} $1" | tee -a "$LOG_DIR/chm13_mapping.log"
}

# Check prerequisites
check_prerequisites() {
    log "Checking prerequisites for CHM13 mapping..."
    
    # Check conda environment
    if [[ "$CONDA_DEFAULT_ENV" != "$CONDA_ENV" ]]; then
        error "Please activate conda environment: conda activate $CONDA_ENV"
        exit 1
    fi
    
    # Check CHM13 reference
    if [[ ! -f "$REFERENCE" ]]; then
        error "CHM13 reference not found: $REFERENCE"
        error "Please download and extract CHM13 reference genome first"
        exit 1
    fi
    
    # Check BWA index
    if [[ ! -f "$REFERENCE.bwt" ]]; then
        error "BWA index not found for CHM13 reference"
        error "Please run: bwa index $REFERENCE"
        exit 1
    fi
    
    # Check samtools index
    if [[ ! -f "$REFERENCE.fai" ]]; then
        error "samtools index not found for CHM13 reference"
        error "Please run: samtools faidx $REFERENCE"
        exit 1
    fi
    
    # Check cleaned reads
    if [[ ! -f "$CLEANED_R1" ]]; then
        error "Cleaned forward reads not found: $CLEANED_R1"
        error "Please run data cleaning with fastp first"
        exit 1
    fi
    
    if [[ ! -f "$CLEANED_R2" ]]; then
        error "Cleaned reverse reads not found: $CLEANED_R2"
        error "Please run data cleaning with fastp first"
        exit 1
    fi
    
    # Check tools
    if ! command -v bwa &> /dev/null; then
        error "BWA not found. Please install BWA."
        exit 1
    fi
    
    if ! command -v samtools &> /dev/null; then
        error "samtools not found. Please install samtools."
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
    
    # Reference genome info
    if [[ -f "$REFERENCE" ]]; then
        ref_size=$(du -h "$REFERENCE" | cut -f1)
        log "CHM13 reference: $ref_size"
    fi
    
    # Cleaned reads info
    if [[ -f "$CLEANED_R1" ]]; then
        r1_size=$(du -h "$CLEANED_R1" | cut -f1)
        log "Cleaned R1 reads: $r1_size"
    fi
    
    if [[ -f "$CLEANED_R2" ]]; then
        r2_size=$(du -h "$CLEANED_R2" | cut -f1)
        log "Cleaned R2 reads: $r2_size"
    fi
}

# Run BWA alignment
run_bwa_alignment() {
    log "Starting BWA alignment to CHM13 reference..."
    
    local sam_output="$OUTPUT_DIR/${SAMPLE_ID}_chm13_aligned.sam"
    local bwa_log="$LOG_DIR/bwa_chm13_alignment.log"
    
    # Check if output already exists
    if [[ -f "$sam_output" ]]; then
        warning "BWA alignment output already exists: $sam_output"
        log "Skipping BWA alignment step"
        return 0
    fi
    
    # Run BWA-MEM with read groups
    log "Running BWA-MEM alignment (this may take 2-4 hours)..."
    
    if bwa mem -t "$THREADS" -M \
        -R "@RG\tID:$SAMPLE_ID\tSM:$SAMPLE_ID\tPL:ILLUMINA\tLB:WGS\tPU:$SAMPLE_ID" \
        "$REFERENCE" \
        "$CLEANED_R1" \
        "$CLEANED_R2" \
        > "$sam_output" 2> "$bwa_log"; then
        log "✓ BWA alignment completed successfully"
        
        # Get alignment file size
        sam_size=$(du -h "$sam_output" | cut -f1)
        log "Alignment file size: $sam_size"
        
        return 0
    else
        error "BWA alignment failed"
        return 1
    fi
}

# Convert SAM to BAM and sort
convert_and_sort() {
    log "Converting SAM to BAM and sorting..."
    
    local sam_input="$OUTPUT_DIR/${SAMPLE_ID}_chm13_aligned.sam"
    local bam_output="$OUTPUT_DIR/${SAMPLE_ID}_chm13_sorted.bam"
    local conversion_log="$LOG_DIR/sam_to_bam_conversion.log"
    
    # Check if output already exists
    if [[ -f "$bam_output" ]]; then
        warning "Sorted BAM file already exists: $bam_output"
        log "Skipping SAM to BAM conversion"
        return 0
    fi
    
    # Convert SAM to BAM and sort
    if samtools view -Sb "$sam_input" 2> "$conversion_log" | \
       samtools sort -@ "$THREADS" -o "$bam_output" 2>> "$conversion_log"; then
        log "✓ SAM to BAM conversion and sorting completed"
        
        # Get BAM file size
        bam_size=$(du -h "$bam_output" | cut -f1)
        log "Sorted BAM file size: $bam_size"
        
        return 0
    else
        error "SAM to BAM conversion failed"
        return 1
    fi
}

# Index BAM file
index_bam() {
    log "Indexing BAM file..."
    
    local bam_input="$OUTPUT_DIR/${SAMPLE_ID}_chm13_sorted.bam"
    local bam_index="$OUTPUT_DIR/${SAMPLE_ID}_chm13_sorted.bam.bai"
    
    # Check if index already exists
    if [[ -f "$bam_index" ]]; then
        warning "BAM index already exists: $bam_index"
        log "Skipping BAM indexing"
        return 0
    fi
    
    # Index BAM file
    if samtools index "$bam_input" 2> "$LOG_DIR/bam_indexing.log"; then
        log "✓ BAM indexing completed"
        return 0
    else
        error "BAM indexing failed"
        return 1
    fi
}

# Generate alignment statistics
generate_statistics() {
    log "Generating alignment statistics..."
    
    local bam_input="$OUTPUT_DIR/${SAMPLE_ID}_chm13_sorted.bam"
    local stats_output="$OUTPUT_DIR/${SAMPLE_ID}_chm13_stats.txt"
    local detailed_stats="$OUTPUT_DIR/${SAMPLE_ID}_chm13_detailed_stats.txt"
    
    # Basic alignment statistics
    if samtools flagstat "$bam_input" > "$stats_output" 2> "$LOG_DIR/flagstat.log"; then
        log "✓ Basic alignment statistics generated"
        
        # Extract key metrics
        local total_reads=$(grep "in total" "$stats_output" | cut -d' ' -f1)
        local mapped_reads=$(grep "mapped (" "$stats_output" | head -1 | cut -d' ' -f1)
        local mapping_rate=$(grep "mapped (" "$stats_output" | head -1 | sed 's/.*(//' | sed 's/%.*//')
        
        log "Total reads: $total_reads"
        log "Mapped reads: $mapped_reads"
        log "Mapping rate: $mapping_rate%"
        
        # Check mapping quality
        if (( $(echo "$mapping_rate > 90" | bc -l) )); then
            log "✓ Excellent mapping rate (>90%)"
        elif (( $(echo "$mapping_rate > 80" | bc -l) )); then
            log "✓ Good mapping rate (>80%)"
        else
            warning "Low mapping rate (<80%) - investigate data quality"
        fi
    else
        error "Failed to generate basic alignment statistics"
        return 1
    fi
    
    # Detailed alignment statistics
    if samtools stats "$bam_input" > "$detailed_stats" 2> "$LOG_DIR/samtools_stats.log"; then
        log "✓ Detailed alignment statistics generated"
    else
        warning "Failed to generate detailed alignment statistics"
    fi
}

# Cleanup function
cleanup() {
    log "Cleaning up temporary files..."
    
    # Remove large SAM file if BAM conversion successful
    local sam_file="$OUTPUT_DIR/${SAMPLE_ID}_chm13_aligned.sam"
    local bam_file="$OUTPUT_DIR/${SAMPLE_ID}_chm13_sorted.bam"
    
    if [[ -f "$bam_file" && -f "$sam_file" ]]; then
        log "Removing temporary SAM file to save space..."
        rm "$sam_file"
        log "✓ Temporary SAM file removed"
    fi
    
    log "✓ Cleanup completed"
}

# Generate summary report
generate_summary() {
    log "Generating mapping summary report..."
    
    local summary_file="$OUTPUT_DIR/chm13_mapping_summary.txt"
    
    cat > "$summary_file" << EOF
# CHM13 T2T Mapping Summary - $SAMPLE_ID
# Analysis Date: $(date)
# Reference: CHM13 T2T v2.0

## Input Files:
- Reference: $REFERENCE
- Forward reads: $CLEANED_R1
- Reverse reads: $CLEANED_R2

## Analysis Parameters:
- BWA version: $(bwa 2>&1 | head -3 | tail -1)
- samtools version: $(samtools --version | head -1)
- Threads used: $THREADS
- Sample ID: $SAMPLE_ID

## Output Files:
- Sorted BAM: $OUTPUT_DIR/${SAMPLE_ID}_chm13_sorted.bam
- BAM index: $OUTPUT_DIR/${SAMPLE_ID}_chm13_sorted.bam.bai
- Basic stats: $OUTPUT_DIR/${SAMPLE_ID}_chm13_stats.txt
- Detailed stats: $OUTPUT_DIR/${SAMPLE_ID}_chm13_detailed_stats.txt

## File Sizes:
EOF
    
    # Add file sizes
    for file in "$OUTPUT_DIR"/*.bam "$OUTPUT_DIR"/*.txt; do
        if [[ -f "$file" ]]; then
            filename=$(basename "$file")
            filesize=$(du -h "$file" | cut -f1)
            echo "- $filename: $filesize" >> "$summary_file"
        fi
    done
    
    cat >> "$summary_file" << EOF

## Mapping Quality:
EOF
    
    # Add mapping statistics if available
    local stats_file="$OUTPUT_DIR/${SAMPLE_ID}_chm13_stats.txt"
    if [[ -f "$stats_file" ]]; then
        local total_reads=$(grep "in total" "$stats_file" | cut -d' ' -f1)
        local mapped_reads=$(grep "mapped (" "$stats_file" | head -1 | cut -d' ' -f1)
        local mapping_rate=$(grep "mapped (" "$stats_file" | head -1 | sed 's/.*(//' | sed 's/%.*//')
        
        echo "- Total reads: $total_reads" >> "$summary_file"
        echo "- Mapped reads: $mapped_reads" >> "$summary_file"
        echo "- Mapping rate: $mapping_rate%" >> "$summary_file"
    fi
    
    cat >> "$summary_file" << EOF

## Next Steps:
1. Review alignment statistics
2. Proceed with variant calling using bcftools
3. Generate coverage analysis
4. Perform quality score recalibration (optional)

## Analysis completed: $(date)
EOF
    
    log "✓ Summary report generated: $summary_file"
}

# Main execution
main() {
    log "Starting CHM13 T2T Mapping Pipeline for WGS Data"
    log "=============================================="
    
    # Run all steps
    check_prerequisites
    setup_directories
    get_file_info
    run_bwa_alignment || exit 1
    convert_and_sort || exit 1
    index_bam || exit 1
    generate_statistics || exit 1
    cleanup
    generate_summary
    
    log "CHM13 T2T Mapping Pipeline Completed Successfully!"
    log "=============================================="
    log "Review the alignment statistics and summary report"
    log "Summary file: $OUTPUT_DIR/chm13_mapping_summary.txt"
    log "Next: Proceed with variant calling or additional analysis"
}

# Handle interrupts gracefully
trap 'error "Script interrupted by user"; exit 1' INT TERM

# Run main function
main "$@"