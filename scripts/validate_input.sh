#!/bin/bash

# Input Validation and Safety Utilities
# Comprehensive validation for WGS pipeline inputs

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Validation results
VALIDATION_PASSED=0
VALIDATION_WARNINGS=0
VALIDATION_ERRORS=0

# Logging functions
log_pass() {
    echo -e "${GREEN}✓${NC} $1"
    ((VALIDATION_PASSED++))
}

log_warn() {
    echo -e "${YELLOW}⚠${NC} $1"
    ((VALIDATION_WARNINGS++))
}

log_error() {
    echo -e "${RED}✗${NC} $1"
    ((VALIDATION_ERRORS++))
}

log_info() {
    echo -e "${BLUE}ℹ${NC} $1"
}

# Check if file exists and is readable
validate_file_exists() {
    local file="$1"
    local description="$2"
    
    if [[ ! -f "$file" ]]; then
        log_error "$description not found: $file"
        return 1
    fi
    
    if [[ ! -r "$file" ]]; then
        log_error "$description not readable: $file"
        return 1
    fi
    
    log_pass "$description found and readable: $file"
    return 0
}

# Check file size
validate_file_size() {
    local file="$1"
    local min_size_mb="${2:-1}"  # Default minimum 1MB
    local max_size_gb="${3:-1000}"  # Default maximum 1TB
    
    if [[ ! -f "$file" ]]; then
        return 1  # File doesn't exist, already handled elsewhere
    fi
    
    # Get file size in MB
    local size_bytes=$(stat -c%s "$file" 2>/dev/null || stat -f%z "$file" 2>/dev/null || echo "0")
    local size_mb=$((size_bytes / 1024 / 1024))
    local size_gb=$((size_mb / 1024))
    
    if [[ $size_mb -lt $min_size_mb ]]; then
        log_warn "File may be too small: $file (${size_mb}MB < ${min_size_mb}MB minimum)"
        return 1
    fi
    
    if [[ $size_gb -gt $max_size_gb ]]; then
        log_error "File too large: $file (${size_gb}GB > ${max_size_gb}GB maximum)"
        return 1
    fi
    
    # Display file size
    if [[ $size_gb -gt 0 ]]; then
        log_pass "File size OK: $file (${size_gb}GB)"
    else
        log_pass "File size OK: $file (${size_mb}MB)"
    fi
    
    return 0
}

# Validate FASTQ file format
validate_fastq_format() {
    local fastq_file="$1"
    local sample_lines="${2:-1000}"  # Check first 1000 lines
    
    if [[ ! -f "$fastq_file" ]]; then
        return 1
    fi
    
    log_info "Validating FASTQ format: $fastq_file"
    
    # Check if file is gzipped
    local cmd="head -n $sample_lines"
    if [[ "$fastq_file" =~ \.gz$ ]]; then
        cmd="zhead -n $sample_lines"
        if ! command -v zhead &> /dev/null; then
            cmd="gunzip -c | head -n $sample_lines"
        fi
    fi
    
    # Read sample lines
    local temp_sample="/tmp/fastq_sample_$$"
    if [[ "$fastq_file" =~ \.gz$ ]]; then
        gunzip -c "$fastq_file" | head -n "$sample_lines" > "$temp_sample" 2>/dev/null || {
            log_error "Cannot read gzipped FASTQ file: $fastq_file"
            return 1
        }
    else
        head -n "$sample_lines" "$fastq_file" > "$temp_sample" 2>/dev/null || {
            log_error "Cannot read FASTQ file: $fastq_file"
            return 1
        }
    fi
    
    # Validate FASTQ format
    local line_num=0
    local read_count=0
    local format_errors=0
    
    while IFS= read -r line; do
        local pos_in_read=$((line_num % 4))
        
        case $pos_in_read in
            0)  # Header line
                if [[ ! "$line" =~ ^@.+$ ]]; then
                    ((format_errors++))
                    if [[ $format_errors -le 3 ]]; then
                        log_error "Invalid header line $((line_num + 1)): $line"
                    fi
                fi
                ((read_count++))
                ;;
            1)  # Sequence line
                if [[ ! "$line" =~ ^[ATCGN]+$ ]]; then
                    ((format_errors++))
                    if [[ $format_errors -le 3 ]]; then
                        log_error "Invalid sequence line $((line_num + 1)): contains non-ATCGN characters"
                    fi
                fi
                ;;
            2)  # Plus line
                if [[ ! "$line" =~ ^\+.*$ ]]; then
                    ((format_errors++))
                    if [[ $format_errors -le 3 ]]; then
                        log_error "Invalid plus line $((line_num + 1)): $line"
                    fi
                fi
                ;;
            3)  # Quality line
                local seq_length=$(sed -n "$((line_num))p" "$temp_sample" | wc -c)
                seq_length=$((seq_length - 1))  # Remove newline
                local qual_length=$(echo "$line" | wc -c)
                qual_length=$((qual_length - 1))  # Remove newline
                
                if [[ $seq_length -ne $qual_length ]]; then
                    ((format_errors++))
                    if [[ $format_errors -le 3 ]]; then
                        log_error "Sequence/quality length mismatch at read $read_count"
                    fi
                fi
                ;;
        esac
        
        ((line_num++))
    done < "$temp_sample"
    
    rm -f "$temp_sample"
    
    # Check if we have complete reads (line count should be multiple of 4)
    if [[ $((line_num % 4)) -ne 0 ]]; then
        log_error "Incomplete FASTQ records (line count not multiple of 4)"
        ((format_errors++))
    fi
    
    if [[ $format_errors -eq 0 ]]; then
        log_pass "FASTQ format valid ($read_count reads sampled)"
        return 0
    else
        log_error "FASTQ format errors found: $format_errors"
        if [[ $format_errors -gt 3 ]]; then
            log_error "... and $((format_errors - 3)) more errors (first 3 shown)"
        fi
        return 1
    fi
}

# Validate FASTQ pair consistency
validate_fastq_pair() {
    local fastq_r1="$1"
    local fastq_r2="$2"
    local sample_reads="${3:-100}"  # Sample 100 reads for pairing check
    
    if [[ ! -f "$fastq_r1" ]] || [[ ! -f "$fastq_r2" ]]; then
        return 1
    fi
    
    log_info "Validating FASTQ pair consistency"
    
    # Get read count from both files (approximate)
    local count_cmd="wc -l"
    if [[ "$fastq_r1" =~ \.gz$ ]]; then
        local r1_lines=$(gunzip -c "$fastq_r1" | wc -l)
        local r2_lines=$(gunzip -c "$fastq_r2" | wc -l)
    else
        local r1_lines=$(wc -l < "$fastq_r1")
        local r2_lines=$(wc -l < "$fastq_r2")
    fi
    
    local r1_reads=$((r1_lines / 4))
    local r2_reads=$((r2_lines / 4))
    
    if [[ $r1_reads -ne $r2_reads ]]; then
        log_error "Read count mismatch: R1 has $r1_reads reads, R2 has $r2_reads reads"
        return 1
    else
        log_pass "Read counts match: $r1_reads paired reads"
    fi
    
    # Sample read headers to check pairing
    local temp_r1="/tmp/r1_headers_$$"
    local temp_r2="/tmp/r2_headers_$$"
    
    # Extract headers from sample reads
    if [[ "$fastq_r1" =~ \.gz$ ]]; then
        gunzip -c "$fastq_r1" | awk 'NR%4==1' | head -n "$sample_reads" > "$temp_r1"
        gunzip -c "$fastq_r2" | awk 'NR%4==1' | head -n "$sample_reads" > "$temp_r2"
    else
        awk 'NR%4==1' "$fastq_r1" | head -n "$sample_reads" > "$temp_r1"
        awk 'NR%4==1' "$fastq_r2" | head -n "$sample_reads" > "$temp_r2"
    fi
    
    # Check header pairing
    local pairing_errors=0
    local line_num=1
    
    while IFS= read -r r1_header <&3 && IFS= read -r r2_header <&4; do
        # Remove /1 and /2 suffixes and compare
        local r1_base=$(echo "$r1_header" | sed 's|/[12]$||' | sed 's| [12]:.*||')
        local r2_base=$(echo "$r2_header" | sed 's|/[12]$||' | sed 's| [12]:.*||')
        
        if [[ "$r1_base" != "$r2_base" ]]; then
            ((pairing_errors++))
            if [[ $pairing_errors -le 3 ]]; then
                log_error "Pairing mismatch at read $line_num:"
                log_error "  R1: $r1_header"
                log_error "  R2: $r2_header"
            fi
        fi
        
        ((line_num++))
    done 3< "$temp_r1" 4< "$temp_r2"
    
    rm -f "$temp_r1" "$temp_r2"
    
    if [[ $pairing_errors -eq 0 ]]; then
        log_pass "FASTQ pairing is consistent ($sample_reads reads checked)"
        return 0
    else
        log_error "FASTQ pairing errors: $pairing_errors"
        return 1
    fi
}

# Validate reference genome
validate_reference_genome() {
    local reference="$1"
    
    if [[ ! -f "$reference" ]]; then
        log_error "Reference genome not found: $reference"
        return 1
    fi
    
    log_info "Validating reference genome: $reference"
    
    # Check if it's a FASTA file
    local first_char=$(head -c 1 "$reference" 2>/dev/null)
    if [[ "$first_char" != ">" ]]; then
        log_error "Reference does not appear to be a FASTA file (doesn't start with '>')"
        return 1
    fi
    
    # Check basic FASTA format
    local line_count=0
    local header_count=0
    local seq_length=0
    
    while IFS= read -r line && [[ $line_count -lt 1000 ]]; do
        if [[ "$line" =~ ^'>'.*$ ]]; then
            ((header_count++))
        elif [[ "$line" =~ ^[ATCGN]+$ ]]; then
            seq_length=$((seq_length + ${#line}))
        elif [[ -n "$line" ]]; then
            log_error "Invalid characters in reference sequence"
            return 1
        fi
        ((line_count++))
    done < "$reference"
    
    if [[ $header_count -eq 0 ]]; then
        log_error "No FASTA headers found in reference"
        return 1
    fi
    
    log_pass "Reference genome format appears valid ($header_count contigs, ~${seq_length}bp sampled)"
    
    # Check for required index files
    local indexes_exist=true
    if [[ ! -f "$reference.fai" ]]; then
        log_warn "samtools index missing: $reference.fai"
        indexes_exist=false
    fi
    
    if [[ ! -f "$reference.bwt" ]]; then
        log_warn "BWA index missing: $reference.bwt"
        indexes_exist=false
    fi
    
    if [[ "$indexes_exist" == "true" ]]; then
        log_pass "Reference genome indexes found"
    else
        log_warn "Some reference indexes missing - will need to be created"
    fi
    
    return 0
}

# Check disk space requirements
validate_disk_space() {
    local required_gb="$1"
    local path="${2:-.}"
    
    # Get available disk space
    local available_kb=$(df "$path" | awk 'NR==2 {print $4}')
    local available_gb=$((available_kb / 1024 / 1024))
    
    if [[ $available_gb -lt $required_gb ]]; then
        log_error "Insufficient disk space: ${available_gb}GB available, ${required_gb}GB required"
        return 1
    else
        log_pass "Sufficient disk space: ${available_gb}GB available (${required_gb}GB required)"
        return 0
    fi
}

# Check memory requirements
validate_memory() {
    local required_gb="$1"
    
    # Get total system memory
    local total_mem_kb
    if [[ "$(uname)" == "Darwin" ]]; then
        # macOS
        total_mem_kb=$(sysctl -n hw.memsize)
        total_mem_kb=$((total_mem_kb / 1024))
    else
        # Linux
        total_mem_kb=$(grep MemTotal /proc/meminfo | awk '{print $2}')
    fi
    
    local total_mem_gb=$((total_mem_kb / 1024 / 1024))
    
    if [[ $total_mem_gb -lt $required_gb ]]; then
        log_error "Insufficient memory: ${total_mem_gb}GB available, ${required_gb}GB required"
        return 1
    else
        log_pass "Sufficient memory: ${total_mem_gb}GB available (${required_gb}GB required)"
        return 0
    fi
}

# Validate sample ID format
validate_sample_id() {
    local sample_id="$1"
    
    # Check for valid characters (alphanumeric, underscore, hyphen)
    if [[ ! "$sample_id" =~ ^[A-Za-z0-9_-]+$ ]]; then
        log_error "Invalid sample ID format: $sample_id (use only letters, numbers, underscore, hyphen)"
        return 1
    fi
    
    # Check length
    if [[ ${#sample_id} -lt 1 ]]; then
        log_error "Sample ID cannot be empty"
        return 1
    fi
    
    if [[ ${#sample_id} -gt 50 ]]; then
        log_error "Sample ID too long: $sample_id (maximum 50 characters)"
        return 1
    fi
    
    log_pass "Sample ID format valid: $sample_id"
    return 0
}

# Create backup before processing
create_backup() {
    local file="$1"
    local backup_dir="${2:-backups}"
    
    if [[ ! -f "$file" ]]; then
        return 1
    fi
    
    mkdir -p "$backup_dir"
    local backup_file="$backup_dir/$(basename "$file").backup.$(date +%Y%m%d_%H%M%S)"
    
    if cp "$file" "$backup_file"; then
        log_pass "Backup created: $backup_file"
        echo "$backup_file"  # Return backup path
        return 0
    else
        log_error "Failed to create backup for: $file"
        return 1
    fi
}

# Main validation function
validate_wgs_inputs() {
    local fastq_r1="$1"
    local fastq_r2="$2"
    local reference="$3"
    local sample_id="$4"
    local output_dir="$5"
    
    echo "🔍 Validating WGS Pipeline Inputs"
    echo "================================="
    
    # Reset counters
    VALIDATION_PASSED=0
    VALIDATION_WARNINGS=0
    VALIDATION_ERRORS=0
    
    # Validate required parameters
    if [[ -z "$fastq_r1" ]] || [[ -z "$fastq_r2" ]] || [[ -z "$reference" ]] || [[ -z "$sample_id" ]]; then
        log_error "Missing required parameters"
        echo "Usage: validate_wgs_inputs <R1.fastq> <R2.fastq> <reference.fa> <sample_id> [output_dir]"
        return 1
    fi
    
    # Validate sample ID
    validate_sample_id "$sample_id"
    
    # Validate input files exist
    validate_file_exists "$fastq_r1" "Forward reads (R1)"
    validate_file_exists "$fastq_r2" "Reverse reads (R2)" 
    validate_file_exists "$reference" "Reference genome"
    
    # Validate file sizes
    validate_file_size "$fastq_r1" 10 200000  # 10MB to 200GB
    validate_file_size "$fastq_r2" 10 200000  # 10MB to 200GB
    validate_file_size "$reference" 1 50000   # 1MB to 50GB
    
    # Validate FASTQ format
    validate_fastq_format "$fastq_r1"
    validate_fastq_format "$fastq_r2"
    
    # Validate FASTQ pairing
    validate_fastq_pair "$fastq_r1" "$fastq_r2"
    
    # Validate reference genome
    validate_reference_genome "$reference"
    
    # Validate system resources (16GB RAM, 400GB disk)
    validate_memory 16
    validate_disk_space 400 "${output_dir:-.}"
    
    # Summary
    echo
    echo "📋 Validation Summary"
    echo "===================="
    printf "✅ Passed: %d\n" "$VALIDATION_PASSED"
    printf "⚠️  Warnings: %d\n" "$VALIDATION_WARNINGS"
    printf "❌ Errors: %d\n" "$VALIDATION_ERRORS"
    
    if [[ $VALIDATION_ERRORS -eq 0 ]]; then
        echo
        log_pass "🎉 Input validation completed successfully!"
        if [[ $VALIDATION_WARNINGS -gt 0 ]]; then
            echo "⚠️ Some warnings were found - please review before proceeding"
        fi
        return 0
    else
        echo
        log_error "❌ Input validation failed - please fix errors before running pipeline"
        return 1
    fi
}

# Command line interface
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    # Script is being run directly, not sourced
    if [[ $# -lt 4 ]]; then
        echo "WGS Pipeline Input Validator"
        echo
        echo "Usage: $0 <R1.fastq> <R2.fastq> <reference.fa> <sample_id> [output_dir]"
        echo
        echo "Examples:"
        echo "  $0 data/raw/sample_R1.fastq.gz data/raw/sample_R2.fastq.gz data/reference/GRCh38.fa MySample"
        echo "  $0 reads_1.fq reads_2.fq ref.fasta SAMPLE001 results/"
        exit 1
    fi
    
    validate_wgs_inputs "$@"
    exit $?
fi