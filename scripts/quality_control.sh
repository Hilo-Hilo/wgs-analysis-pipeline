#!/bin/bash

# WGS Quality Control Analysis Script
# This script performs comprehensive quality control analysis on raw FASTQ files
# Usage: ./scripts/quality_control.sh [OPTIONS]

set -e
set -o pipefail

# Script version and info
SCRIPT_VERSION="2.1"
SCRIPT_NAME="quality_control.sh"

# Exit codes
EX_USAGE=2
EX_PREREQ=10
EX_INPUT=11
EX_PAIRING=12
EX_TINY_FILE=13
EX_FASTQC_FAIL=14
EX_OUTPUT_MISSING=15

# Default configuration (16GB RAM optimized)
CONDA_ENV="wgs_analysis"
RAW_DATA_DIR="data/raw"
OUTPUT_DIR="results/quality_control"
LOG_DIR="logs"
THREADS=4
DRY_RUN=false
VERBOSE=false
FORCE_OVERWRITE=false

# Edge-case guardrails
MIN_FASTQ_BYTES=4096          # <4KB is almost certainly invalid for WGS FASTQ
WARN_SMALL_FASTQ_BYTES=1000000 # warn under 1MB

# Runtime state
FASTQ_FILES=()
CRITICAL_QUALITY_FILES=()
FASTQ_FAILURES=0

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

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

EXIT CODES:
    0   Success (warnings possible)
    2   Invalid command-line usage
    10  Missing prerequisite/tool
    11  Invalid input files or directories
    12  Broken FASTQ pairing detected
    13  Tiny/corrupt FASTQ detected
    14  FastQC execution failure
    15  Missing expected output artifacts

EXAMPLES:
    $0
    $0 --input-dir /path/to/fastq/files
    $0 --threads 16
    $0 --dry-run
    $0 --force --verbose
EOF
}

# Version function
show_version() {
    echo "$SCRIPT_NAME version $SCRIPT_VERSION"
    exit 0
}

safe_mkdir() {
    local dir="$1"
    [[ -n "$dir" ]] && mkdir -p "$dir" 2>/dev/null || true
}

log_to_file() {
    local line="$1"
    safe_mkdir "$LOG_DIR"
    echo -e "$line" >> "$LOG_DIR/quality_control.log" 2>/dev/null || true
}

# Logging functions
log() {
    local message="$1"
    local timestamp="[$(date '+%Y-%m-%d %H:%M:%S')]"
    echo -e "${GREEN}${timestamp}${NC} $message"
    log_to_file "${timestamp} $message"
}

error() {
    local message="$1"
    local timestamp="[$(date '+%Y-%m-%d %H:%M:%S')]"
    echo -e "${RED}${timestamp} ERROR:${NC} $message" >&2
    log_to_file "${timestamp} ERROR: $message"
}

warning() {
    local message="$1"
    local timestamp="[$(date '+%Y-%m-%d %H:%M:%S')]"
    echo -e "${YELLOW}${timestamp} WARNING:${NC} $message"
    log_to_file "${timestamp} WARNING: $message"
}

info() {
    local message="$1"
    local timestamp="[$(date '+%Y-%m-%d %H:%M:%S')]"
    if [[ "$VERBOSE" == "true" ]]; then
        echo -e "${BLUE}${timestamp} INFO:${NC} $message"
    fi
    log_to_file "${timestamp} INFO: $message"
}

# Progress indicator
show_progress() {
    local current=$1
    local total=$2
    local desc="$3"

    if [[ "$total" -le 0 ]]; then
        return 0
    fi

    local percent=$((current * 100 / total))
    local bars=$((percent / 5))
    local bar_text=""
    local i
    for ((i=0; i<bars; i++)); do
        bar_text+="#"
    done

    printf "\rProgress: [%-20s] %d%% %s" "$bar_text" "$percent" "$desc"
    if [[ $current -eq $total ]]; then
        echo ""
    fi
}

require_value() {
    local flag="$1"
    local value="${2:-}"
    if [[ -z "$value" || "$value" == -* ]]; then
        error "Option $flag requires a value"
        exit $EX_USAGE
    fi
}

# Parse command line arguments
parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -h|--help)
                show_help
                exit 0
                ;;
            --version)
                show_version
                ;;
            -i|--input-dir)
                require_value "$1" "${2:-}"
                RAW_DATA_DIR="$2"
                shift 2
                ;;
            -o|--output-dir)
                require_value "$1" "${2:-}"
                OUTPUT_DIR="$2"
                shift 2
                ;;
            -t|--threads)
                require_value "$1" "${2:-}"
                THREADS="$2"
                shift 2
                ;;
            -e|--env)
                require_value "$1" "${2:-}"
                CONDA_ENV="$2"
                shift 2
                ;;
            -l|--log-dir)
                require_value "$1" "${2:-}"
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
                error "Unknown option: $1"
                echo "Use --help for usage information"
                exit $EX_USAGE
                ;;
        esac
    done
}

# Validate input parameters
validate_parameters() {
    info "Validating parameters..."

    if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -le 0 ]]; then
        error "Invalid threads value: $THREADS. Must be a positive integer."
        exit $EX_USAGE
    fi

    # Convert relative paths to absolute paths when possible
    RAW_DATA_DIR=$(realpath "$RAW_DATA_DIR" 2>/dev/null || echo "$RAW_DATA_DIR")
    OUTPUT_DIR=$(realpath "$OUTPUT_DIR" 2>/dev/null || echo "$OUTPUT_DIR")
    LOG_DIR=$(realpath "$LOG_DIR" 2>/dev/null || echo "$LOG_DIR")

    info "Using input directory: $RAW_DATA_DIR"
    info "Using output directory: $OUTPUT_DIR"
    info "Using log directory: $LOG_DIR"
    info "Using threads: $THREADS"
}

get_file_size_bytes() {
    local file="$1"
    stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null || echo 0
}

collect_fastq_files() {
    FASTQ_FILES=()
    while IFS= read -r -d '' file; do
        FASTQ_FILES+=("$file")
    done < <(find "$RAW_DATA_DIR" -maxdepth 1 -type f \( -name "*.fq.gz" -o -name "*.fastq.gz" \) -print0 2>/dev/null)
}

pair_key_and_mate() {
    local basename_file="$1"

    if [[ "$basename_file" =~ ^(.+)_R1\.(fastq|fq)\.gz$ ]]; then
        echo "${BASH_REMATCH[1]}|R1"
        return 0
    elif [[ "$basename_file" =~ ^(.+)_R2\.(fastq|fq)\.gz$ ]]; then
        echo "${BASH_REMATCH[1]}|R2"
        return 0
    elif [[ "$basename_file" =~ ^(.+)_1\.(fastq|fq)\.gz$ ]]; then
        echo "${BASH_REMATCH[1]}|R1"
        return 0
    elif [[ "$basename_file" =~ ^(.+)_2\.(fastq|fq)\.gz$ ]]; then
        echo "${BASH_REMATCH[1]}|R2"
        return 0
    fi

    return 1
}

increment_pair_count() {
    local key="$1"
    local mate="$2"

    local idx=-1
    local i
    for i in "${!PAIR_KEYS[@]}"; do
        if [[ "${PAIR_KEYS[$i]}" == "$key" ]]; then
            idx=$i
            break
        fi
    done

    if [[ $idx -lt 0 ]]; then
        PAIR_KEYS+=("$key")
        PAIR_R1_COUNTS+=(0)
        PAIR_R2_COUNTS+=(0)
        idx=$((${#PAIR_KEYS[@]} - 1))
    fi

    if [[ "$mate" == "R1" ]]; then
        PAIR_R1_COUNTS[$idx]=$((PAIR_R1_COUNTS[$idx] + 1))
    else
        PAIR_R2_COUNTS[$idx]=$((PAIR_R2_COUNTS[$idx] + 1))
    fi
}

validate_pairing_and_input_sizes() {
    PAIR_KEYS=()
    PAIR_R1_COUNTS=()
    PAIR_R2_COUNTS=()
    local recognized=0
    local pairing_errors=0

    local file basename_file size_bytes pair_info key mate
    for file in "${FASTQ_FILES[@]}"; do
        basename_file=$(basename "$file")
        size_bytes=$(get_file_size_bytes "$file")

        if [[ "$size_bytes" -lt "$MIN_FASTQ_BYTES" ]]; then
            error "FASTQ file appears too small (${size_bytes} bytes): $basename_file"
            error "Action: verify download/decompression; expected WGS FASTQ files are much larger."
            return $EX_TINY_FILE
        elif [[ "$size_bytes" -lt "$WARN_SMALL_FASTQ_BYTES" ]]; then
            warning "Small FASTQ file detected (${size_bytes} bytes): $basename_file"
            warning "This may be a tiny test dataset or truncated upload."
        fi

        if ! gzip -t "$file" 2>/dev/null; then
            error "Corrupt gzip FASTQ detected: $basename_file"
            error "Action: re-transfer/re-generate this file."
            return $EX_INPUT
        fi

        if pair_info=$(pair_key_and_mate "$basename_file"); then
            recognized=$((recognized + 1))
            key="${pair_info%%|*}"
            mate="${pair_info##*|}"
            increment_pair_count "$key" "$mate"
        fi
    done

    if [[ "$recognized" -gt 0 ]]; then
        local i
        for i in "${!PAIR_KEYS[@]}"; do
            if [[ "${PAIR_R1_COUNTS[$i]}" -ne 1 || "${PAIR_R2_COUNTS[$i]}" -ne 1 ]]; then
                pairing_errors=$((pairing_errors + 1))
                error "Broken FASTQ pairing for sample '${PAIR_KEYS[$i]}': R1=${PAIR_R1_COUNTS[$i]}, R2=${PAIR_R2_COUNTS[$i]}"
            fi
        done

        if [[ "$pairing_errors" -gt 0 ]]; then
            error "FASTQ pairing validation failed."
            error "Action: ensure one and only one mate per pair (*_R1/*_R2 or *_1/*_2)."
            return $EX_PAIRING
        fi
    fi

    return 0
}

check_disk_space() {
    local check_path="$OUTPUT_DIR"
    while [[ ! -d "$check_path" && "$check_path" != "/" ]]; do
        check_path=$(dirname "$check_path")
    done

    local available_space
    available_space=$(df "$check_path" 2>/dev/null | awk 'NR==2 {print $4}' || echo "0")

    if [[ "$available_space" -lt 1048576 ]]; then # <1GB
        warning "Low disk space available (<1GB) at $check_path. FastQC may fail."
    fi
}

# Check prerequisites
check_prerequisites() {
    log "Checking prerequisites..."

    if ! command -v fastqc >/dev/null 2>&1; then
        error "FastQC not found in PATH."
        echo "Install FastQC (for example):"
        echo "  conda install -c bioconda fastqc"
        echo "Or ensure your environment/container PATH includes fastqc."
        exit $EX_PREREQ
    fi

    if [[ ! -d "$RAW_DATA_DIR" ]]; then
        error "Input directory not found: $RAW_DATA_DIR"
        echo "Create the directory or specify a different path with --input-dir"
        exit $EX_INPUT
    fi

    collect_fastq_files

    if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
        error "No FASTQ files found in: $RAW_DATA_DIR"
        echo "Expected files with extensions: .fq.gz or .fastq.gz"
        echo "Available files in directory:"
        ls -la "$RAW_DATA_DIR" 2>/dev/null || echo "  Directory is empty or unreadable"
        exit $EX_INPUT
    fi

    local pairing_rc=0
    set +e
    validate_pairing_and_input_sizes
    pairing_rc=$?
    set -e
    if [[ "$pairing_rc" -ne 0 ]]; then
        exit "$pairing_rc"
    fi

    log "✓ Found ${#FASTQ_FILES[@]} FASTQ file(s)"
    local file
    for file in "${FASTQ_FILES[@]}"; do
        info "  $(basename "$file")"
    done

    check_disk_space

    log "✓ Prerequisites check passed"
}

# Setup directories
setup_directories() {
    safe_mkdir "$OUTPUT_DIR"
    safe_mkdir "$LOG_DIR"
}

# Get file information
get_file_info() {
    log "Getting file information..."

    local file filename filesize
    for file in "${FASTQ_FILES[@]}"; do
        filename=$(basename "$file")
        filesize=$(du -h "$file" | cut -f1)
        log "File: $filename, Size: $filesize"
    done
}

# Run FastQC analysis
run_fastqc() {
    local input_file="$1"
    local filename basename_clean fastqc_log

    filename=$(basename "$input_file")
    basename_clean="${filename%.fq.gz}"
    basename_clean="${basename_clean%.fastq.gz}"
    fastqc_log="$LOG_DIR/fastqc_${basename_clean}.log"

    log "Starting FastQC analysis for $filename..."

    # Check if output already exists
    if [[ -f "$OUTPUT_DIR/${basename_clean}_fastqc.html" ]] && [[ "$FORCE_OVERWRITE" != "true" ]]; then
        warning "FastQC output already exists for $filename. Use --force to overwrite."
        return 0
    fi

    if [[ "$DRY_RUN" == "true" ]]; then
        echo "Would run: fastqc --threads $THREADS --quiet --extract --outdir $OUTPUT_DIR $input_file"
        return 0
    fi

    if [[ ! -f "$input_file" ]]; then
        error "Input file not found: $input_file"
        return $EX_INPUT
    fi

    info "Running FastQC with $THREADS thread(s) on $filename"

    if fastqc --threads "$THREADS" --quiet --extract --outdir "$OUTPUT_DIR" "$input_file" >> "$fastqc_log" 2>&1; then
        local html_report="$OUTPUT_DIR/${basename_clean}_fastqc.html"
        local data_report="$OUTPUT_DIR/${basename_clean}_fastqc/fastqc_data.txt"

        if [[ ! -f "$html_report" || ! -f "$data_report" ]]; then
            error "FastQC reported success but expected outputs are missing for $filename"
            error "Expected: $html_report and $data_report"
            return $EX_OUTPUT_MISSING
        fi

        log "✓ FastQC analysis completed for $filename"
        return 0
    else
        error "FastQC analysis failed for $filename"
        error "Check log file: $fastqc_log"
        return $EX_FASTQC_FAIL
    fi
}

# Analyze FastQC results
analyze_results() {
    local fastqc_dir="$1"
    local filename data_file

    filename=$(basename "$fastqc_dir" _fastqc)
    data_file="$fastqc_dir/fastqc_data.txt"

    log "Analyzing results for $filename..."

    if [[ ! -f "$data_file" ]]; then
        error "FastQC data file not found: $data_file"
        return $EX_OUTPUT_MISSING
    fi

    local total_sequences poor_quality sequence_length gc_content warnings failures
    total_sequences=$(awk -F'\t' '$1=="Total Sequences"{print $2}' "$data_file" | head -1)
    poor_quality=$(awk -F'\t' '$1=="Sequences flagged as poor quality"{print $2}' "$data_file" | head -1)
    sequence_length=$(awk -F'\t' '$1=="Sequence length"{print $2}' "$data_file" | head -1)
    gc_content=$(awk -F'\t' '$1=="%GC"{print $2}' "$data_file" | head -1)

    # Use grep -c with proper tab matching; sanitize to integer to avoid arithmetic errors
    warnings=$(grep -c $'^>>.*\tWARN$' "$data_file" 2>/dev/null || echo "0")
    failures=$(grep -c $'^>>.*\tFAIL$' "$data_file" 2>/dev/null || echo "0")
    # Strip any whitespace/newlines and ensure numeric
    warnings="${warnings//[^0-9]/}"
    failures="${failures//[^0-9]/}"
    warnings="${warnings:-0}"
    failures="${failures:-0}"

    total_sequences=${total_sequences:-unknown}
    poor_quality=${poor_quality:-unknown}
    sequence_length=${sequence_length:-unknown}
    gc_content=${gc_content:-unknown}

    log "Results for $filename:"
    log "  Total sequences: $total_sequences"
    log "  Poor quality sequences: $poor_quality"
    log "  Sequence length: $sequence_length"
    log "  GC content: ${gc_content}%"

    if [[ "$failures" -gt 0 ]]; then
        warning "  $failures quality check(s) FAILED"
    fi

    if [[ "$warnings" -gt 0 ]]; then
        warning "  $warnings quality check(s) generated WARNINGS"
    fi

    if [[ "$failures" -eq 0 && "$warnings" -eq 0 ]]; then
        log "  ✓ All quality checks PASSED"
    fi

    # Guardrail: extremely poor quality input should be explicit but not fatal to QC stage.
    if [[ "$failures" -ge 6 ]]; then
        CRITICAL_QUALITY_FILES+=("$filename ($failures FAIL modules)")
        warning "  Extreme low-quality profile detected."
        warning "  Action: consider less aggressive downstream thresholds or resequencing review."
    fi

    return 0
}

# Generate summary report
generate_summary() {
    local summary_file="$OUTPUT_DIR/quality_control_summary.txt"

    log "Generating summary report..."

    {
        echo "# Quality Control Summary"
        echo "# Analysis Date: $(date)"
        echo "# Script: $SCRIPT_NAME v$SCRIPT_VERSION"
        echo ""
        echo "## Input Files"
        local file filename filesize
        for file in "${FASTQ_FILES[@]}"; do
            filename=$(basename "$file")
            filesize=$(du -h "$file" | cut -f1)
            echo "- $filename: $filesize"
        done
        echo ""
        echo "## Analysis Parameters"
        echo "- FastQC version: $(fastqc --version 2>/dev/null | head -1 || echo unknown)"
        echo "- Threads: $THREADS"
        echo "- Output directory: $OUTPUT_DIR"
        echo ""
        echo "## Per-file Metrics"

        local fastqc_dir metric_file out_name total_sequences poor_quality sequence_length gc_content warnings failures
        for fastqc_dir in "$OUTPUT_DIR"/*_fastqc; do
            [[ -d "$fastqc_dir" ]] || continue

            out_name=$(basename "$fastqc_dir" _fastqc)
            metric_file="$fastqc_dir/fastqc_data.txt"
            echo "### $out_name"

            if [[ -f "$metric_file" ]]; then
                total_sequences=$(awk -F'\t' '$1=="Total Sequences"{print $2}' "$metric_file" | head -1)
                poor_quality=$(awk -F'\t' '$1=="Sequences flagged as poor quality"{print $2}' "$metric_file" | head -1)
                sequence_length=$(awk -F'\t' '$1=="Sequence length"{print $2}' "$metric_file" | head -1)
                gc_content=$(awk -F'\t' '$1=="%GC"{print $2}' "$metric_file" | head -1)
                warnings=$(grep -c $'^>>.*\tWARN$' "$metric_file" 2>/dev/null || echo "0")
                failures=$(grep -c $'^>>.*\tFAIL$' "$metric_file" 2>/dev/null || echo "0")
                # Sanitize to integers
                warnings="${warnings//[^0-9]/}"
                failures="${failures//[^0-9]/}"

                echo "- Total sequences: ${total_sequences:-unknown}"
                echo "- Poor quality sequences: ${poor_quality:-unknown}"
                echo "- Sequence length: ${sequence_length:-unknown}"
                echo "- GC content: ${gc_content:-unknown}%"
                echo "- WARN modules: ${warnings:-0}"
                echo "- FAIL modules: ${failures:-0}"
            else
                echo "- Missing metrics file: $metric_file"
            fi
            echo ""
        done

        if [[ ${#CRITICAL_QUALITY_FILES[@]} -gt 0 ]]; then
            echo "## Critical Quality Alerts"
            local q
            for q in "${CRITICAL_QUALITY_FILES[@]}"; do
                echo "- $q"
            done
            echo ""
        fi

        echo "## Recommendations"
        echo "- Review FastQC HTML reports for adapter content and per-base quality."
        echo "- If many FAIL modules are present, inspect trimming aggressiveness before alignment."
        echo "- Confirm pairing integrity in input delivery pipelines to avoid downstream crashes."
        echo ""

        echo "## Generated HTML Reports"
        local html_found=false
        local html
        for html in "$OUTPUT_DIR"/*_fastqc.html; do
            if [[ -f "$html" ]]; then
                echo "- $html"
                html_found=true
            fi
        done

        if [[ "$html_found" != "true" ]]; then
            echo "- No HTML reports found"
        fi
    } > "$summary_file"

    if [[ ! -s "$summary_file" ]]; then
        error "Failed to generate summary report: $summary_file"
        return $EX_OUTPUT_MISSING
    fi

    log "✓ Summary report generated: $summary_file"
    return 0
}

# Open results for review
open_results() {
    log "Opening results for review..."

    local html_file
    for html_file in "$OUTPUT_DIR"/*.html; do
        [[ -f "$html_file" ]] || continue
        log "Opening: $html_file"
        open "$html_file" 2>/dev/null || {
            warning "Could not open $html_file automatically"
            log "Please manually open: $html_file"
        }
    done
}

# Cleanup function
cleanup() {
    info "Cleanup completed"
}

# Main execution
main() {
    parse_arguments "$@"
    validate_parameters
    setup_directories

    log "Starting Quality Control Analysis for WGS Data"
    log "============================================="

    if [[ "$DRY_RUN" == "true" ]]; then
        log "DRY RUN MODE - No files will be processed"
    fi

    check_prerequisites
    get_file_info

    if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
        error "No FASTQ files found to process"
        exit $EX_INPUT
    fi

    log "Processing ${#FASTQ_FILES[@]} FASTQ file(s)..."

    local failed_files=0
    local first_failure_code=0
    local i current file rc
    for i in "${!FASTQ_FILES[@]}"; do
        file="${FASTQ_FILES[$i]}"
        current=$((i + 1))
        show_progress "$current" "${#FASTQ_FILES[@]}" "Processing $(basename "$file")"

        set +e
        run_fastqc "$file"
        rc=$?
        set -e
        if [[ "$rc" -ne 0 ]]; then
            failed_files=$((failed_files + 1))
            if [[ "$first_failure_code" -eq 0 ]]; then
                first_failure_code=$rc
            fi
            warning "Failed to process: $file (exit code $rc)"
        fi
    done

    echo "" # newline after progress indicator

    if [[ "$failed_files" -gt 0 ]]; then
        warning "$failed_files file(s) failed to process"
    fi

    if [[ "$DRY_RUN" == "true" ]]; then
        log "Dry run completed. No analysis performed."
        exit 0
    fi

    local analysis_count=0
    local analysis_failures=0
    local fastqc_dir analysis_rc

    log "Analyzing FastQC results..."
    for fastqc_dir in "$OUTPUT_DIR"/*_fastqc; do
        [[ -d "$fastqc_dir" ]] || continue
        analysis_count=$((analysis_count + 1))

        set +e
        analyze_results "$fastqc_dir"
        analysis_rc=$?
        set -e

        if [[ "$analysis_rc" -ne 0 ]]; then
            analysis_failures=$((analysis_failures + 1))
            if [[ "$first_failure_code" -eq 0 ]]; then
                first_failure_code=$analysis_rc
            fi
        fi
    done

    if [[ "$analysis_count" -eq 0 ]]; then
        error "No FastQC output directories were generated."
        exit $EX_OUTPUT_MISSING
    fi

    if [[ "$analysis_failures" -gt 0 ]]; then
        warning "$analysis_failures FastQC result directory(ies) could not be analyzed"
    fi

    set +e
    generate_summary
    rc=$?
    set -e
    if [[ "$rc" -ne 0 ]]; then
        exit "$rc"
    fi

    if [[ -t 1 ]]; then
        open_results
    else
        log "Batch mode detected - skipping automatic file opening"
    fi

    cleanup

    if [[ "$failed_files" -eq 0 && "$analysis_failures" -eq 0 ]]; then
        log "✓ Quality Control Analysis Completed Successfully!"
        log "=============================================="
        log "Summary file: $OUTPUT_DIR/quality_control_summary.txt"
        exit 0
    fi

    warning "Quality Control Analysis completed with errors"
    log "=============================================="
    log "Summary file: $OUTPUT_DIR/quality_control_summary.txt"

    if [[ "$first_failure_code" -eq 0 ]]; then
        first_failure_code=$EX_FASTQC_FAIL
    fi

    exit "$first_failure_code"
}

trap 'error "Script interrupted by user"; cleanup; exit 130' INT TERM

main "$@"
