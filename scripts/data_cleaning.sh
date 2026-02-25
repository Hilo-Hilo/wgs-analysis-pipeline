#!/bin/bash

# WGS Data Cleaning Pipeline - 16GB RAM Optimized
# This script performs quality-based trimming and filtering of WGS FASTQ files
# Usage: ./scripts/data_cleaning.sh [OPTIONS]

set -e
set -o pipefail

# Script version and info
SCRIPT_VERSION="2.1-16GB"
SCRIPT_NAME="data_cleaning.sh"

# Exit codes
EX_USAGE=2
EX_PREREQ=10
EX_INPUT=11
EX_PAIRING=12
EX_TINY_FILE=13
EX_FASTP_FAIL=14
EX_OUTPUT_MISSING=15
EX_OVERTRIM=16

# Script paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
CONFIG_FILE="$ROOT_DIR/config/default.conf"

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Missing config file: $CONFIG_FILE" >&2
    exit $EX_INPUT
fi

# shellcheck disable=SC1090
source "$CONFIG_FILE"

# Track CLI overrides explicitly
CLI_SAMPLE_ID=""

# Edge-case guardrails
MIN_FASTQ_BYTES=4096            # <4KB likely invalid for WGS FASTQ
WARN_SMALL_FASTQ_BYTES=1000000  # warn if <1MB
MIN_KEEP_RATIO=0.05             # fail if <5% reads survive trimming
WARN_KEEP_RATIO=0.20            # warn if <20% reads survive trimming

# Runtime metrics
FASTP_BEFORE_READS=""
FASTP_AFTER_READS=""
FASTP_Q30_BEFORE=""
FASTP_Q30_AFTER=""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

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

EDGE-CASE GUARDRAILS:
    - Detects broken FASTQ pairing and duplicate mates
    - Detects tiny/corrupt FASTQ inputs before fastp
    - Validates required outputs are created
    - Flags over-aggressive trimming when <5% reads remain

EXIT CODES:
    0   Success (warnings possible)
    2   Invalid command-line usage
    10  Missing prerequisite/tool
    11  Invalid input files or directories
    12  Broken FASTQ pairing detected
    13  Tiny/corrupt FASTQ detected
    14  fastp execution failure
    15  Missing expected output artifacts
    16  Over-aggressive trimming outcome
EOF
}

safe_mkdir() {
    local dir="$1"
    [[ -n "$dir" ]] && mkdir -p "$dir" 2>/dev/null || true
}

# Flag to track whether logging to file is ready (LOG_DIR resolved to absolute path)
_LOG_TO_FILE_READY=false

log_to_file() {
    local line="$1"
    # Only write to file if logging is ready (LOG_DIR is absolute and set up)
    # This prevents early writes to relative paths before set_defaults runs
    if [[ "$_LOG_TO_FILE_READY" == "true" && -n "$LOG_DIR" ]]; then
        safe_mkdir "$LOG_DIR"
        echo -e "$line" >> "$LOG_DIR/data_cleaning.log" 2>/dev/null || true
    fi
}

enable_file_logging() {
    _LOG_TO_FILE_READY=true
}

# Logging functions
log() {
    local msg="$1"
    local timestamp="[$(date '+%Y-%m-%d %H:%M:%S')]"
    echo -e "${GREEN}${timestamp}${NC} $msg"
    log_to_file "${timestamp} $msg"
}

error() {
    local msg="$1"
    local timestamp="[$(date '+%Y-%m-%d %H:%M:%S')]"
    echo -e "${RED}${timestamp} ERROR:${NC} $msg" >&2
    log_to_file "${timestamp} ERROR: $msg"
}

warning() {
    local msg="$1"
    local timestamp="[$(date '+%Y-%m-%d %H:%M:%S')]"
    echo -e "${YELLOW}${timestamp} WARNING:${NC} $msg"
    log_to_file "${timestamp} WARNING: $msg"
}

info() {
    local msg="$1"
    local timestamp="[$(date '+%Y-%m-%d %H:%M:%S')]"
    if [[ "$VERBOSE" == "true" ]]; then
        echo -e "${BLUE}${timestamp} INFO:${NC} $msg"
    fi
    log_to_file "${timestamp} INFO: $msg"
}

require_value() {
    local flag="$1"
    local value="${2:-}"
    if [[ -z "$value" || "$value" == -* ]]; then
        error "Option $flag requires a value"
        exit $EX_USAGE
    fi
}

get_file_size_bytes() {
    local file="$1"
    stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null || echo 0
}

is_float_lt() {
    local a="$1"
    local b="$2"
    awk -v a="$a" -v b="$b" 'BEGIN { exit !(a < b) }'
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
                echo "$SCRIPT_NAME version $SCRIPT_VERSION"
                exit 0
                ;;
            -i|--input-dir)
                require_value "$1" "${2:-}"
                INPUT_DIR="$2"
                shift 2
                ;;
            -o|--output-dir)
                require_value "$1" "${2:-}"
                OUTPUT_DIR="$2"
                shift 2
                ;;
            -s|--sample-id)
                require_value "$1" "${2:-}"
                CLI_SAMPLE_ID="$2"
                SAMPLE_ID="$2"
                shift 2
                ;;
            -t|--threads)
                require_value "$1" "${2:-}"
                THREADS="$2"
                shift 2
                ;;
            --min-length)
                require_value "$1" "${2:-}"
                MIN_LENGTH="$2"
                shift 2
                ;;
            --quality)
                require_value "$1" "${2:-}"
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
                error "Unknown option: $1"
                echo "Use --help for usage information"
                exit $EX_USAGE
                ;;
        esac
    done
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

    INPUT_DIR=$(realpath "$INPUT_DIR" 2>/dev/null || echo "$INPUT_DIR")
    OUTPUT_DIR=$(realpath "$OUTPUT_DIR" 2>/dev/null || echo "$OUTPUT_DIR")
    LOG_DIR=$(realpath "$LOG_DIR" 2>/dev/null || echo "$LOG_DIR")
}

validate_parameters() {
    if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -le 0 ]]; then
        error "Invalid --threads value: $THREADS"
        exit $EX_USAGE
    fi

    if ! [[ "$MIN_LENGTH" =~ ^[0-9]+$ ]] || [[ "$MIN_LENGTH" -le 0 ]]; then
        error "Invalid --min-length value: $MIN_LENGTH"
        exit $EX_USAGE
    fi

    if ! [[ "$QUALITY_THRESHOLD" =~ ^[0-9]+$ ]] || [[ "$QUALITY_THRESHOLD" -lt 0 || "$QUALITY_THRESHOLD" -gt 40 ]]; then
        error "Invalid --quality value: $QUALITY_THRESHOLD (expected 0-40)"
        exit $EX_USAGE
    fi
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

# Auto-detect input FASTQ files with strict pairing checks
detect_input_files() {
    info "Detecting FASTQ files in: $INPUT_DIR"

    if [[ ! -d "$INPUT_DIR" ]]; then
        error "Input directory not found: $INPUT_DIR"
        exit $EX_INPUT
    fi

    local all_fastq=()
    while IFS= read -r -d '' file; do
        all_fastq+=("$file")
    done < <(find "$INPUT_DIR" -maxdepth 1 -type f \( -name "*.fq.gz" -o -name "*.fastq.gz" \) -print0 2>/dev/null)

    if [[ ${#all_fastq[@]} -eq 0 ]]; then
        error "No FASTQ files found in: $INPUT_DIR"
        echo "Expected files with extensions: .fq.gz or .fastq.gz"
        exit $EX_INPUT
    fi

    PAIR_KEYS=()
    PAIR_R1_FILES=()
    PAIR_R2_FILES=()

    local file base pair_info key mate idx i
    local recognized=0

    for file in "${all_fastq[@]}"; do
        base=$(basename "$file")

        if ! pair_info=$(pair_key_and_mate "$base"); then
            warning "Ignoring FASTQ with unrecognized paired-end naming: $base"
            warning "Expected *_R1/*_R2 or *_1/*_2 naming pattern"
            continue
        fi

        recognized=$((recognized + 1))
        key="${pair_info%%|*}"
        mate="${pair_info##*|}"

        idx=-1
        for i in "${!PAIR_KEYS[@]}"; do
            if [[ "${PAIR_KEYS[$i]}" == "$key" ]]; then
                idx=$i
                break
            fi
        done

        if [[ $idx -lt 0 ]]; then
            PAIR_KEYS+=("$key")
            PAIR_R1_FILES+=("")
            PAIR_R2_FILES+=("")
            idx=$((${#PAIR_KEYS[@]} - 1))
        fi

        if [[ "$mate" == "R1" ]]; then
            if [[ -n "${PAIR_R1_FILES[$idx]}" ]]; then
                error "Duplicate R1 detected for sample '$key':"
                error "  ${PAIR_R1_FILES[$idx]}"
                error "  $file"
                exit $EX_PAIRING
            fi
            PAIR_R1_FILES[$idx]="$file"
        else
            if [[ -n "${PAIR_R2_FILES[$idx]}" ]]; then
                error "Duplicate R2 detected for sample '$key':"
                error "  ${PAIR_R2_FILES[$idx]}"
                error "  $file"
                exit $EX_PAIRING
            fi
            PAIR_R2_FILES[$idx]="$file"
        fi
    done

    if [[ "$recognized" -eq 0 ]]; then
        error "No properly named paired-end FASTQ files found in: $INPUT_DIR"
        error "Action: rename inputs to *_R1/*_R2 or *_1/*_2 conventions."
        exit $EX_PAIRING
    fi

    local complete_pairs=0
    for i in "${!PAIR_KEYS[@]}"; do
        if [[ -z "${PAIR_R1_FILES[$i]}" || -z "${PAIR_R2_FILES[$i]}" ]]; then
            error "Broken FASTQ pairing for sample '${PAIR_KEYS[$i]}':"
            error "  R1: ${PAIR_R1_FILES[$i]:-(missing)}"
            error "  R2: ${PAIR_R2_FILES[$i]:-(missing)}"
            error "Action: ensure one and only one mate file exists for each pair."
            exit $EX_PAIRING
        fi
        complete_pairs=$((complete_pairs + 1))
    done

    local selected_index=0
    if [[ -n "$CLI_SAMPLE_ID" ]]; then
        local found=false
        for i in "${!PAIR_KEYS[@]}"; do
            if [[ "${PAIR_KEYS[$i]}" == "$CLI_SAMPLE_ID" ]]; then
                selected_index=$i
                found=true
                break
            fi
        done

        if [[ "$found" != "true" ]]; then
            warning "Requested --sample-id '$CLI_SAMPLE_ID' did not match detected pair names."
            warning "Proceeding with first detected pair: ${PAIR_KEYS[0]}"
        fi
    elif [[ "$complete_pairs" -gt 1 ]]; then
        warning "Multiple FASTQ pairs found. Processing first pair only."
        for i in "${!PAIR_KEYS[@]}"; do
            warning "  Pair $((i + 1)): ${PAIR_KEYS[$i]}"
        done
        warning "Use --sample-id to disambiguate." 
    fi

    RAW_R1="${PAIR_R1_FILES[$selected_index]}"
    RAW_R2="${PAIR_R2_FILES[$selected_index]}"

    info "Selected input files:"
    info "  R1: $RAW_R1"
    info "  R2: $RAW_R2"
}

validate_fastq_input_file() {
    local file="$1"
    local label="$2"
    local size

    if [[ ! -f "$file" ]]; then
        error "$label file not found: $file"
        return $EX_INPUT
    fi

    size=$(get_file_size_bytes "$file")
    if [[ "$size" -lt "$MIN_FASTQ_BYTES" ]]; then
        error "$label file appears too small (${size} bytes): $file"
        error "Action: verify transfer/decompression; file may be truncated."
        return $EX_TINY_FILE
    elif [[ "$size" -lt "$WARN_SMALL_FASTQ_BYTES" ]]; then
        warning "$label file is very small (${size} bytes): $file"
        warning "This may be a tiny test dataset or truncated upload."
    fi

    if ! gzip -t "$file" 2>/dev/null; then
        error "$label file is not a valid gzip stream: $file"
        return $EX_INPUT
    fi

    return 0
}

# Check prerequisites
check_prerequisites() {
    log "Checking prerequisites for data cleaning..."

    if ! command -v fastp >/dev/null 2>&1; then
        error "fastp not found. Install with: conda install -c bioconda fastp"
        exit $EX_PREREQ
    fi

    local rc=0
    set +e
    validate_fastq_input_file "$RAW_R1" "Raw R1"
    rc=$?
    set -e
    if [[ "$rc" -ne 0 ]]; then
        exit "$rc"
    fi

    set +e
    validate_fastq_input_file "$RAW_R2" "Raw R2"
    rc=$?
    set -e
    if [[ "$rc" -ne 0 ]]; then
        exit "$rc"
    fi

    # Cross-platform available memory check
    local available_mem_gb=""
    if [[ -r /proc/meminfo ]]; then
        local available_mem_kb
        available_mem_kb=$(awk '/MemAvailable/ {print $2}' /proc/meminfo 2>/dev/null || echo "0")
        available_mem_gb=$((available_mem_kb / 1024 / 1024))
    elif command -v vm_stat >/dev/null 2>&1; then
        local page_size free_pages inactive_pages speculative_pages
        page_size=$(vm_stat 2>/dev/null | awk '/page size of/ {gsub("\\.", "", $8); print $8}' | head -1)
        free_pages=$(vm_stat 2>/dev/null | awk '/Pages free/ {gsub("\\.", "", $3); print $3}' | head -1)
        inactive_pages=$(vm_stat 2>/dev/null | awk '/Pages inactive/ {gsub("\\.", "", $3); print $3}' | head -1)
        speculative_pages=$(vm_stat 2>/dev/null | awk '/Pages speculative/ {gsub("\\.", "", $3); print $3}' | head -1)
        page_size=${page_size:-4096}
        free_pages=${free_pages:-0}
        inactive_pages=${inactive_pages:-0}
        speculative_pages=${speculative_pages:-0}
        local free_bytes=$(( (free_pages + inactive_pages + speculative_pages) * page_size ))
        available_mem_gb=$((free_bytes / 1024 / 1024 / 1024))
    fi

    if [[ -n "$available_mem_gb" && "$available_mem_gb" -lt 6 ]]; then
        warning "Low available memory: ${available_mem_gb}GB. Data cleaning may be slow."
    fi

    log "Prerequisites check passed"
}

# Setup directories
setup_directories() {
    # Now that LOG_DIR is resolved to absolute path, enable file logging
    safe_mkdir "$LOG_DIR"
    enable_file_logging
    
    log "Setting up directories..."
    safe_mkdir "$OUTPUT_DIR"
    safe_mkdir "$OUTPUT_DIR/reports"
    info "Output directory: $OUTPUT_DIR"
}

validate_fastp_outputs() {
    local cleaned_r1="$1"
    local cleaned_r2="$2"
    local report_html="$3"
    local report_json="$4"

    local missing=0

    if [[ ! -s "$cleaned_r1" ]]; then
        error "Missing or empty cleaned R1 output: $cleaned_r1"
        missing=1
    fi

    if [[ ! -s "$cleaned_r2" ]]; then
        error "Missing or empty cleaned R2 output: $cleaned_r2"
        missing=1
    fi

    if [[ ! -s "$report_html" ]]; then
        error "Missing or empty fastp HTML report: $report_html"
        missing=1
    fi

    if [[ ! -s "$report_json" ]]; then
        error "Missing or empty fastp JSON report: $report_json"
        missing=1
    fi

    if [[ "$missing" -ne 0 ]]; then
        return $EX_OUTPUT_MISSING
    fi

    if ! gzip -t "$cleaned_r1" 2>/dev/null; then
        error "Cleaned R1 is not valid gzip: $cleaned_r1"
        return $EX_OUTPUT_MISSING
    fi

    if ! gzip -t "$cleaned_r2" 2>/dev/null; then
        error "Cleaned R2 is not valid gzip: $cleaned_r2"
        return $EX_OUTPUT_MISSING
    fi

    return 0
}

extract_fastp_metrics() {
    local json_file="$1"

    if ! command -v python3 >/dev/null 2>&1; then
        return 1
    fi

    python3 - "$json_file" << 'PY'
import json
import sys

path = sys.argv[1]
try:
    with open(path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    before = int(data.get('summary', {}).get('before_filtering', {}).get('total_reads', -1))
    after = int(data.get('summary', {}).get('after_filtering', {}).get('total_reads', -1))

    q30_before = data.get('summary', {}).get('before_filtering', {}).get('q30_rate', -1)
    q30_after = data.get('summary', {}).get('after_filtering', {}).get('q30_rate', -1)

    if q30_before is None:
        q30_before = -1
    if q30_after is None:
        q30_after = -1

    print(f"{before}\t{after}\t{q30_before}\t{q30_after}")
except Exception:
    sys.exit(1)
PY
}

evaluate_cleaning_quality() {
    local report_json="$1"

    local metrics
    if ! metrics=$(extract_fastp_metrics "$report_json" 2>/dev/null); then
        warning "Unable to parse fastp JSON metrics from $report_json"
        warning "Skipping quantitative over-trimming checks."
        return 0
    fi

    FASTP_BEFORE_READS=$(echo "$metrics" | awk -F'\t' '{print $1}')
    FASTP_AFTER_READS=$(echo "$metrics" | awk -F'\t' '{print $2}')
    FASTP_Q30_BEFORE=$(echo "$metrics" | awk -F'\t' '{print $3}')
    FASTP_Q30_AFTER=$(echo "$metrics" | awk -F'\t' '{print $4}')

    if [[ "$FASTP_BEFORE_READS" -le 0 || "$FASTP_AFTER_READS" -lt 0 ]]; then
        warning "fastp report contains unexpected read-count metrics: before=$FASTP_BEFORE_READS after=$FASTP_AFTER_READS"
        return 0
    fi

    if [[ "$FASTP_AFTER_READS" -eq 0 ]]; then
        error "All reads were removed during cleaning (after_filtering.total_reads = 0)."
        error "Action: relax --quality or --min-length thresholds and rerun."
        return $EX_OVERTRIM
    fi

    local keep_ratio
    keep_ratio=$(awk -v a="$FASTP_AFTER_READS" -v b="$FASTP_BEFORE_READS" 'BEGIN { if (b==0) print 0; else printf "%.4f", (a/b) }')
    log "Read retention ratio: $keep_ratio (after=$FASTP_AFTER_READS / before=$FASTP_BEFORE_READS)"

    if is_float_lt "$keep_ratio" "$MIN_KEEP_RATIO"; then
        error "Over-aggressive trimming detected: only $keep_ratio of reads retained (<$MIN_KEEP_RATIO)."
        error "Action: lower --quality and/or --min-length; inspect fastp report before proceeding."
        return $EX_OVERTRIM
    fi

    if is_float_lt "$keep_ratio" "$WARN_KEEP_RATIO"; then
        warning "Low read retention ratio: $keep_ratio (<$WARN_KEEP_RATIO)."
        warning "Downstream coverage may be insufficient; review report and consider threshold tuning."
    fi

    if is_float_lt "$FASTP_Q30_BEFORE" "0.30"; then
        warning "Input appears extremely low-quality (Q30 before filtering=${FASTP_Q30_BEFORE})."
        warning "Action: expect significant read loss and monitor downstream coverage."
    fi

    if is_float_lt "$FASTP_Q30_AFTER" "$FASTP_Q30_BEFORE"; then
        warning "Post-cleaning Q30 did not improve (before=${FASTP_Q30_BEFORE}, after=${FASTP_Q30_AFTER})."
        warning "Action: inspect adapter/quality settings and raw data integrity."
    fi

    return 0
}

# Run fastp cleaning (16GB optimized)
run_fastp_cleaning() {
    log "Starting fastp data cleaning (16GB optimized)..."

    CLEANED_R1="$OUTPUT_DIR/${SAMPLE_ID}_clean_R1.fq.gz"
    CLEANED_R2="$OUTPUT_DIR/${SAMPLE_ID}_clean_R2.fq.gz"
    REPORT_DIR="$OUTPUT_DIR/reports"
    REPORT_HTML="${REPORT_DIR}/${SAMPLE_ID}_cleaning_report.html"
    REPORT_JSON="${REPORT_DIR}/${SAMPLE_ID}_cleaning_report.json"

    if [[ -f "$CLEANED_R1" && -f "$CLEANED_R2" && "$FORCE_OVERWRITE" != "true" ]]; then
        warning "Cleaned reads already exist:"
        echo "  R1: $CLEANED_R1"
        echo "  R2: $CLEANED_R2"
        echo "Use --force to overwrite"
        return 0
    fi

    if [[ "$DRY_RUN" == "true" ]]; then
        echo "Would run fastp cleaning with:"
        echo "  Input R1: $RAW_R1"
        echo "  Input R2: $RAW_R2"
        echo "  Output R1: $CLEANED_R1"
        echo "  Output R2: $CLEANED_R2"
        echo "  Report HTML: $REPORT_HTML"
        echo "  Report JSON: $REPORT_JSON"
        echo "  Threads: $THREADS"
        echo "  Min length: $MIN_LENGTH"
        echo "  Quality: $QUALITY_THRESHOLD"
        return 0
    fi

    info "Running fastp with 16GB optimized settings..."
    info "Threads: $THREADS"
    info "Quality threshold: $QUALITY_THRESHOLD, Min length: $MIN_LENGTH"

    local fastp_cmd=(
        fastp
        -i "$RAW_R1" -I "$RAW_R2"
        -o "$CLEANED_R1" -O "$CLEANED_R2"
        --html "$REPORT_HTML" --json "$REPORT_JSON"
        --thread "$THREADS"
        --qualified_quality_phred "$QUALITY_THRESHOLD"
        --unqualified_percent_limit 20
        --length_required "$MIN_LENGTH"
        --detect_adapter_for_pe
        --correction
        --cut_tail --cut_tail_window_size 4
        --cut_tail_mean_quality "$QUALITY_THRESHOLD"
        --overrepresentation_analysis
    )

    if "${fastp_cmd[@]}" >> "$LOG_DIR/data_cleaning.log" 2>&1; then
        log "Data cleaning completed successfully"

        local rc=0
        set +e
        validate_fastp_outputs "$CLEANED_R1" "$CLEANED_R2" "$REPORT_HTML" "$REPORT_JSON"
        rc=$?
        set -e
        if [[ "$rc" -ne 0 ]]; then
            return "$rc"
        fi

        set +e
        evaluate_cleaning_quality "$REPORT_JSON"
        rc=$?
        set -e
        if [[ "$rc" -ne 0 ]]; then
            return "$rc"
        fi

        local r1_size_raw r2_size_raw r1_size_clean r2_size_clean
        r1_size_raw=$(du -h "$RAW_R1" | cut -f1)
        r2_size_raw=$(du -h "$RAW_R2" | cut -f1)
        r1_size_clean=$(du -h "$CLEANED_R1" | cut -f1)
        r2_size_clean=$(du -h "$CLEANED_R2" | cut -f1)

        log "File size comparison:"
        log "  Raw R1: $r1_size_raw -> Cleaned R1: $r1_size_clean"
        log "  Raw R2: $r2_size_raw -> Cleaned R2: $r2_size_clean"

        return 0
    else
        error "Data cleaning failed"
        error "Check log file: $LOG_DIR/data_cleaning.log"
        return $EX_FASTP_FAIL
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

    {
        echo "# Data Cleaning Summary - $SAMPLE_ID"
        echo "# Generated: $(date)"
        echo "# Tool: fastp with 16GB RAM optimizations"
        echo ""
        echo "## Input Files"
        echo "- Forward reads: $RAW_R1"
        echo "- Reverse reads: $RAW_R2"
        echo ""
        echo "## Cleaning Parameters"
        echo "- Tool: fastp $(fastp --version 2>&1 | head -1 || echo unknown)"
        echo "- Threads: $THREADS"
        echo "- Quality threshold: Q$QUALITY_THRESHOLD"
        echo "- Minimum length: ${MIN_LENGTH}bp"
        echo "- Adapter detection: Automatic paired-end"
        echo "- Tail cutting: Window size 4, mean quality $QUALITY_THRESHOLD"
        echo "- Correction: Enabled"
        echo "- Overrepresentation analysis: Enabled"
        echo ""

        if [[ -n "$FASTP_BEFORE_READS" && -n "$FASTP_AFTER_READS" ]]; then
            local keep_ratio
            keep_ratio=$(awk -v a="$FASTP_AFTER_READS" -v b="$FASTP_BEFORE_READS" 'BEGIN { if (b==0) print 0; else printf "%.4f", (a/b) }')
            echo "## Cleaning Metrics"
            echo "- Reads before filtering: $FASTP_BEFORE_READS"
            echo "- Reads after filtering: $FASTP_AFTER_READS"
            echo "- Retention ratio: $keep_ratio"
            echo "- Q30 before: ${FASTP_Q30_BEFORE:-unknown}"
            echo "- Q30 after: ${FASTP_Q30_AFTER:-unknown}"
            echo ""
        fi

        echo "## Output Files"
        echo "- Cleaned forward reads: $CLEANED_R1"
        echo "- Cleaned reverse reads: $CLEANED_R2"
        echo "- HTML report: $REPORT_HTML"
        echo "- JSON report: $REPORT_JSON"
        echo ""

        echo "## File Sizes"
        if [[ -f "$RAW_R1" ]]; then
            echo "- Raw R1: $(du -h "$RAW_R1" | cut -f1)"
        fi
        if [[ -f "$RAW_R2" ]]; then
            echo "- Raw R2: $(du -h "$RAW_R2" | cut -f1)"
        fi
        if [[ -f "$CLEANED_R1" ]]; then
            echo "- Cleaned R1: $(du -h "$CLEANED_R1" | cut -f1)"
        fi
        if [[ -f "$CLEANED_R2" ]]; then
            echo "- Cleaned R2: $(du -h "$CLEANED_R2" | cut -f1)"
        fi
        echo ""

        echo "## Next Steps"
        echo "1. Review cleaning report: $REPORT_HTML"
        echo "2. Proceed with alignment: ./scripts/alignment.sh"
        echo "3. Monitor system resources during alignment"
        echo ""
        echo "## Analysis completed: $(date)"
    } > "$summary_file"

    if [[ ! -s "$summary_file" ]]; then
        error "Failed to write summary report: $summary_file"
        return $EX_OUTPUT_MISSING
    fi

    log "Summary report generated: $summary_file"
    return 0
}

# Cleanup temporary files
cleanup() {
    info "Cleanup completed"
}

# Main execution
main() {
    parse_arguments "$@"
    set_defaults
    validate_parameters
    setup_directories

    log "Starting 16GB Optimized Data Cleaning Pipeline"
    log "============================================="

    if [[ "$DRY_RUN" == "true" ]]; then
        log "DRY RUN MODE - No files will be processed"
    fi

    info "Sample ID: $SAMPLE_ID"
    info "Input directory: $INPUT_DIR"
    info "Output directory: $OUTPUT_DIR"
    info "Threads: $THREADS (16GB optimized)"

    detect_input_files
    check_prerequisites

    if [[ "$DRY_RUN" == "true" ]]; then
        log "Dry run completed"
        exit 0
    fi

    local rc=0

    set +e
    run_fastp_cleaning
    rc=$?
    set -e
    if [[ "$rc" -ne 0 ]]; then
        exit "$rc"
    fi

    set +e
    generate_summary
    rc=$?
    set -e
    if [[ "$rc" -ne 0 ]]; then
        exit "$rc"
    fi

    cleanup

    log "16GB Optimized Data Cleaning Pipeline Completed!"
    log "=============================================="
    log "Results:"
    log "  Cleaned reads: $OUTPUT_DIR/${SAMPLE_ID}_clean_R*.fq.gz"
    log "  Cleaning report: $REPORT_HTML"
    log "  Summary: $OUTPUT_DIR/${SAMPLE_ID}_cleaning_summary.txt"
    log ""
    log "Next step: Alignment with ./scripts/alignment.sh"
}

trap 'error "Script interrupted by user"; cleanup; exit 130' INT TERM

main "$@"
