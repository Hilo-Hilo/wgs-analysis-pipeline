#!/usr/bin/env bash

# WGS Pipeline Runner with Progress Monitoring
# Enhanced pipeline execution with real-time progress tracking and resource monitoring

set -e

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source hardware detection (safe on Bash 3+)
source "$SCRIPT_DIR/scripts/detect_hardware.sh"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Default configuration
DEFAULT_THREADS=4
DEFAULT_SAMPLE_ID="WGS_SAMPLE"
DEFAULT_USE_GPU=false
DEFAULT_GPU_ALIGNER="parabricks"
DEFAULT_GPU_COUNT=1
ENABLE_NOTIFICATIONS=false
EMAIL_ADDRESS=""
SLACK_WEBHOOK=""
RESUME_FROM=""
ASSUME_YES=false
SKIP_REQUIREMENTS_CHECK=false
REGISTRY_DB=""

# Auto-detection mode (enabled by default, disabled with explicit overrides)
AUTO_DETECT=true
DETECTED_MODE=""
CLI_THREADS=""
CLI_USE_GPU=""
CLI_GPU_COUNT=""

# Valid CLI options for suggestions/validation
VALID_OPTIONS=(
    -h --help
    -s --sample-id
    -t --threads
    -i --input-dir
    -o --output-dir
    --use-gpu
    --no-gpu
    --gpu-aligner
    --gpu-count
    --no-auto
    --show-hardware
    --resume-from
    --notify-email
    --notify-slack
    --steps
    -y --yes
    --skip-requirements-check
    --dry-run
    --verbose
)

AVAILABLE_STEPS=(
    "quality-control"
    "data-cleaning"
    "alignment"
    "variant-calling"
    "annotation"
)

trim() {
    local value="$1"
    # shellcheck disable=SC2001
    value="$(echo "$value" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
    printf '%s' "$value"
}

suggest_closest() {
    local needle="$1"
    shift || true

    if [[ $# -eq 0 ]]; then
        return 0
    fi

    if command -v python3 >/dev/null 2>&1; then
        python3 - "$needle" "$@" <<'PY'
import difflib
import sys

needle = sys.argv[1]
candidates = sys.argv[2:]
match = difflib.get_close_matches(needle, candidates, n=1, cutoff=0.45)
if match:
    print(match[0])
PY
    fi
}

print_cli_error() {
    local message="$1"
    echo -e "${RED}❌ ${message}${NC}" >&2
    echo "Run '$0 --help' for usage and examples." >&2
}

unknown_option_error() {
    local bad_opt="$1"
    local suggestion
    suggestion="$(suggest_closest "$bad_opt" "${VALID_OPTIONS[@]}")"

    if [[ -n "$suggestion" ]]; then
        print_cli_error "Unknown option: $bad_opt (did you mean '$suggestion'?)"
    else
        print_cli_error "Unknown option: $bad_opt"
    fi
    exit 2
}

missing_value_error() {
    local opt="$1"
    print_cli_error "Option '$opt' requires a value."
    exit 2
}

normalize_step_name() {
    local step
    step="$(trim "$1")"
    case "$step" in
        quality_control) echo "quality-control" ;;
        data_cleaning) echo "data-cleaning" ;;
        variant_calling) echo "variant-calling" ;;
        *) echo "$step" ;;
    esac
}

require_bash4_or_die() {
    local major="${BASH_VERSINFO[0]:-0}"

    if [[ -z "${BASH_VERSION:-}" || "$major" -lt 4 ]]; then
        echo -e "${RED}❌ run_pipeline.sh requires GNU Bash 4.0 or newer.${NC}" >&2
        echo "Detected shell: ${BASH_VERSION:-unknown}" >&2
        echo >&2
        echo "Why this is required:" >&2
        echo "  - Progress monitoring uses associative arrays (Bash 4+ feature)." >&2
        echo >&2
        echo "How to fix on macOS:" >&2
        echo "  1) Install modern bash (e.g. Homebrew): brew install bash" >&2
        echo "  2) Run with that binary, for example:" >&2
        echo "     /opt/homebrew/bin/bash run_pipeline.sh --help" >&2
        exit 1
    fi
}

validate_positive_integer() {
    local name="$1"
    local value="$2"

    if [[ ! "$value" =~ ^[1-9][0-9]*$ ]]; then
        print_cli_error "$name must be a positive integer (got: '$value')."
        exit 2
    fi
}

validate_step_name() {
    local step="$1"
    local context="$2"

    if [[ -z "${PIPELINE_STEPS[$step]:-}" ]]; then
        local suggestion
        suggestion="$(suggest_closest "$step" "${AVAILABLE_STEPS[@]}")"

        if [[ -n "$suggestion" ]]; then
            print_cli_error "Invalid ${context}: '$step' (did you mean '$suggestion'?)"
        else
            print_cli_error "Invalid ${context}: '$step'"
        fi

        echo "Available steps: ${AVAILABLE_STEPS[*]}" >&2
        exit 2
    fi
}

# Help function
show_help() {
    cat << EOF
WGS Pipeline Runner with Progress Monitoring

USAGE:
    $0 [OPTIONS]

CORE OPTIONS:
    -h, --help                  Show this help message
    -s, --sample-id ID          Sample identifier (default: WGS_SAMPLE)
    -t, --threads NUM           Number of threads (overrides auto-detection)
    -i, --input-dir DIR         Input directory with raw FASTQ files
    -o, --output-dir DIR        Output directory (default: results)
    --use-gpu                   Use GPU alignment (auto-enabled if GPUs detected)
    --no-gpu                    Disable GPU even if available
    --gpu-aligner NAME          GPU aligner backend (default: parabricks)
    --gpu-count NUM             Number of GPUs (overrides auto-detection)
    --no-auto                   Disable hardware auto-detection entirely
    --show-hardware             Show detected hardware and exit
    --resume-from STEP          Resume pipeline from specific step
    --notify-email EMAIL        Send email notifications
    --notify-slack WEBHOOK      Send Slack notifications
    --steps STEPS               Run specific steps only (comma-separated)

    -t, --threads NUM           Number of threads (default: 4)
    -i, --input-dir DIR         Input directory with raw FASTQ files (default: data/raw)
    -o, --output-dir DIR        Output directory (default: results)
    --steps STEPS               Run specific steps only (comma-separated)
    --resume-from STEP          Resume pipeline from a specific step

GPU OPTIONS (alignment step):
    --use-gpu                   Enable GPU alignment mode
    --gpu-aligner NAME          GPU aligner backend (default: parabricks)
    --gpu-count NUM             Number of GPUs for alignment step (default: 1)

SAFETY / EXECUTION OPTIONS:
    -y, --yes                   Continue on requirement warnings without prompting
    --skip-requirements-check   Skip scripts/check_requirements.sh

    --skip-requirements-check   Skip running scripts/check_requirements.sh
    --registry-db PATH          Optional SQLite sample registry path
    --dry-run                   Show what would be done without executing
    --verbose                   Enable verbose output

HARDWARE AUTO-DETECTION:
    The pipeline automatically detects available CPU cores, RAM, and GPUs
    to select optimal settings. Manual overrides (--threads, --use-gpu,
    --gpu-count) take precedence over auto-detected values.
    
    Performance modes (auto-selected based on hardware):
    - laptop:      2-4 threads, CPU-only (8-16GB RAM)
    - workstation: 4-8 threads, GPU if available (32GB+ RAM, 8+ cores)
    - server:      16-24 threads, GPU if available (64GB+ RAM, 16+ cores)
    - dgx:         32+ threads, multi-GPU (8+ GPUs or 128GB+ RAM)

NOTIFICATION OPTIONS:
    --notify-email EMAIL        Send email notifications
    --notify-slack WEBHOOK      Send Slack notifications

PIPELINE STEPS:
    quality-control             FastQC analysis of raw reads
    data-cleaning               Adapter trimming and quality filtering
    alignment                   BWA/Parabricks mapping to reference genome
    variant-calling             bcftools variant calling
    annotation                  VEP variant annotation

SHELL COMPATIBILITY:
    This script requires GNU Bash 4+.
    macOS ships Bash 3.2 by default, which is not sufficient.
    Install newer bash (e.g. Homebrew) and run with that binary.

EXAMPLES:
    # Run complete pipeline (auto-detects optimal settings)
    $0 --sample-id MySample --input-dir data/raw

    # Show detected hardware without running
    $0 --show-hardware

    # Override auto-detected threads
    $0 --sample-id MySample --threads 8

    # Force CPU mode even if GPUs available
    $0 --sample-id MySample --no-gpu

    # Force GPU mode with specific GPU count
    $0 --sample-id MySample --use-gpu --gpu-count 2

    # Disable auto-detection entirely (use defaults)
    $0 --sample-id MySample --no-auto

    # Run with notifications
    $0 --sample-id MySample --notify-email user@example.com

    # Resume from alignment step
    $0 --resume-from alignment

    # Run specific steps only

    # Run with GPU alignment
    $0 --sample-id MySample --use-gpu --gpu-aligner parabricks --gpu-count 1

    # Run only QC + cleaning
    $0 --steps quality-control,data-cleaning

    # Resume from alignment
    $0 --resume-from alignment

    # Dry run for configuration validation
    $0 --dry-run --verbose
EOF
}

# Parse command line arguments
parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -h|--help)
                show_help
                exit 0
                ;;

            -s|--sample-id)
                [[ $# -lt 2 ]] && missing_value_error "$1"
                SAMPLE_ID="$2"
                shift 2
                ;;
            --sample-id=*)
                SAMPLE_ID="${1#*=}"
                [[ -z "$SAMPLE_ID" ]] && missing_value_error "--sample-id"
                shift
                ;;

            -t|--threads)
                CLI_THREADS="$2"

                [[ $# -lt 2 ]] && missing_value_error "$1"
                THREADS="$2"
                shift 2
                ;;
            --threads=*)
                THREADS="${1#*=}"
                [[ -z "$THREADS" ]] && missing_value_error "--threads"
                shift
                ;;

            -i|--input-dir)
                [[ $# -lt 2 ]] && missing_value_error "$1"
                INPUT_DIR="$2"
                shift 2
                ;;
            --input-dir=*)
                INPUT_DIR="${1#*=}"
                [[ -z "$INPUT_DIR" ]] && missing_value_error "--input-dir"
                shift
                ;;

            -o|--output-dir)
                [[ $# -lt 2 ]] && missing_value_error "$1"
                OUTPUT_DIR="$2"
                shift 2
                ;;
            --output-dir=*)
                OUTPUT_DIR="${1#*=}"
                [[ -z "$OUTPUT_DIR" ]] && missing_value_error "--output-dir"
                shift
                ;;

            --use-gpu)
                CLI_USE_GPU="true"
                USE_GPU=true
                shift
                ;;
            --no-gpu)
                CLI_USE_GPU="false"
                USE_GPU=false
                shift
                ;;

            --gpu-aligner)
                [[ $# -lt 2 ]] && missing_value_error "$1"
                GPU_ALIGNER="$2"
                shift 2
                ;;
            --gpu-aligner=*)
                GPU_ALIGNER="${1#*=}"
                [[ -z "$GPU_ALIGNER" ]] && missing_value_error "--gpu-aligner"
                shift
                ;;

            --gpu-count)
                [[ $# -lt 2 ]] && missing_value_error "$1"
                CLI_GPU_COUNT="$2"
                GPU_COUNT="$2"
                shift 2
                ;;
            --gpu-count=*)
                GPU_COUNT="${1#*=}"
                [[ -z "$GPU_COUNT" ]] && missing_value_error "--gpu-count"
                CLI_GPU_COUNT="$GPU_COUNT"
                shift
                ;;

            --no-auto)
                AUTO_DETECT=false
                shift
                ;;
            --show-hardware)
                print_summary
                exit 0
                ;;

            --resume-from)
                [[ $# -lt 2 ]] && missing_value_error "$1"
                RESUME_FROM="$2"
                shift 2
                ;;
            --resume-from=*)
                RESUME_FROM="${1#*=}"
                [[ -z "$RESUME_FROM" ]] && missing_value_error "--resume-from"
                shift
                ;;

            --notify-email)
                [[ $# -lt 2 ]] && missing_value_error "$1"
                EMAIL_ADDRESS="$2"
                ENABLE_NOTIFICATIONS=true
                shift 2
                ;;
            --notify-email=*)
                EMAIL_ADDRESS="${1#*=}"
                [[ -z "$EMAIL_ADDRESS" ]] && missing_value_error "--notify-email"
                ENABLE_NOTIFICATIONS=true
                shift
                ;;

            --notify-slack)
                [[ $# -lt 2 ]] && missing_value_error "$1"
                SLACK_WEBHOOK="$2"
                ENABLE_NOTIFICATIONS=true
                shift 2
                ;;
            --notify-slack=*)
                SLACK_WEBHOOK="${1#*=}"
                [[ -z "$SLACK_WEBHOOK" ]] && missing_value_error "--notify-slack"
                ENABLE_NOTIFICATIONS=true
                shift
                ;;

            --steps)
                [[ $# -lt 2 ]] && missing_value_error "$1"
                SPECIFIC_STEPS="$2"
                shift 2
                ;;
            --steps=*)
                SPECIFIC_STEPS="${1#*=}"
                [[ -z "$SPECIFIC_STEPS" ]] && missing_value_error "--steps"
                shift
                ;;

            -y|--yes)
                ASSUME_YES=true
                shift
                ;;

            --skip-requirements-check)
                SKIP_REQUIREMENTS_CHECK=true
                shift
                ;;

            --registry-db)
                REGISTRY_DB="$2"
                shift 2
                ;;
            --dry-run)
                DRY_RUN=true
                shift
                ;;

            --verbose)
                VERBOSE=true
                shift
                ;;

            --)
                # No positional arguments currently supported.
                shift
                if [[ $# -gt 0 ]]; then
                    print_cli_error "Unexpected positional arguments: $*"
                    exit 2
                fi
                ;;

            -*)
                unknown_option_error "$1"
                ;;

            *)
                print_cli_error "Unexpected positional argument: $1"
                exit 2
                ;;
        esac
    done
}

# Set defaults with hardware auto-detection
set_defaults() {
    SAMPLE_ID="${SAMPLE_ID:-$DEFAULT_SAMPLE_ID}"
    INPUT_DIR="${INPUT_DIR:-data/raw}"
    OUTPUT_DIR="${OUTPUT_DIR:-results}"
    DRY_RUN="${DRY_RUN:-false}"
    VERBOSE="${VERBOSE:-false}"
    ASSUME_YES="${ASSUME_YES:-false}"
    SKIP_REQUIREMENTS_CHECK="${SKIP_REQUIREMENTS_CHECK:-false}"
    
    # Apply hardware auto-detection if enabled
    if [[ "$AUTO_DETECT" == "true" ]]; then
        DETECTED_MODE=$(select_performance_mode)
        
        # Auto-detect threads (unless CLI override provided)
        if [[ -z "$CLI_THREADS" ]]; then
            THREADS=$(recommend_threads)
        else
            THREADS="$CLI_THREADS"
        fi
        
        # Auto-detect GPU settings (unless CLI override provided)
        if [[ -z "$CLI_USE_GPU" ]]; then
            local auto_gpu_settings
            read -r auto_use_gpu auto_aligner auto_gpu_count <<< "$(recommend_gpu_settings)"
            USE_GPU="$auto_use_gpu"
            GPU_ALIGNER="${GPU_ALIGNER:-$auto_aligner}"
            if [[ -z "$CLI_GPU_COUNT" ]]; then
                GPU_COUNT="$auto_gpu_count"
            else
                GPU_COUNT="$CLI_GPU_COUNT"
            fi
        else
            USE_GPU="${USE_GPU:-$DEFAULT_USE_GPU}"
            GPU_ALIGNER="${GPU_ALIGNER:-$DEFAULT_GPU_ALIGNER}"
            GPU_COUNT="${GPU_COUNT:-$DEFAULT_GPU_COUNT}"
        fi
    else
        # Auto-detection disabled: use defaults
        DETECTED_MODE="manual"
        THREADS="${THREADS:-$DEFAULT_THREADS}"
        USE_GPU="${USE_GPU:-$DEFAULT_USE_GPU}"
        GPU_ALIGNER="${GPU_ALIGNER:-$DEFAULT_GPU_ALIGNER}"
        GPU_COUNT="${GPU_COUNT:-$DEFAULT_GPU_COUNT}"
    fi

    validate_positive_integer "threads" "$THREADS"

    if [[ "$USE_GPU" == "true" ]]; then
        validate_positive_integer "gpu-count" "$GPU_COUNT"
    else
        # CPU mode: allow zero GPUs from auto-detection; only enforce numeric shape.
        if [[ -z "$GPU_COUNT" ]]; then
            GPU_COUNT=0
        fi
        if [[ ! "$GPU_COUNT" =~ ^[0-9]+$ ]]; then
            print_cli_error "gpu-count must be a non-negative integer in CPU mode (got: '$GPU_COUNT')."
            exit 2
        fi
    fi

    REGISTRY_DB="${REGISTRY_DB:-}"
}

registry_enabled() {
    [[ -n "${REGISTRY_DB:-}" ]]
}

registry_script_path() {
    echo "$SCRIPT_DIR/scripts/sample_registry.py"
}

ensure_registry_record() {
    if ! registry_enabled; then
        return 0
    fi

    local registry_script
    registry_script="$(registry_script_path)"

    if [[ ! -f "$registry_script" ]]; then
        echo "⚠️ Registry enabled but script not found: $registry_script"
        return 0
    fi

    local default_r1="$INPUT_DIR/${SAMPLE_ID}_R1.fastq.gz"
    local default_r2="$INPUT_DIR/${SAMPLE_ID}_R2.fastq.gz"

    if ! python3 "$registry_script" init --db "$REGISTRY_DB" >/dev/null 2>&1; then
        echo "⚠️ Could not initialize sample registry at $REGISTRY_DB (continuing without registry updates)"
        return 0
    fi

    if ! python3 "$registry_script" add \
        --db "$REGISTRY_DB" \
        --sample-id "$SAMPLE_ID" \
        --fastq-r1 "$default_r1" \
        --fastq-r2 "$default_r2" \
        --output-dir "$OUTPUT_DIR" \
        --status pending \
        --allow-existing \
        --quiet >/dev/null 2>&1; then
        echo "⚠️ Failed to ensure registry record for sample '$SAMPLE_ID' (continuing)"
    fi
}

update_registry_status() {
    local status="$1"
    local note="$2"

    if ! registry_enabled; then
        return 0
    fi

    local registry_script
    registry_script="$(registry_script_path)"

    if [[ ! -f "$registry_script" ]]; then
        return 0
    fi

    if ! python3 "$registry_script" update \
        --db "$REGISTRY_DB" \
        --sample-id "$SAMPLE_ID" \
        --output-dir "$OUTPUT_DIR" \
        --status "$status" \
        --notes "$note" \
        --append-notes \
        --quiet >/dev/null 2>&1; then
        echo "⚠️ Registry update failed for sample '$SAMPLE_ID' (status=$status). Continuing."
    fi
}

# Define pipeline steps
define_pipeline_steps() {
    # Define all available steps
    declare -gA PIPELINE_STEPS
    PIPELINE_STEPS["quality-control"]="scripts/quality_control.sh"
    PIPELINE_STEPS["data-cleaning"]="scripts/data_cleaning.sh"
    PIPELINE_STEPS["alignment"]="scripts/alignment.sh"
    PIPELINE_STEPS["variant-calling"]="scripts/variant_calling.sh"
    PIPELINE_STEPS["annotation"]="scripts/vep_annotation.sh"

    # Step descriptions
    declare -gA STEP_DESCRIPTIONS
    STEP_DESCRIPTIONS["quality-control"]="FastQC quality control analysis"
    STEP_DESCRIPTIONS["data-cleaning"]="Adapter trimming and quality filtering"
    STEP_DESCRIPTIONS["alignment"]="BWA alignment to reference genome"
    STEP_DESCRIPTIONS["variant-calling"]="bcftools variant calling"
    STEP_DESCRIPTIONS["annotation"]="VEP variant annotation"

    # Default step order / user-selected step order
    if [[ -n "${SPECIFIC_STEPS:-}" ]]; then
        local raw_steps=()
        IFS=',' read -r -a raw_steps <<< "$SPECIFIC_STEPS"

        STEP_ORDER=()
        local raw_step=""
        for raw_step in "${raw_steps[@]}"; do
            local step
            step="$(normalize_step_name "$raw_step")"

            if [[ -z "$step" ]]; then
                continue
            fi

            validate_step_name "$step" "step name in --steps"
            STEP_ORDER+=("$step")
        done

        if [[ ${#STEP_ORDER[@]} -eq 0 ]]; then
            print_cli_error "--steps was provided, but no valid steps were parsed."
            echo "Example: --steps quality-control,data-cleaning" >&2
            exit 2
        fi
    else
        STEP_ORDER=("${AVAILABLE_STEPS[@]}")
    fi

    # Handle resume functionality
    if [[ -n "${RESUME_FROM:-}" ]]; then
        local resume_step
        resume_step="$(normalize_step_name "$RESUME_FROM")"
        validate_step_name "$resume_step" "resume step"

        local resume_index=-1
        local i
        for i in "${!STEP_ORDER[@]}"; do
            if [[ "${STEP_ORDER[$i]}" == "$resume_step" ]]; then
                resume_index=$i
                break
            fi
        done

        if [[ $resume_index -ge 0 ]]; then
            STEP_ORDER=("${STEP_ORDER[@]:$resume_index}")
            echo "🔄 Resuming pipeline from: $resume_step"
        else
            print_cli_error "Resume step '$resume_step' is not included in selected --steps list."
            echo "Selected steps: ${STEP_ORDER[*]}" >&2
            exit 2
        fi
    fi
}

# Run a single pipeline step
run_pipeline_step() {
    local step_name="$1"
    local script_path="${PIPELINE_STEPS[$step_name]}"
    local description="${STEP_DESCRIPTIONS[$step_name]}"

    if [[ ! -f "$SCRIPT_DIR/$script_path" ]]; then
        echo -e "${RED}❌ Script not found for step '$step_name': $script_path${NC}" >&2
        return 1
    fi

    # Start step monitoring
    start_step "$step_name" "$description"

    # Build command arguments
    local cmd_args=()
    cmd_args+=(--sample-id "$SAMPLE_ID")
    cmd_args+=(--threads "$THREADS")

    # Add step-specific arguments
    case "$step_name" in
        "quality-control"|"data-cleaning")
            cmd_args+=(--input-dir "$INPUT_DIR")
            ;;
        "alignment")
            cmd_args+=(--output-dir "$OUTPUT_DIR")
            if [[ "$USE_GPU" == "true" ]]; then
                cmd_args+=(--use-gpu)
                cmd_args+=(--gpu-aligner "$GPU_ALIGNER")
                cmd_args+=(--gpu-count "$GPU_COUNT")
            fi
            ;;
        "variant-calling"|"annotation")
            cmd_args+=(--output-dir "$OUTPUT_DIR")
            ;;
    esac

    if [[ "$DRY_RUN" == "true" ]]; then
        cmd_args+=(--dry-run)
    fi

    if [[ "$VERBOSE" == "true" ]]; then
        cmd_args+=(--verbose)
    fi

    # Execute the step
    local success=true
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "🧪 Would run: $script_path ${cmd_args[*]}"
        sleep 1  # Simulate execution time while keeping tests fast
    else
        if "$SCRIPT_DIR/$script_path" "${cmd_args[@]}"; then
            :
        else
            local step_exit_code=$?
            echo -e "${RED}❌ Step '$step_name' failed with exit code $step_exit_code.${NC}" >&2
            echo "   Command: $script_path ${cmd_args[*]}" >&2
            success=false
        fi
    fi

    # Complete step monitoring
    complete_step "$step_name" "$success"

    [[ "$success" == "true" ]]
}

# Validate environment
validate_environment() {
    echo "🔍 Validating environment..."

    # Check if scripts directory exists
    if [[ ! -d "$SCRIPT_DIR/scripts" ]]; then
        echo -e "${RED}❌ Scripts directory not found: $SCRIPT_DIR/scripts${NC}" >&2
        exit 1
    fi

    # Check if required scripts exist
    local missing_scripts=()
    local step=""
    for step in "${STEP_ORDER[@]}"; do
        local script_path="${PIPELINE_STEPS[$step]}"
        if [[ ! -f "$SCRIPT_DIR/$script_path" ]]; then
            missing_scripts+=("$script_path")
        fi
    done

    if [[ ${#missing_scripts[@]} -gt 0 ]]; then
        echo -e "${RED}❌ Missing scripts:${NC}" >&2
        printf '  %s\n' "${missing_scripts[@]}" >&2
        exit 1
    fi

    # Validate input directory if early steps are requested.
    local requires_raw_input=false
    for step in "${STEP_ORDER[@]}"; do
        if [[ "$step" == "quality-control" || "$step" == "data-cleaning" ]]; then
            requires_raw_input=true
            break
        fi
    done

    if [[ "$requires_raw_input" == "true" && ! -d "$INPUT_DIR" ]]; then
        if [[ "$DRY_RUN" == "true" ]]; then
            echo -e "${YELLOW}⚠️ Input directory does not exist (dry-run continues): $INPUT_DIR${NC}"
        else
            echo -e "${RED}❌ Input directory not found: $INPUT_DIR${NC}" >&2
            echo "   Provide --input-dir with FASTQ files, or run --dry-run to validate config only." >&2
            exit 1
        fi
    fi

    # Run system requirements check if available
    if [[ "$SKIP_REQUIREMENTS_CHECK" == "true" ]]; then
        echo "⚠️ Skipping system requirements check (--skip-requirements-check)"
    elif [[ -f "$SCRIPT_DIR/scripts/check_requirements.sh" ]]; then
        local req_args=(--min-ram 16 --min-disk 400)

        # Dry-run should not require full bioinformatics stack.
        if [[ "$DRY_RUN" == "true" ]]; then
            req_args+=(--skip-conda --skip-env --skip-tools)
        fi

        if ! "$SCRIPT_DIR/scripts/check_requirements.sh" "${req_args[@]}"; then
            echo -e "${YELLOW}⚠️ System requirements check failed${NC}"

            if [[ "$ASSUME_YES" == "true" ]]; then
                echo "⚠️ Continuing due to --yes"
            elif [[ ! -t 0 ]]; then
                echo -e "${RED}❌ Non-interactive session:${NC} rerun with --yes or --skip-requirements-check to continue." >&2
                exit 1
            else
                read -r -p "Continue anyway? (y/N): "
                if [[ ! "$REPLY" =~ ^[Yy]$ ]]; then
                    exit 1
                fi
            fi
        fi
    fi

    echo "✅ Environment validation passed"
}

# Show pipeline summary
show_pipeline_summary() {
    echo
    echo "📋 Pipeline Configuration Summary:"
    echo "=================================="
    
    # Show hardware detection status
    if [[ "$AUTO_DETECT" == "true" ]]; then
        local hw_cores hw_ram hw_gpus
        hw_cores=$(detect_cpu_cores)
        hw_ram=$(detect_ram_gb)
        hw_gpus=$(detect_gpu_count)
        
        printf "🔍 Hardware Detected:\n"
        printf "   CPU: %s cores | RAM: %sGB | GPUs: %s\n" "$hw_cores" "$hw_ram" "$hw_gpus"
        printf "   Performance Mode: %s\n" "$DETECTED_MODE"
        
        # Show which settings are auto vs manual
        local thread_source="auto"
        local gpu_source="auto"
        [[ -n "$CLI_THREADS" ]] && thread_source="manual"
        [[ -n "$CLI_USE_GPU" ]] && gpu_source="manual"
        printf "   Threads: %s (%s) | GPU: %s (%s)\n" "$THREADS" "$thread_source" "$USE_GPU" "$gpu_source"
        echo
    fi
    
    printf "Sample ID: %s\n" "$SAMPLE_ID"
    printf "Threads: %s\n" "$THREADS"
    printf "Input Directory: %s\n" "$INPUT_DIR"
    printf "Output Directory: %s\n" "$OUTPUT_DIR"
    printf "Steps to Run: %s\n" "${STEP_ORDER[*]}"

    if [[ "$USE_GPU" == "true" ]]; then
        printf "Alignment Mode: GPU (%s, %d GPUs)\n" "$GPU_ALIGNER" "$GPU_COUNT"
    else
        printf "Alignment Mode: CPU (BWA)\n"
    fi

    if [[ "$ENABLE_NOTIFICATIONS" == "true" ]]; then
        printf "Notifications: Enabled"
        [[ -n "$EMAIL_ADDRESS" ]] && printf " (Email: %s)" "$EMAIL_ADDRESS"
        [[ -n "$SLACK_WEBHOOK" ]] && printf " (Slack: Enabled)"
        printf "\n"
    fi

    if [[ "$DRY_RUN" == "true" ]]; then
        printf "Mode: DRY RUN (no actual processing)\n"
    fi

    if [[ "$SKIP_REQUIREMENTS_CHECK" == "true" ]]; then
        printf "Requirements Check: Skipped\n"
    fi

    if registry_enabled; then
        printf "Sample Registry DB: %s\n" "$REGISTRY_DB"
    else
        printf "Sample Registry DB: disabled\n"
    fi

    echo
}

handle_interrupt() {
    echo
    echo "⚠️ Pipeline interrupted by user"
    if declare -F cleanup_progress_monitor >/dev/null 2>&1; then
        cleanup_progress_monitor
    fi
    exit 1
}

# Main execution
main() {
    parse_arguments "$@"
    set_defaults

    # Require bash4+ only for actual execution path (help works everywhere).
    require_bash4_or_die

    # Source progress monitoring after bash compatibility check.
    # shellcheck disable=SC1091
    source "$SCRIPT_DIR/scripts/progress_monitor.sh"

    define_pipeline_steps

    # Show summary
    show_pipeline_summary

    # Validate environment
    validate_environment

    # Optional sample-registry setup (non-blocking)
    ensure_registry_record
    update_registry_status "running" "Pipeline started"

    # Initialize progress monitoring
    init_progress_monitor "WGS_Pipeline_${SAMPLE_ID}" "${#STEP_ORDER[@]}" "logs"

    # Enable notifications if configured
    if [[ "$ENABLE_NOTIFICATIONS" == "true" ]]; then
        enable_notifications "$EMAIL_ADDRESS" "$SLACK_WEBHOOK"
        send_notification "Started" "WGS Pipeline started for sample: $SAMPLE_ID" "info"
    fi

    # Run pipeline steps
    local pipeline_success=true
    local failed_step=""

    local step=""
    for step in "${STEP_ORDER[@]}"; do
        if ! run_pipeline_step "$step"; then
            pipeline_success=false
            failed_step="$step"
            break
        fi
    done

    # Show resource summary
    show_resource_summary

    # Finish monitoring
    finish_progress_monitor "$pipeline_success"

    # Final status
    if [[ "$pipeline_success" == "true" ]]; then
        update_registry_status "completed" "Pipeline finished successfully"
        echo "🎉 Pipeline completed successfully!"
        echo "📁 Results available in: $OUTPUT_DIR"
        exit 0
    else
        update_registry_status "failed" "Pipeline failed at step: $failed_step"
        echo "💥 Pipeline failed at step: $failed_step"
        echo "📋 Check logs for details: logs/WGS_Pipeline_${SAMPLE_ID}_progress.log"
        exit 1
    fi
}

# Handle interrupts
trap handle_interrupt INT TERM

trap 'echo; echo "⚠️ Pipeline interrupted by user"; update_registry_status "failed" "Pipeline interrupted by user"; cleanup_progress_monitor; exit 1' INT TERM

# Run main function
main "$@"
