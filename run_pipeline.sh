#!/bin/bash

# WGS Pipeline Runner with Progress Monitoring
# Enhanced pipeline execution with real-time progress tracking and resource monitoring

set -e

# Require Bash 4+ (associative arrays used by progress monitor)
if [[ -z "${BASH_VERSINFO:-}" || "${BASH_VERSINFO[0]}" -lt 4 ]]; then
    echo "❌ run_pipeline.sh requires Bash 4 or newer."
    echo "   Detected: ${BASH_VERSION:-unknown}"
    echo "   On macOS, install newer bash (e.g. via Homebrew) and run with that binary."
    exit 1
fi

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source progress monitoring
source "$SCRIPT_DIR/scripts/progress_monitor.sh"

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

# Help function
show_help() {
    cat << EOF
WGS Pipeline Runner with Progress Monitoring

USAGE:
    $0 [OPTIONS]

OPTIONS:
    -h, --help                  Show this help message
    -s, --sample-id ID          Sample identifier (default: WGS_SAMPLE)
    -t, --threads NUM           Number of threads (default: 4)
    -i, --input-dir DIR         Input directory with raw FASTQ files
    -o, --output-dir DIR        Output directory (default: results)
    --use-gpu                   Use GPU alignment for alignment step
    --gpu-aligner NAME          GPU aligner backend (default: parabricks)
    --gpu-count NUM             Number of GPUs for alignment step (default: 1)
    --resume-from STEP          Resume pipeline from specific step
    --notify-email EMAIL        Send email notifications
    --notify-slack WEBHOOK      Send Slack notifications
    --steps STEPS               Run specific steps only (comma-separated)
    -y, --yes                   Continue on requirement warnings without prompting
    --skip-requirements-check   Skip running scripts/check_requirements.sh
    --registry-db PATH          Optional SQLite sample registry path
    --dry-run                   Show what would be done without executing
    --verbose                   Enable verbose output

PIPELINE STEPS:
    1. quality-control          FastQC analysis of raw reads
    2. data-cleaning           Adapter trimming and quality filtering
    3. alignment               BWA mapping to reference genome
    4. variant-calling         bcftools variant calling
    5. annotation              VEP variant annotation

EXAMPLES:
    # Run complete pipeline
    $0 --sample-id MySample --input-dir data/raw

    # Run with notifications
    $0 --sample-id MySample --notify-email user@example.com

    # Resume from alignment step
    $0 --resume-from alignment

    # Use GPU alignment on DGX (alignment step only)
    $0 --sample-id MySample --input-dir data/raw --use-gpu --gpu-aligner parabricks --gpu-count 1

    # Run specific steps only
    $0 --steps quality-control,data-cleaning

    # Dry run to check configuration
    $0 --dry-run --verbose
EOF
}

# Parse command line arguments
parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                show_help
                exit 0
                ;;
            -s|--sample-id)
                SAMPLE_ID="$2"
                shift 2
                ;;
            -t|--threads)
                THREADS="$2"
                shift 2
                ;;
            -i|--input-dir)
                INPUT_DIR="$2"
                shift 2
                ;;
            -o|--output-dir)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            --use-gpu)
                USE_GPU=true
                shift
                ;;
            --gpu-aligner)
                GPU_ALIGNER="$2"
                shift 2
                ;;
            --gpu-count)
                GPU_COUNT="$2"
                shift 2
                ;;
            --resume-from)
                RESUME_FROM="$2"
                shift 2
                ;;
            --notify-email)
                EMAIL_ADDRESS="$2"
                ENABLE_NOTIFICATIONS=true
                shift 2
                ;;
            --notify-slack)
                SLACK_WEBHOOK="$2"
                ENABLE_NOTIFICATIONS=true
                shift 2
                ;;
            --steps)
                SPECIFIC_STEPS="$2"
                shift 2
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
            *)
                echo "Unknown option: $1"
                echo "Use --help for usage information"
                exit 1
                ;;
        esac
    done
}

# Set defaults
set_defaults() {
    SAMPLE_ID="${SAMPLE_ID:-$DEFAULT_SAMPLE_ID}"
    THREADS="${THREADS:-$DEFAULT_THREADS}"
    INPUT_DIR="${INPUT_DIR:-data/raw}"
    OUTPUT_DIR="${OUTPUT_DIR:-results}"
    USE_GPU="${USE_GPU:-$DEFAULT_USE_GPU}"
    GPU_ALIGNER="${GPU_ALIGNER:-$DEFAULT_GPU_ALIGNER}"
    GPU_COUNT="${GPU_COUNT:-$DEFAULT_GPU_COUNT}"
    DRY_RUN="${DRY_RUN:-false}"
    VERBOSE="${VERBOSE:-false}"
    ASSUME_YES="${ASSUME_YES:-false}"
    SKIP_REQUIREMENTS_CHECK="${SKIP_REQUIREMENTS_CHECK:-false}"
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
    
    # Default step order
    if [[ -n "$SPECIFIC_STEPS" ]]; then
        IFS=',' read -ra STEP_ORDER <<< "$SPECIFIC_STEPS"
    else
        STEP_ORDER=("quality-control" "data-cleaning" "alignment" "variant-calling" "annotation")
    fi
    
    # Handle resume functionality
    if [[ -n "$RESUME_FROM" ]]; then
        local resume_index=-1
        for i in "${!STEP_ORDER[@]}"; do
            if [[ "${STEP_ORDER[$i]}" == "$RESUME_FROM" ]]; then
                resume_index=$i
                break
            fi
        done
        
        if [[ $resume_index -ge 0 ]]; then
            STEP_ORDER=("${STEP_ORDER[@]:$resume_index}")
            echo "🔄 Resuming pipeline from: $RESUME_FROM"
        else
            echo "❌ Invalid resume step: $RESUME_FROM"
            echo "Available steps: ${!PIPELINE_STEPS[*]}"
            exit 1
        fi
    fi
}

# Run a single pipeline step
run_pipeline_step() {
    local step_name="$1"
    local script_path="${PIPELINE_STEPS[$step_name]}"
    local description="${STEP_DESCRIPTIONS[$step_name]}"
    
    if [[ ! -f "$SCRIPT_DIR/$script_path" ]]; then
        echo "❌ Script not found: $script_path"
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
        sleep 2  # Simulate execution time
    else
        if ! "$SCRIPT_DIR/$script_path" "${cmd_args[@]}"; then
            success=false
        fi
    fi
    
    # Complete step monitoring
    complete_step "$step_name" "$success"
    
    return $([ "$success" = true ] && echo 0 || echo 1)
}

# Validate environment
validate_environment() {
    echo "🔍 Validating environment..."
    
    # Check if scripts directory exists
    if [[ ! -d "$SCRIPT_DIR/scripts" ]]; then
        echo "❌ Scripts directory not found: $SCRIPT_DIR/scripts"
        exit 1
    fi
    
    # Check if required scripts exist
    local missing_scripts=()
    for step in "${STEP_ORDER[@]}"; do
        local script_path="${PIPELINE_STEPS[$step]}"
        if [[ ! -f "$SCRIPT_DIR/$script_path" ]]; then
            missing_scripts+=("$script_path")
        fi
    done
    
    if [[ ${#missing_scripts[@]} -gt 0 ]]; then
        echo "❌ Missing scripts:"
        printf '  %s\n' "${missing_scripts[@]}"
        exit 1
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
            echo "⚠️ System requirements check failed"

            if [[ "$ASSUME_YES" == "true" ]]; then
                echo "⚠️ Continuing due to --yes"
            elif [[ ! -t 0 ]]; then
                echo "❌ Non-interactive session: rerun with --yes or --skip-requirements-check to continue."
                exit 1
            else
                read -p "Continue anyway? (y/N): " -r
                if [[ ! $REPLY =~ ^[Yy]$ ]]; then
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
    printf "Sample ID: %s\n" "$SAMPLE_ID"
    printf "Threads: %s\n" "$THREADS"
    printf "Input Directory: %s\n" "$INPUT_DIR"
    printf "Output Directory: %s\n" "$OUTPUT_DIR"
    printf "Steps to Run: %s\n" "${STEP_ORDER[*]}"
    if [[ "$USE_GPU" == "true" ]]; then
        printf "Alignment Mode: GPU (%s, %s GPU)\n" "$GPU_ALIGNER" "$GPU_COUNT"
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

# Main execution
main() {
    parse_arguments "$@"
    set_defaults
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
trap 'echo; echo "⚠️ Pipeline interrupted by user"; update_registry_status "failed" "Pipeline interrupted by user"; cleanup_progress_monitor; exit 1' INT TERM

# Run main function
main "$@"
