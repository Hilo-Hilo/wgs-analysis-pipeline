#!/bin/bash

# Configuration Loader for WGS Analysis Pipeline
# This utility script loads and validates configuration files
# Source this script from other pipeline scripts to load configuration

# DO NOT RUN THIS SCRIPT DIRECTLY - It should be sourced by other scripts
# Usage: source scripts/load_config.sh [CONFIG_FILE]

# =============================================================================
# CONFIGURATION LOADING FUNCTIONS
# =============================================================================

# Default configuration file
DEFAULT_CONFIG_FILE="config/default.conf"

# Load configuration from file
load_config() {
    local config_file="$1"
    
    # Use default config if none specified
    if [[ -z "$config_file" ]]; then
        config_file="$DEFAULT_CONFIG_FILE"
    fi
    
    # Check if config file exists
    if [[ ! -f "$config_file" ]]; then
        echo "ERROR: Configuration file not found: $config_file" >&2
        echo "Available configurations:" >&2
        ls -la config/*.conf 2>/dev/null >&2 || echo "No configuration files found in config/" >&2
        return 1
    fi
    
    # Validate config file before sourcing
    if ! validate_config_syntax "$config_file"; then
        echo "ERROR: Invalid configuration syntax in: $config_file" >&2
        return 1
    fi
    
    # Source the configuration file
    echo "Loading configuration: $config_file" >&2
    source "$config_file"
    
    # Set derived variables
    set_derived_variables
    
    # Validate configuration values
    if ! validate_config_values; then
        echo "ERROR: Configuration validation failed" >&2
        return 1
    fi
    
    echo "Configuration loaded successfully" >&2
    return 0
}

# Validate configuration file syntax
validate_config_syntax() {
    local config_file="$1"
    
    # Check for common syntax errors
    if grep -E '^[^#]*[[:space:]]+=' "$config_file" >/dev/null 2>&1; then
        echo "ERROR: Spaces around '=' found in config file (should be VAR=value)" >&2
        grep -n -E '^[^#]*[[:space:]]+=' "$config_file" >&2
        return 1
    fi
    
    # Try to source in a subshell to catch syntax errors
    if ! bash -n "$config_file" 2>/dev/null; then
        echo "ERROR: Bash syntax errors in config file" >&2
        bash -n "$config_file" 2>&1 >&2
        return 1
    fi
    
    return 0
}

# Set derived variables based on configuration
set_derived_variables() {
    # Set reference genome path (GRCh38 only)
    case "${PRIMARY_REFERENCE:-grch38}" in
        "grch38")
            REFERENCE_GENOME="${GRCH38_REFERENCE:-data/reference/GRCh38/GRCh38_latest_genomic.fna}"
            REFERENCE_GTF="${GRCH38_GTF:-data/reference/GRCh38/GRCh38_latest_genomic.gtf}"
            ;;
        *)
            echo "ERROR: Only GRCh38 reference is supported, got: ${PRIMARY_REFERENCE}" >&2
            return 1
            ;;
    esac
    
    # Set sample-specific output paths
    if [[ -n "${SAMPLE_ID:-}" ]]; then
        SAMPLE_OUTPUT_DIR="${RESULTS_DIR}/${SAMPLE_ID}"
        SAMPLE_LOG_DIR="${LOG_DIR}/${SAMPLE_ID}"
    fi
    
    # Set tool-specific memory limits with defaults
    BWA_MEMORY="${MAX_MEMORY_BWA:-16}G"
    SAMTOOLS_MEMORY="${MAX_MEMORY_SAMTOOLS:-8}G"
    BCFTOOLS_MEMORY="${MAX_MEMORY_BCFTOOLS:-8}G"
    VEP_MEMORY="${MAX_MEMORY_VEP:-16}G"
    
    # Convert boolean strings to bash-friendly format
    normalize_boolean_vars
}

# Normalize boolean variables
normalize_boolean_vars() {
    local bool_vars=(
        "KEEP_INTERMEDIATE_FILES"
        "GENERATE_DETAILED_STATS"
        "DEFAULT_DRY_RUN"
        "DEFAULT_VERBOSE"
        "DEFAULT_FORCE_OVERWRITE"
        "AUTO_CLEANUP"
        "ENABLE_CHECKPOINTING"
        "VALIDATE_INPUT_FILES"
        "CHECK_FILE_FORMATS"
        "VERIFY_TOOL_VERSIONS"
        "ENABLE_EMAIL_NOTIFICATIONS"
    )
    
    for var in "${bool_vars[@]}"; do
        local value="${!var:-false}"
        # Convert to lowercase for bash compatibility
        local lower_value=$(echo "$value" | tr '[:upper:]' '[:lower:]')
        case "$lower_value" in
            "true"|"yes"|"1"|"on")
                eval "$var=true"
                ;;
            "false"|"no"|"0"|"off"|"")
                eval "$var=false"
                ;;
            *)
                echo "WARNING: Invalid boolean value for $var: $value (defaulting to false)" >&2
                eval "$var=false"
                ;;
        esac
    done
}

# Validate configuration values
validate_config_values() {
    local validation_failed=false
    
    # Check required variables
    local required_vars=(
        "CONDA_ENV"
        "THREADS"
        "RAW_DATA_DIR"
        "RESULTS_DIR"
        "LOG_DIR"
    )
    
    for var in "${required_vars[@]}"; do
        if [[ -z "${!var:-}" ]]; then
            echo "ERROR: Required configuration variable not set: $var" >&2
            validation_failed=true
        fi
    done
    
    # Validate numeric values
    if ! [[ "${THREADS:-}" =~ ^[0-9]+$ ]] || [[ "${THREADS:-0}" -le 0 ]]; then
        echo "ERROR: THREADS must be a positive integer, got: ${THREADS:-}" >&2
        validation_failed=true
    fi
    
    if ! [[ "${MIN_FREE_SPACE_GB:-100}" =~ ^[0-9]+$ ]]; then
        echo "ERROR: MIN_FREE_SPACE_GB must be a positive integer, got: ${MIN_FREE_SPACE_GB:-}" >&2
        validation_failed=true
    fi
    
    # Validate reference genome choice
    case "${PRIMARY_REFERENCE:-}" in
        "grch38"|"")
            ;;
        *)
            echo "ERROR: PRIMARY_REFERENCE must be 'grch38', got: ${PRIMARY_REFERENCE}" >&2
            validation_failed=true
            ;;
    esac
    
    # Check if critical directories exist or can be created
    local dirs_to_check=(
        "${RESULTS_DIR}"
        "${LOG_DIR}"
    )
    
    for dir in "${dirs_to_check[@]}"; do
        if [[ -n "$dir" ]]; then
            if ! mkdir -p "$dir" 2>/dev/null; then
                echo "ERROR: Cannot create directory: $dir" >&2
                validation_failed=true
            fi
        fi
    done
    
    if [[ "$validation_failed" == "true" ]]; then
        return 1
    fi
    
    return 0
}

# Show current configuration
show_config() {
    echo "Current WGS Analysis Configuration:"
    echo "=================================="
    echo ""
    echo "Environment:"
    echo "  Conda Environment: ${CONDA_ENV:-not set}"
    echo "  Threads: ${THREADS:-not set}"
    echo "  Primary Reference: ${PRIMARY_REFERENCE:-not set}"
    echo "  Sample ID: ${SAMPLE_ID:-not set}"
    echo ""
    echo "Directories:"
    echo "  Raw Data: ${RAW_DATA_DIR:-not set}"
    echo "  Processed Data: ${PROCESSED_DATA_DIR:-not set}"
    echo "  Reference: ${REFERENCE_DIR:-not set}"
    echo "  Results: ${RESULTS_DIR:-not set}"
    echo "  Logs: ${LOG_DIR:-not set}"
    echo ""
    echo "Reference Genome:"
    echo "  File: ${REFERENCE_GENOME:-not set}"
    echo "  GTF: ${REFERENCE_GTF:-not set}"
    echo ""
    echo "Resource Limits:"
    echo "  BWA Memory: ${BWA_MEMORY:-not set}"
    echo "  VEP Memory: ${VEP_MEMORY:-not set}"
    echo "  Min Free Space: ${MIN_FREE_SPACE_GB:-not set}GB"
    echo ""
    echo "Workflow Options:"
    echo "  Keep Intermediate Files: ${KEEP_INTERMEDIATE_FILES:-not set}"
    echo "  Generate Detailed Stats: ${GENERATE_DETAILED_STATS:-not set}"
    echo "  Auto Cleanup: ${AUTO_CLEANUP:-not set}"
    echo "  Enable Checkpointing: ${ENABLE_CHECKPOINTING:-not set}"
    echo ""
}

# List available configuration files
list_configs() {
    echo "Available Configuration Files:"
    echo "=============================="
    
    if [[ -d "config" ]]; then
        for config_file in config/*.conf; do
            if [[ -f "$config_file" ]]; then
                local basename
                basename=$(basename "$config_file")
                echo "  $basename"
                
                # Show description if available
                local description
                description=$(grep "^# .*Configuration" "$config_file" | head -1 | sed 's/^# *//')
                if [[ -n "$description" ]]; then
                    echo "    Description: $description"
                fi
                
                # Show key settings
                local sample_id threads
                sample_id=$(grep "^SAMPLE_ID=" "$config_file" | cut -d'=' -f2 | tr -d '"' | head -1)
                threads=$(grep "^THREADS=" "$config_file" | cut -d'=' -f2 | tr -d '"' | head -1)
                
                if [[ -n "$sample_id" || -n "$threads" ]]; then
                    echo -n "    Settings:"
                    [[ -n "$sample_id" ]] && echo -n " Sample=$sample_id"
                    [[ -n "$threads" ]] && echo -n " Threads=$threads"
                    echo ""
                fi
                echo ""
            fi
        done
    else
        echo "  No config directory found"
    fi
}

# Create a new configuration file from template
create_config() {
    local config_name="$1"
    
    if [[ -z "$config_name" ]]; then
        echo "ERROR: Please specify a configuration name" >&2
        echo "Usage: create_config <name>" >&2
        echo "Example: create_config my_analysis" >&2
        return 1
    fi
    
    local new_config="config/${config_name}.conf"
    
    if [[ -f "$new_config" ]]; then
        echo "ERROR: Configuration already exists: $new_config" >&2
        echo "Remove existing file or choose a different name" >&2
        return 1
    fi
    
    # Copy from example configuration
    if [[ -f "config/example.conf" ]]; then
        cp "config/example.conf" "$new_config"
        echo "Created new configuration: $new_config"
        echo "Edit this file to customize your analysis settings"
    else
        echo "ERROR: Template configuration not found: config/example.conf" >&2
        return 1
    fi
}

# =============================================================================
# COMMAND LINE INTERFACE
# =============================================================================

# If script is run directly (not sourced), provide CLI interface
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    case "${1:-}" in
        "show"|"display")
            if [[ -n "${2:-}" ]]; then
                load_config "$2" && show_config
            else
                load_config && show_config
            fi
            ;;
        "list"|"ls")
            list_configs
            ;;
        "create")
            create_config "${2:-}"
            ;;
        "validate")
            config_file="${2:-$DEFAULT_CONFIG_FILE}"
            echo "Validating configuration: $config_file"
            if load_config "$config_file"; then
                echo "✓ Configuration is valid"
            else
                echo "✗ Configuration validation failed"
                exit 1
            fi
            ;;
        "help"|"-h"|"--help"|"")
            cat << EOF
Configuration Loader Utility

USAGE:
    source scripts/load_config.sh [CONFIG_FILE]  # Load config in current shell
    scripts/load_config.sh COMMAND [OPTIONS]     # Run utility commands

COMMANDS:
    show [CONFIG_FILE]     Show current or specified configuration
    list                   List available configuration files
    create <NAME>          Create new configuration from template
    validate [CONFIG_FILE] Validate configuration syntax and values
    help                   Show this help message

EXAMPLES:
    # Load default configuration in a script
    source scripts/load_config.sh

    # Load specific configuration
    source scripts/load_config.sh config/my_analysis.conf

    # Show current configuration
    scripts/load_config.sh show

    # Create new configuration
    scripts/load_config.sh create my_project

    # Validate configuration
    scripts/load_config.sh validate config/my_analysis.conf

FILES:
    config/default.conf    Default pipeline settings
    config/example.conf    Template for custom configurations
    config/*.conf          User-specific configurations

EOF
            ;;
        *)
            echo "Unknown command: $1" >&2
            echo "Use 'help' for usage information" >&2
            exit 1
            ;;
    esac
fi