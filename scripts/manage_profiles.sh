#!/bin/bash

# WGS Pipeline Profile Manager
# Manages configuration profiles for different system types

set -e

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
PROFILES_DIR="$PROJECT_ROOT/config/profiles"
CONFIG_DIR="$PROJECT_ROOT/config"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Help function
show_help() {
    cat << EOF
WGS Pipeline Profile Manager

DESCRIPTION:
    Manages configuration profiles for different system types.
    Profiles optimize pipeline settings for specific hardware configurations.

USAGE:
    $0 <command> [options]

COMMANDS:
    list                    List all available profiles
    show <profile>          Show profile details
    set <profile>           Set active profile
    create <name>           Create new profile from template
    validate <profile>      Validate profile configuration
    compare <prof1> <prof2> Compare two profiles
    reset                   Reset to default configuration
    help                    Show this help message

AVAILABLE PROFILES:
    laptop         - 8-16GB RAM, 2-4 cores (conservative settings)
    workstation    - 32-64GB RAM, 8-16 cores (balanced settings)
    server         - 128GB+ RAM, 32+ cores (high performance)
    cloud          - Variable resources (cost-optimized)

EXAMPLES:
    # List available profiles
    $0 list

    # Show profile details
    $0 show laptop

    # Set active profile
    $0 set workstation

    # Validate profile
    $0 validate server

    # Compare profiles
    $0 compare laptop server

    # Create custom profile
    $0 create my_system

EOF
}

# Logging functions
log() {
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1"
}

error() {
    echo -e "${RED}[$(date '+%Y-%m-%d %H:%M:%S')] ERROR:${NC} $1" >&2
}

warning() {
    echo -e "${YELLOW}[$(date '+%Y-%m-%d %H:%M:%S')] WARNING:${NC} $1"
}

info() {
    echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')] INFO:${NC} $1"
}

# List available profiles
list_profiles() {
    log "Available configuration profiles:"
    echo
    
    if [[ ! -d "$PROFILES_DIR" ]]; then
        error "Profiles directory not found: $PROFILES_DIR"
        return 1
    fi
    
    local profiles=($(find "$PROFILES_DIR" -name "*.conf" -exec basename {} .conf \; | sort))
    
    if [[ ${#profiles[@]} -eq 0 ]]; then
        warning "No profiles found in $PROFILES_DIR"
        return 0
    fi
    
    printf "%-15s %-40s %-15s\n" "Profile" "Description" "Resources"
    printf "%s\n" "$(printf '=%.0s' {1..75})"
    
    for profile in "${profiles[@]}"; do
        local profile_file="$PROFILES_DIR/${profile}.conf"
        local description=""
        local resources=""
        
        if [[ -f "$profile_file" ]]; then
            description=$(grep "^PROFILE_DESCRIPTION=" "$profile_file" 2>/dev/null | cut -d'"' -f2 || echo "No description")
            local max_mem=$(grep "^MAX_MEMORY_GB=" "$profile_file" 2>/dev/null | cut -d'=' -f2 || echo "?")
            local max_threads=$(grep "^MAX_THREADS=" "$profile_file" 2>/dev/null | cut -d'=' -f2 || echo "?")
            resources="${max_mem}GB, ${max_threads} threads"
        fi
        
        printf "%-15s %-40s %-15s\n" "$profile" "$description" "$resources"
    done
    
    echo
    
    # Show current active profile
    local current_profile=""
    if [[ -f "$CONFIG_DIR/current_profile.conf" ]]; then
        current_profile=$(cat "$CONFIG_DIR/current_profile.conf" 2>/dev/null || echo "none")
        info "Current active profile: $current_profile"
    else
        info "No profile currently active (using default configuration)"
    fi
}

# Show profile details
show_profile() {
    local profile="$1"
    
    if [[ -z "$profile" ]]; then
        error "Profile name required"
        return 1
    fi
    
    local profile_file="$PROFILES_DIR/${profile}.conf"
    
    if [[ ! -f "$profile_file" ]]; then
        error "Profile not found: $profile"
        echo "Available profiles: $(find "$PROFILES_DIR" -name "*.conf" -exec basename {} .conf \; | tr '\n' ' ')"
        return 1
    fi
    
    log "Profile details: $profile"
    echo
    
    # Show key settings
    local settings=(
        "PROFILE_DESCRIPTION:Description"
        "MAX_MEMORY_GB:Max Memory (GB)"
        "MAX_THREADS:Max Threads"
        "FASTP_THREADS:Fastp Threads"
        "BWA_THREADS:BWA Threads"
        "BWA_MEM_LIMIT:BWA Memory Limit"
        "FASTP_QUALITY_THRESHOLD:Quality Threshold"
        "ENABLE_COMPRESSION:Compression"
        "CLEANUP_INTERMEDIATE:Cleanup Intermediate"
    )
    
    printf "%-25s %s\n" "Setting" "Value"
    printf "%s\n" "$(printf '=%.0s' {1..50})"
    
    for setting_info in "${settings[@]}"; do
        local key="${setting_info%:*}"
        local label="${setting_info#*:}"
        local value=$(grep "^${key}=" "$profile_file" 2>/dev/null | cut -d'=' -f2- | sed 's/"//g' || echo "not set")
        printf "%-25s %s\n" "$label" "$value"
    done
    
    echo
    info "Full configuration file: $profile_file"
}

# Set active profile
set_profile() {
    local profile="$1"
    
    if [[ -z "$profile" ]]; then
        error "Profile name required"
        return 1
    fi
    
    local profile_file="$PROFILES_DIR/${profile}.conf"
    
    if [[ ! -f "$profile_file" ]]; then
        error "Profile not found: $profile"
        return 1
    fi
    
    # Validate profile first
    if ! validate_profile "$profile" --quiet; then
        error "Profile validation failed. Not setting as active."
        return 1
    fi
    
    # Create config directory if needed
    mkdir -p "$CONFIG_DIR"
    
    # Copy profile to active configuration
    cp "$profile_file" "$CONFIG_DIR/active.conf"
    echo "$profile" > "$CONFIG_DIR/current_profile.conf"
    
    log "Active profile set to: $profile"
    info "Configuration applied: $CONFIG_DIR/active.conf"
}

# Validate profile
validate_profile() {
    local profile="$1"
    local quiet="$2"
    
    if [[ -z "$profile" ]]; then
        error "Profile name required"
        return 1
    fi
    
    local profile_file="$PROFILES_DIR/${profile}.conf"
    
    if [[ ! -f "$profile_file" ]]; then
        error "Profile not found: $profile"
        return 1
    fi
    
    [[ "$quiet" != "--quiet" ]] && log "Validating profile: $profile"
    
    local validation_errors=0
    
    # Required settings
    local required_settings=(
        "PROFILE_NAME"
        "MAX_MEMORY_GB"
        "MAX_THREADS"
        "FASTP_THREADS"
        "BWA_THREADS"
    )
    
    for setting in "${required_settings[@]}"; do
        if ! grep -q "^${setting}=" "$profile_file"; then
            [[ "$quiet" != "--quiet" ]] && error "Missing required setting: $setting"
            ((validation_errors++))
        fi
    done
    
    # Validate numeric values
    local numeric_settings=(
        "MAX_MEMORY_GB"
        "MAX_THREADS"
        "FASTP_THREADS"
        "BWA_THREADS"
    )
    
    for setting in "${numeric_settings[@]}"; do
        local value=$(grep "^${setting}=" "$profile_file" 2>/dev/null | cut -d'=' -f2)
        if [[ -n "$value" && ! "$value" =~ ^[0-9]+$ ]]; then
            [[ "$quiet" != "--quiet" ]] && error "Invalid numeric value for $setting: $value"
            ((validation_errors++))
        fi
    done
    
    # Check logical constraints
    local max_threads=$(grep "^MAX_THREADS=" "$profile_file" 2>/dev/null | cut -d'=' -f2 || echo "0")
    local bwa_threads=$(grep "^BWA_THREADS=" "$profile_file" 2>/dev/null | cut -d'=' -f2 || echo "0")
    
    if [[ "$bwa_threads" -gt "$max_threads" ]]; then
        [[ "$quiet" != "--quiet" ]] && warning "BWA_THREADS ($bwa_threads) exceeds MAX_THREADS ($max_threads)"
    fi
    
    if [[ $validation_errors -eq 0 ]]; then
        [[ "$quiet" != "--quiet" ]] && log "Profile validation passed: $profile"
        return 0
    else
        [[ "$quiet" != "--quiet" ]] && error "Profile validation failed: $validation_errors errors found"
        return 1
    fi
}

# Compare profiles
compare_profiles() {
    local profile1="$1"
    local profile2="$2"
    
    if [[ -z "$profile1" || -z "$profile2" ]]; then
        error "Two profile names required"
        return 1
    fi
    
    local profile1_file="$PROFILES_DIR/${profile1}.conf"
    local profile2_file="$PROFILES_DIR/${profile2}.conf"
    
    if [[ ! -f "$profile1_file" ]]; then
        error "Profile not found: $profile1"
        return 1
    fi
    
    if [[ ! -f "$profile2_file" ]]; then
        error "Profile not found: $profile2"
        return 1
    fi
    
    log "Comparing profiles: $profile1 vs $profile2"
    echo
    
    # Key settings to compare
    local settings=(
        "MAX_MEMORY_GB"
        "MAX_THREADS"
        "FASTP_THREADS"
        "BWA_THREADS"
        "BWA_MEM_LIMIT"
        "FASTP_QUALITY_THRESHOLD"
        "ENABLE_COMPRESSION"
        "CLEANUP_INTERMEDIATE"
    )
    
    printf "%-25s %-15s %-15s %s\n" "Setting" "$profile1" "$profile2" "Difference"
    printf "%s\n" "$(printf '=%.0s' {1..65})"
    
    for setting in "${settings[@]}"; do
        local value1=$(grep "^${setting}=" "$profile1_file" 2>/dev/null | cut -d'=' -f2- | sed 's/"//g' || echo "unset")
        local value2=$(grep "^${setting}=" "$profile2_file" 2>/dev/null | cut -d'=' -f2- | sed 's/"//g' || echo "unset")
        
        local diff=""
        if [[ "$value1" != "$value2" ]]; then
            diff="DIFFERENT"
        else
            diff="same"
        fi
        
        printf "%-25s %-15s %-15s %s\n" "$setting" "$value1" "$value2" "$diff"
    done
}

# Create new profile
create_profile() {
    local profile_name="$1"
    
    if [[ -z "$profile_name" ]]; then
        error "Profile name required"
        return 1
    fi
    
    if [[ "$profile_name" =~ [^a-zA-Z0-9_-] ]]; then
        error "Profile name can only contain letters, numbers, underscores, and hyphens"
        return 1
    fi
    
    local new_profile_file="$PROFILES_DIR/${profile_name}.conf"
    
    if [[ -f "$new_profile_file" ]]; then
        error "Profile already exists: $profile_name"
        return 1
    fi
    
    log "Creating new profile: $profile_name"
    
    # Use workstation profile as template
    local template_file="$PROFILES_DIR/workstation.conf"
    
    if [[ ! -f "$template_file" ]]; then
        error "Template profile not found: workstation.conf"
        return 1
    fi
    
    # Copy template and update profile name
    cp "$template_file" "$new_profile_file"
    sed -i.bak "s/PROFILE_NAME=\"workstation\"/PROFILE_NAME=\"$profile_name\"/" "$new_profile_file"
    sed -i.bak "s/PROFILE_DESCRIPTION=\".*\"/PROFILE_DESCRIPTION=\"Custom profile: $profile_name\"/" "$new_profile_file"
    rm -f "${new_profile_file}.bak"
    
    log "Profile created: $new_profile_file"
    info "Edit the file to customize settings for your system"
    info "Validate with: $0 validate $profile_name"
}

# Reset to default
reset_configuration() {
    log "Resetting to default configuration..."
    
    rm -f "$CONFIG_DIR/active.conf"
    rm -f "$CONFIG_DIR/current_profile.conf"
    
    log "Configuration reset complete"
    info "Pipeline will use default settings from config/default.conf"
}

# Main function
main() {
    local command="$1"
    shift
    
    case "$command" in
        list)
            list_profiles
            ;;
        show)
            show_profile "$1"
            ;;
        set)
            set_profile "$1"
            ;;
        create)
            create_profile "$1"
            ;;
        validate)
            validate_profile "$1"
            ;;
        compare)
            compare_profiles "$1" "$2"
            ;;
        reset)
            reset_configuration
            ;;
        help|--help|-h)
            show_help
            ;;
        *)
            if [[ -z "$command" ]]; then
                show_help
            else
                error "Unknown command: $command"
                echo "Use 'help' for usage information"
                exit 1
            fi
            ;;
    esac
}

# Run main function
main "$@"