#!/bin/bash

# WGS Analysis Requirements Checker
# This script validates system requirements before running WGS analysis
# Usage: ./scripts/check_requirements.sh [OPTIONS]

set -e  # Exit on any error

# Script version and info
SCRIPT_VERSION="1.0"
SCRIPT_NAME="check_requirements.sh"

# Default configuration (16GB optimized)
CONDA_ENV="wgs_analysis"
MIN_RAM_GB=16
MIN_DISK_GB=400
REQUIRED_TOOLS=("fastqc" "fastp" "bwa" "samtools" "bcftools" "vep")
VERBOSE=false
SKIP_CONDA=false
SKIP_ENV=false
SKIP_TOOLS=false

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Help function
show_help() {
    cat << EOF
$SCRIPT_NAME v$SCRIPT_VERSION - WGS Analysis Requirements Checker

DESCRIPTION:
    Validates system requirements for 16GB RAM optimized WGS analysis pipeline:
    - 16GB RAM minimum with efficient memory usage
    - 400GB disk space for complete analysis
    - Required bioinformatics tools
    - Local system optimization checks

USAGE:
    $0 [OPTIONS]

OPTIONS:
    -h, --help              Show this help message and exit
    -e, --env ENV           Conda environment to check (default: $CONDA_ENV)
    --min-ram GB            Minimum RAM required in GB (default: $MIN_RAM_GB)
    --min-disk GB           Minimum disk space in GB (default: $MIN_DISK_GB)
    --skip-conda            Skip conda installation check (smoke testing only)
    --skip-env              Skip conda environment check (smoke testing only)
    --skip-tools            Skip bioinformatics tool checks (smoke testing only)
    --verbose               Enable verbose output
    --version               Show version information

EXAMPLES:
    # Basic requirements check
    $0

    # Check specific environment with verbose output
    $0 --env my_wgs_env --verbose

    # Check with custom resource requirements
    $0 --min-ram 32 --min-disk 1000

    # CI/local smoke check without conda/toolchain dependencies
    $0 --min-ram 1 --min-disk 1 --skip-conda --skip-env --skip-tools

REQUIREMENTS CHECKED:
    - System RAM (default: >= ${MIN_RAM_GB}GB)
    - Available disk space (default: >= ${MIN_DISK_GB}GB)
    - Conda environment existence and activation
    - Required bioinformatics tools installation
    - File system permissions

OUTPUT:
    - Color-coded status for each requirement
    - Resource usage estimates for WGS analysis
    - Recommendations for system optimization
    - Installation commands for missing tools

AUTHOR:
    WGS Analysis Pipeline
    https://github.com/Hilo-Hilo/wgs-analysis-pipeline

EOF
}

# Version function
show_version() {
    echo "$SCRIPT_NAME version $SCRIPT_VERSION"
    exit 0
}

# Parse command line arguments
parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                show_help
                exit 0
                ;;
            --version)
                show_version
                ;;
            -e|--env)
                CONDA_ENV="$2"
                shift 2
                ;;
            --min-ram)
                MIN_RAM_GB="$2"
                shift 2
                ;;
            --min-disk)
                MIN_DISK_GB="$2"
                shift 2
                ;;
            --skip-conda)
                SKIP_CONDA=true
                shift
                ;;
            --skip-env)
                SKIP_ENV=true
                shift
                ;;
            --skip-tools)
                SKIP_TOOLS=true
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
print_status() {
    local status="$1"
    local message="$2"
    case "$status" in
        "PASS")
            echo -e "${GREEN}✓ PASS${NC} $message"
            ;;
        "FAIL")
            echo -e "${RED}✗ FAIL${NC} $message"
            ;;
        "WARN")
            echo -e "${YELLOW}⚠ WARN${NC} $message"
            ;;
        "INFO")
            if [[ "$VERBOSE" == "true" ]]; then
                echo -e "${BLUE}ℹ INFO${NC} $message"
            fi
            ;;
    esac
}

# Convert bytes to human readable format
bytes_to_human() {
    local bytes=$1
    if (( bytes >= 1073741824 )); then
        echo "$(( bytes / 1073741824 )) GB"
    elif (( bytes >= 1048576 )); then
        echo "$(( bytes / 1048576 )) MB"
    else
        echo "$bytes bytes"
    fi
}

# Resolve sysctl binary for macOS in non-interactive shells
resolve_sysctl_bin() {
    local sysctl_bin
    sysctl_bin=$(command -v sysctl 2>/dev/null || true)
    if [[ -z "$sysctl_bin" && -x "/usr/sbin/sysctl" ]]; then
        sysctl_bin="/usr/sbin/sysctl"
    fi
    echo "$sysctl_bin"
}

# Check system RAM
check_system_ram() {
    print_status "INFO" "Checking system RAM..."
    
    local total_ram_kb
    if [[ "$(uname)" == "Darwin" ]]; then
        # macOS
        local sysctl_bin
        local total_ram_bytes
        sysctl_bin=$(resolve_sysctl_bin)
        if [[ -z "$sysctl_bin" ]] || ! total_ram_bytes=$($sysctl_bin -n hw.memsize 2>/dev/null); then
            print_status "FAIL" "Unable to determine system RAM on macOS (sysctl not available)"
            return 1
        fi
        total_ram_kb=$((total_ram_bytes / 1024))
    else
        # Linux
        total_ram_kb=$(grep MemTotal /proc/meminfo | awk '{print $2}')
    fi
    
    local total_ram_gb=$((total_ram_kb / 1024 / 1024))
    
    if [[ $total_ram_gb -ge $MIN_RAM_GB ]]; then
        print_status "PASS" "System RAM: ${total_ram_gb}GB (>= ${MIN_RAM_GB}GB required)"
    else
        print_status "FAIL" "System RAM: ${total_ram_gb}GB (< ${MIN_RAM_GB}GB required)"
        echo "  Recommendation: WGS analysis requires at least ${MIN_RAM_GB}GB RAM"
        echo "  For optimal performance, consider 32GB+ RAM"
        return 1
    fi
    
    # Check available RAM
    local available_ram_kb
    if [[ "$(uname)" == "Darwin" ]]; then
        local vm_stat_bin
        vm_stat_bin=$(command -v vm_stat 2>/dev/null || true)
        if [[ -z "$vm_stat_bin" && -x "/usr/bin/vm_stat" ]]; then
            vm_stat_bin="/usr/bin/vm_stat"
        fi

        if [[ -n "$vm_stat_bin" ]]; then
            available_ram_kb=$($vm_stat_bin | awk '/free/ {print $3}' | sed 's/\.//')
            available_ram_kb=$((available_ram_kb * 4))  # Convert pages to KB
        else
            available_ram_kb=0
            print_status "WARN" "Unable to determine available RAM on macOS (vm_stat not available)"
        fi
    else
        available_ram_kb=$(grep MemAvailable /proc/meminfo | awk '{print $2}')
    fi
    
    local available_ram_gb=$((available_ram_kb / 1024 / 1024))
    print_status "INFO" "Available RAM: ${available_ram_gb}GB"
    
    return 0
}

# Check disk space
check_disk_space() {
    print_status "INFO" "Checking disk space..."
    
    local current_dir="."
    local available_kb
    available_kb=$(df "$current_dir" | awk 'NR==2 {print $4}')
    local available_gb=$((available_kb / 1024 / 1024))
    
    if [[ $available_gb -ge $MIN_DISK_GB ]]; then
        print_status "PASS" "Available disk space: ${available_gb}GB (>= ${MIN_DISK_GB}GB required)"
    else
        print_status "FAIL" "Available disk space: ${available_gb}GB (< ${MIN_DISK_GB}GB required)"
        echo "  Recommendation: WGS analysis requires at least ${MIN_DISK_GB}GB free space"
        echo "  Typical space usage:"
        echo "    - Raw FASTQ files: 50-100GB"
        echo "    - Reference genomes: 3-10GB"
        echo "    - Alignment files: 100-200GB"
        echo "    - Variant files: 1-10GB"
        return 1
    fi
    
    return 0
}

# Check conda installation
check_conda() {
    print_status "INFO" "Checking conda installation..."
    
    if ! command -v conda &> /dev/null; then
        print_status "FAIL" "Conda not found in PATH"
        echo "  Install conda from: https://docs.conda.io/en/latest/miniconda.html"
        return 1
    fi
    
    local conda_version
    conda_version=$(conda --version 2>/dev/null | awk '{print $2}')
    print_status "PASS" "Conda installed: version $conda_version"
    
    return 0
}

# Check conda environment
check_conda_environment() {
    print_status "INFO" "Checking conda environment: $CONDA_ENV"

    if ! command -v conda &> /dev/null; then
        print_status "FAIL" "Conda not found in PATH (required to check environment '$CONDA_ENV')"
        return 1
    fi
    
    if ! conda env list | grep -q "^$CONDA_ENV\s"; then
        print_status "FAIL" "Conda environment '$CONDA_ENV' not found"
        echo "  Create environment with:"
        echo "    conda create -n $CONDA_ENV -c bioconda -c conda-forge \\"
        echo "      python=3.9 fastqc fastp bwa samtools bcftools"
        return 1
    fi
    
    print_status "PASS" "Conda environment '$CONDA_ENV' exists"
    
    # Check if environment is currently active
    if [[ "$CONDA_DEFAULT_ENV" == "$CONDA_ENV" ]]; then
        print_status "PASS" "Environment '$CONDA_ENV' is currently active"
    else
        print_status "WARN" "Environment '$CONDA_ENV' is not active"
        echo "  Activate with: conda activate $CONDA_ENV"
    fi
    
    return 0
}

# Check required tools
check_tools() {
    print_status "INFO" "Checking required bioinformatics tools..."
    
    local missing_tools=()
    local failed=false
    
    for tool in "${REQUIRED_TOOLS[@]}"; do
        if command -v "$tool" &> /dev/null; then
            local version
            case "$tool" in
                "fastqc")
                    version=$(fastqc --version 2>/dev/null | awk '{print $2}')
                    ;;
                "fastp")
                    version=$(fastp --version 2>&1 | head -1 | awk '{print $2}')
                    ;;
                "bwa")
                    version=$(bwa 2>&1 | grep "Version" | awk '{print $2}')
                    ;;
                "samtools")
                    version=$(samtools --version 2>/dev/null | head -1 | awk '{print $2}')
                    ;;
                "bcftools")
                    version=$(bcftools --version 2>/dev/null | head -1 | awk '{print $2}')
                    ;;
                *)
                    version="unknown"
                    ;;
            esac
            print_status "PASS" "$tool installed (version: $version)"
        else
            print_status "FAIL" "$tool not found"
            missing_tools+=("$tool")
            failed=true
        fi
    done
    
    if [[ "$failed" == "true" ]]; then
        echo "  Install missing tools with:"
        echo "    conda install -c bioconda ${missing_tools[*]}"
        return 1
    fi
    
    return 0
}

# Check file permissions
check_file_permissions() {
    print_status "INFO" "Checking file system permissions..."
    
    # Check write permissions in current directory
    if [[ -w "." ]]; then
        print_status "PASS" "Current directory is writable"
    else
        print_status "FAIL" "Current directory is not writable"
        echo "  Ensure you have write permissions in the working directory"
        return 1
    fi
    
    # Test creating and removing a temporary file
    local temp_file=".wgs_requirements_test_$$"
    if touch "$temp_file" 2>/dev/null && rm "$temp_file" 2>/dev/null; then
        print_status "PASS" "File creation/deletion works"
    else
        print_status "FAIL" "Cannot create/delete files"
        return 1
    fi
    
    return 0
}

# Estimate resource usage
estimate_resource_usage() {
    print_status "INFO" "Resource usage estimates for WGS analysis..."
    
    echo ""
    echo "Expected resource usage for 30x WGS analysis:"
    echo "┌─────────────────────┬─────────────┬─────────────┐"
    echo "│ Analysis Step       │ Time        │ Peak Memory │"
    echo "├─────────────────────┼─────────────┼─────────────┤"
    echo "│ Quality Control     │ 10-30 min   │ 1-2 GB      │"
    echo "│ Read Cleaning       │ 30-60 min   │ 2-4 GB      │"
    echo "│ Alignment (BWA)     │ 2-6 hours   │ 8-16 GB     │"
    echo "│ Variant Calling     │ 1-4 hours   │ 4-8 GB      │"
    echo "│ Annotation (VEP)    │ 30-90 min   │ 8-16 GB     │"
    echo "└─────────────────────┴─────────────┴─────────────┘"
    echo ""
    echo "Storage requirements:"
    echo "  - Raw FASTQ files:     50-100 GB"
    echo "  - Cleaned FASTQ files: 40-80 GB"
    echo "  - Reference genomes:   3-10 GB"
    echo "  - Alignment (BAM):     80-150 GB"
    echo "  - Variants (VCF):      1-10 GB"
    echo "  - Annotations:         5-20 GB"
    echo "  ─────────────────────────────────"
    echo "  Total estimated:       180-370 GB"
    echo ""
}

# Generate recommendations
generate_recommendations() {
    echo ""
    echo "System Optimization Recommendations:"
    echo "────────────────────────────────────"
    
    # RAM recommendations
    local total_ram_kb
    if [[ "$(uname)" == "Darwin" ]]; then
        local sysctl_bin
        local total_ram_bytes
        sysctl_bin=$(resolve_sysctl_bin)
        if [[ -n "$sysctl_bin" ]] && total_ram_bytes=$($sysctl_bin -n hw.memsize 2>/dev/null); then
            total_ram_kb=$((total_ram_bytes / 1024))
        else
            total_ram_kb=0
            print_status "WARN" "Unable to determine total RAM for optimization hints"
        fi
    else
        total_ram_kb=$(grep MemTotal /proc/meminfo | awk '{print $2}')
    fi
    local total_ram_gb=$((total_ram_kb / 1024 / 1024))
    
    if [[ $total_ram_gb -lt 32 ]]; then
        echo "• Consider upgrading to 32GB+ RAM for optimal performance"
    fi
    
    # CPU recommendations
    local cpu_cores
    if [[ "$(uname)" == "Darwin" ]]; then
        local sysctl_bin
        sysctl_bin=$(resolve_sysctl_bin)
        if [[ -n "$sysctl_bin" ]]; then
            cpu_cores=$($sysctl_bin -n hw.ncpu 2>/dev/null || echo "unknown")
        else
            cpu_cores="unknown"
        fi
    else
        cpu_cores=$(nproc)
    fi
    
    echo "• Use $cpu_cores CPU cores for parallel processing (detected on system)"
    
    # Storage recommendations  
    echo "• Use SSD storage for faster I/O operations"
    echo "• Monitor disk space during analysis - clean intermediate files if needed"
    
    # Environment recommendations
    echo "• Run analysis in screen/tmux for long-running jobs"
    echo "• Use cloud instances for large-scale analysis if local resources insufficient"
    
    echo ""
}

# Main execution
main() {
    parse_arguments "$@"
    
    echo "WGS Analysis Requirements Check"
    echo "==============================="
    echo "Environment: $CONDA_ENV"
    echo "Minimum RAM: ${MIN_RAM_GB}GB"
    echo "Minimum Disk: ${MIN_DISK_GB}GB"
    echo ""
    
    local failed_checks=0
    
    # Run all checks
    check_system_ram || failed_checks=$((failed_checks + 1))
    echo ""
    
    check_disk_space || failed_checks=$((failed_checks + 1))
    echo ""
    
    if [[ "$SKIP_CONDA" == "true" ]]; then
        print_status "WARN" "Skipping conda installation check (--skip-conda)"
    else
        check_conda || failed_checks=$((failed_checks + 1))
    fi
    echo ""
    
    if [[ "$SKIP_ENV" == "true" ]]; then
        print_status "WARN" "Skipping conda environment check (--skip-env)"
    else
        check_conda_environment || failed_checks=$((failed_checks + 1))
    fi
    echo ""
    
    if [[ "$SKIP_TOOLS" == "true" ]]; then
        print_status "WARN" "Skipping bioinformatics tool checks (--skip-tools)"
    else
        check_tools || failed_checks=$((failed_checks + 1))
    fi
    echo ""
    
    check_file_permissions || failed_checks=$((failed_checks + 1))
    echo ""
    
    # Show resource estimates
    estimate_resource_usage
    
    # Generate recommendations
    generate_recommendations
    
    # Final summary
    echo "Summary:"
    echo "──────────"
    if [[ $failed_checks -eq 0 ]]; then
        print_status "PASS" "All requirements checks passed!"
        echo ""
        echo "Your system is ready for WGS analysis."
        echo "You can proceed with the analysis pipeline."
    else
        print_status "FAIL" "$failed_checks requirement(s) failed"
        echo ""
        echo "Please address the failed requirements before running WGS analysis."
        echo "Refer to the recommendations above for guidance."
        exit 1
    fi
}

# Run main function
main "$@"