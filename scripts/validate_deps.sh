#!/bin/bash

# WGS Pipeline Dependency Version Validator
# Enforces version ranges defined in DEPENDENCIES.md
# Usage: ./scripts/validate_deps.sh [OPTIONS]

set -e

# Script info
SCRIPT_VERSION="1.0.0"
SCRIPT_NAME="validate_deps.sh"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Options
VERBOSE=false
CI_MODE=false
SPECIFIC_TOOL=""
STRICT=false

# Version ranges (min_major.min_minor max_major.max_minor)
# Format: "tool|min_version|max_version|required"
# Note: max versions set high to allow minor/patch updates; review periodically
declare -a TOOL_VERSIONS=(
    "fastqc|0.11.9|0.12.99|required"
    "fastp|0.20.0|1.99.99|required"
    "bwa-mem2|2.2.1|2.99.99|required"
    "samtools|1.15|1.99.99|required"
    "bcftools|1.15|1.99.99|required"
    "vep|105.0|199.99|optional"
)

# Runtime dependencies
# Note: bash 3.2 is acceptable (macOS ships with 3.2 due to GPL3 licensing)
declare -a RUNTIME_DEPS=(
    "bash|3.2|5.99|required"
    "python3|3.8|3.99|required"
)

# Counters
PASS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0

# Help
show_help() {
    cat << EOF
$SCRIPT_NAME v$SCRIPT_VERSION - WGS Pipeline Dependency Validator

DESCRIPTION:
    Validates that all pipeline dependencies are within supported version ranges.
    Enforces the version policy defined in DEPENDENCIES.md.

USAGE:
    $0 [OPTIONS]

OPTIONS:
    -h, --help          Show this help message
    -v, --verbose       Show detailed version parsing
    --ci-mode           CI-friendly mode: still enforce runtime deps, but tolerate missing bio tools
    --strict            Fail on any warning (missing optional tools)
    --tool TOOL         Check only a specific tool
    --version           Show script version

EXAMPLES:
    # Full validation
    $0

    # CI mode (tolerant of missing optional tools)
    $0 --ci-mode

    # Check single tool
    $0 --tool samtools

EXIT CODES:
    0 - All required versions within range
    1 - Version out of range
    2 - Required tool missing
    3 - Invalid arguments

EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        --ci-mode)
            CI_MODE=true
            shift
            ;;
        --strict)
            STRICT=true
            shift
            ;;
        --tool)
            SPECIFIC_TOOL="$2"
            shift 2
            ;;
        --version)
            echo "$SCRIPT_NAME version $SCRIPT_VERSION"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 3
            ;;
    esac
done

# Logging
log_pass() { echo -e "${GREEN}✓ PASS${NC} $1"; PASS_COUNT=$((PASS_COUNT + 1)); }
log_fail() { echo -e "${RED}✗ FAIL${NC} $1"; FAIL_COUNT=$((FAIL_COUNT + 1)); }
log_warn() { echo -e "${YELLOW}⚠ WARN${NC} $1"; }
log_skip() { echo -e "${BLUE}○ SKIP${NC} $1"; SKIP_COUNT=$((SKIP_COUNT + 1)); }
log_info() { [[ "$VERBOSE" == "true" ]] && echo -e "${BLUE}ℹ INFO${NC} $1"; return 0; }

# Parse version string to comparable format
# Handles: "0.11.9", "1.15", "105", "v0.23.4"
parse_version() {
    local version="$1"
    # Remove leading 'v' if present
    version="${version#v}"
    # Remove any suffix after numeric portion (e.g., "-r123", "_beta")
    version=$(echo "$version" | sed -E 's/[-_].*//')
    # Extract major.minor.patch (fill missing with 0)
    local major minor patch
    IFS='.' read -r major minor patch <<< "$version"
    major=${major:-0}
    minor=${minor:-0}
    patch=${patch:-0}
    # Remove non-numeric characters
    major=$(echo "$major" | tr -cd '0-9')
    minor=$(echo "$minor" | tr -cd '0-9')
    patch=$(echo "$patch" | tr -cd '0-9')
    echo "$major $minor $patch"
}

# Compare versions: returns 0 if v1 >= v2
version_gte() {
    local v1_parts v2_parts
    read -r v1_maj v1_min v1_pat <<< "$(parse_version "$1")"
    read -r v2_maj v2_min v2_pat <<< "$(parse_version "$2")"
    
    log_info "Comparing $1 ($v1_maj.$v1_min.$v1_pat) >= $2 ($v2_maj.$v2_min.$v2_pat)"
    
    if [[ $v1_maj -gt $v2_maj ]]; then return 0; fi
    if [[ $v1_maj -lt $v2_maj ]]; then return 1; fi
    if [[ $v1_min -gt $v2_min ]]; then return 0; fi
    if [[ $v1_min -lt $v2_min ]]; then return 1; fi
    if [[ $v1_pat -ge $v2_pat ]]; then return 0; fi
    return 1
}

# Compare versions: returns 0 if v1 <= v2
version_lte() {
    local v1_parts v2_parts
    read -r v1_maj v1_min v1_pat <<< "$(parse_version "$1")"
    read -r v2_maj v2_min v2_pat <<< "$(parse_version "$2")"
    
    log_info "Comparing $1 ($v1_maj.$v1_min.$v1_pat) <= $2 ($v2_maj.$v2_min.$v2_pat)"
    
    if [[ $v1_maj -lt $v2_maj ]]; then return 0; fi
    if [[ $v1_maj -gt $v2_maj ]]; then return 1; fi
    if [[ $v1_min -lt $v2_min ]]; then return 0; fi
    if [[ $v1_min -gt $v2_min ]]; then return 1; fi
    if [[ $v1_pat -le $v2_pat ]]; then return 0; fi
    return 1
}

# Get tool version
get_tool_version() {
    local tool="$1"
    local version=""
    
    case "$tool" in
        fastqc)
            version=$(fastqc --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1)
            ;;
        fastp)
            version=$(fastp --version 2>&1 | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1)
            ;;
        bwa-mem2)
            version=$(bwa-mem2 version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1)
            ;;
        samtools)
            version=$(samtools --version 2>/dev/null | head -1 | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1)
            ;;
        bcftools)
            version=$(bcftools --version 2>/dev/null | head -1 | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1)
            ;;
        vep)
            version=$(vep --help 2>/dev/null | grep -i "version" | grep -oE '[0-9]+' | head -1)
            ;;
        bash)
            version="${BASH_VERSION%(*}"
            version=$(echo "$version" | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1)
            ;;
        python3)
            version=$(python3 --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1)
            ;;
        conda)
            version=$(conda --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1)
            ;;
    esac
    
    echo "$version"
}

# Validate a single tool
validate_tool() {
    local spec="$1"
    local tool_class="${2:-bio}"
    local tool min_ver max_ver required

    IFS='|' read -r tool min_ver max_ver required <<< "$spec"

    # Skip if checking specific tool and this isn't it
    if [[ -n "$SPECIFIC_TOOL" && "$tool" != "$SPECIFIC_TOOL" ]]; then
        return 0
    fi

    # Check if tool exists
    if ! command -v "$tool" &>/dev/null; then
        if [[ "$required" == "required" ]]; then
            if [[ "$CI_MODE" == "true" && "$tool_class" == "bio" ]]; then
                log_skip "$tool: not found (required bio tool, CI mode)"
                return 0
            fi
            log_fail "$tool: not found (required)"
            return 2
        else
            if [[ "$CI_MODE" == "true" ]]; then
                log_skip "$tool: not found (optional, CI mode)"
            else
                log_warn "$tool: not found (optional)"
                if [[ "$STRICT" == "true" ]]; then
                    FAIL_COUNT=$((FAIL_COUNT + 1))
                    return 1
                fi
            fi
            return 0
        fi
    fi
    
    # Get version
    local version
    version=$(get_tool_version "$tool")
    
    if [[ -z "$version" ]]; then
        log_warn "$tool: unable to parse version"
        if [[ "$STRICT" == "true" && "$required" == "required" ]]; then
            FAIL_COUNT=$((FAIL_COUNT + 1))
            return 1
        fi
        return 0
    fi
    
    log_info "$tool: detected version $version (range: $min_ver - $max_ver)"
    
    # Check minimum version
    if ! version_gte "$version" "$min_ver"; then
        log_fail "$tool: version $version < minimum $min_ver"
        echo "  └─ Upgrade with: conda install -c bioconda $tool>=$min_ver"
        return 1
    fi
    
    # Check maximum version
    if ! version_lte "$version" "$max_ver"; then
        log_fail "$tool: version $version > maximum $max_ver"
        echo "  └─ Downgrade or update DEPENDENCIES.md if version is known-good"
        return 1
    fi
    
    log_pass "$tool: $version (range: $min_ver - $max_ver)"
    return 0
}

# Main validation
main() {
    echo "WGS Pipeline Dependency Validator"
    echo "=================================="
    [[ "$CI_MODE" == "true" ]] && echo "Mode: CI (tolerant)"
    [[ "$STRICT" == "true" ]] && echo "Mode: Strict"
    [[ -n "$SPECIFIC_TOOL" ]] && echo "Checking: $SPECIFIC_TOOL only"
    echo ""
    
    local exit_code=0
    
    echo "Runtime Dependencies:"
    echo "---------------------"
    for spec in "${RUNTIME_DEPS[@]}"; do
        validate_tool "$spec" "runtime" || exit_code=$?
    done
    echo ""
    
    echo "Bioinformatics Tools:"
    echo "---------------------"
    for spec in "${TOOL_VERSIONS[@]}"; do
        validate_tool "$spec" "bio" || exit_code=$?
    done
    echo ""
    
    # Summary
    echo "Summary:"
    echo "--------"
    echo "  Passed: $PASS_COUNT"
    echo "  Failed: $FAIL_COUNT"
    echo "  Skipped: $SKIP_COUNT"
    echo ""
    
    if [[ $FAIL_COUNT -gt 0 ]]; then
        echo -e "${RED}❌ Dependency validation failed${NC}"
        echo "See DEPENDENCIES.md for version requirements."
        exit 1
    else
        echo -e "${GREEN}✅ All dependency versions validated${NC}"
        exit 0
    fi
}

main "$@"
