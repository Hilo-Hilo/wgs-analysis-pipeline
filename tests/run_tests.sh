#!/bin/bash

# WGS Pipeline Test Suite
# Comprehensive testing framework for the 16GB-optimized WGS pipeline

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Test configuration
TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$TEST_DIR")"
SAMPLE_DATA_DIR="$TEST_DIR/sample_data"
EXPECTED_DIR="$TEST_DIR/expected_outputs"
TEMP_TEST_DIR="/tmp/wgs_pipeline_test_$$"
VERBOSE=false
CLEANUP=true
UNIT_ONLY=false
INTEGRATION_ONLY=false
QUICK_TESTS=false

# Logging functions
log() {
    echo -e "${GREEN}[TEST]${NC} $1"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

warning() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

info() {
    if [[ "$VERBOSE" == "true" ]]; then
        echo -e "${BLUE}[INFO]${NC} $1"
    fi
    return 0
}

# Help function
show_help() {
    cat << EOF
WGS Pipeline Test Suite

USAGE:
    $0 [OPTIONS]

OPTIONS:
    -h, --help              Show this help message
    -v, --verbose           Enable verbose output
    --no-cleanup            Don't clean up test files after completion
    --unit-only             Run only unit tests
    --integration-only      Run only integration tests
    --quick                 Run quick tests only (skip long-running tests)

EXAMPLES:
    # Run all tests
    $0

    # Run with verbose output
    $0 --verbose

    # Run only quick tests
    $0 --quick --verbose
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
        --no-cleanup)
            CLEANUP=false
            shift
            ;;
        --unit-only)
            UNIT_ONLY=true
            shift
            ;;
        --integration-only)
            INTEGRATION_ONLY=true
            shift
            ;;
        --quick)
            QUICK_TESTS=true
            shift
            ;;
        *)
            error "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Setup test environment
setup_test_environment() {
    log "Setting up test environment..."
    
    # Create temporary test directory
    mkdir -p "$TEMP_TEST_DIR"
    cd "$TEMP_TEST_DIR"
    
    # Create test directory structure
    mkdir -p {data/raw,data/processed,results,logs,config}
    
    # Copy configuration
    cp "$ROOT_DIR/config/default.conf" config/
    
    info "Test environment created at: $TEMP_TEST_DIR"
}

# Generate minimal test FASTQ data
generate_test_data() {
    log "Generating test FASTQ data..."
    
    local output_r1="$TEMP_TEST_DIR/data/raw/test_sample_R1.fastq.gz"
    local output_r2="$TEMP_TEST_DIR/data/raw/test_sample_R2.fastq.gz"
    
    # Generate 1000 paired-end reads (tiny dataset for testing)
    python3 - << EOF
import gzip
import random

# Simple FASTQ read generator
def generate_read(read_id, length=100):
    bases = 'ATCG'
    sequence = ''.join(random.choice(bases) for _ in range(length))
    quality = ''.join(chr(random.randint(33, 73)) for _ in range(length))
    return f"@read_{read_id}\n{sequence}\n+\n{quality}\n"

# Generate R1 file
with gzip.open('$output_r1', 'wt') as f1:
    for i in range(1000):
        f1.write(generate_read(f"{i}_1"))

# Generate R2 file  
with gzip.open('$output_r2', 'wt') as f2:
    for i in range(1000):
        f2.write(generate_read(f"{i}_2"))
EOF
    
    info "Generated test FASTQ files: $(du -h $output_r1 $output_r2)"
}

# Unit tests
run_unit_tests() {
    log "Running unit tests..."

    local tests_passed=0
    local tests_failed=0

    # Test 1: Requirements checker smoke mode (portable, no conda/toolchain required)
    info "Testing requirements checker (smoke mode)..."
    if "$ROOT_DIR/scripts/check_requirements.sh" \
        --min-ram 1 \
        --min-disk 1 \
        --skip-conda \
        --skip-env \
        --skip-tools > /dev/null 2>&1; then
        echo "  ✓ Requirements checker smoke mode works"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Requirements checker smoke mode failed"
        tests_failed=$((tests_failed + 1))
    fi

    # Test 2: Configuration loading
    info "Testing configuration loading..."
    if bash -c "cd '$TEMP_TEST_DIR' && source '$ROOT_DIR/scripts/load_config.sh' && load_config config/default.conf" > /dev/null 2>&1; then
        echo "  ✓ Configuration loading works"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Configuration loading failed"
        tests_failed=$((tests_failed + 1))
    fi

    # Test 3: Help messages for scripts
    info "Testing help messages..."
    local help_tests_passed=0
    for script in "$ROOT_DIR/scripts"/*.sh; do
        local script_name
        script_name=$(basename "$script")

        if "$script" --help > /dev/null 2>&1; then
            info "  ✓ $script_name help works"
            help_tests_passed=$((help_tests_passed + 1))
        else
            warning "  ✗ $script_name help failed"
        fi
    done

    if [[ $help_tests_passed -gt 0 ]]; then
        echo "  ✓ Help messages work ($help_tests_passed scripts)"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ No help messages work"
        tests_failed=$((tests_failed + 1))
    fi

    # Test 4: Argument handling (invalid flag should be rejected)
    info "Testing script argument handling..."
    local arg_scripts=(
        "$ROOT_DIR/scripts/check_requirements.sh"
        "$ROOT_DIR/scripts/quality_control.sh"
        "$ROOT_DIR/scripts/data_cleaning.sh"
        "$ROOT_DIR/scripts/alignment.sh"
        "$ROOT_DIR/scripts/variant_calling.sh"
        "$ROOT_DIR/scripts/vep_annotation.sh"
    )
    local arg_tests_failed=0

    for script in "${arg_scripts[@]}"; do
        local script_name
        script_name=$(basename "$script")

        if "$script" --definitely-invalid-option > /dev/null 2>&1; then
            warning "  ✗ $script_name accepted an invalid option"
            arg_tests_failed=$((arg_tests_failed + 1))
        else
            info "  ✓ $script_name rejected invalid option"
        fi
    done

    if [[ $arg_tests_failed -eq 0 ]]; then
        echo "  ✓ Script argument handling works"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Script argument handling failed ($arg_tests_failed script(s))"
        tests_failed=$((tests_failed + 1))
    fi

    log "Unit tests completed: $tests_passed passed, $tests_failed failed"
    return $tests_failed
}

# Integration tests
run_integration_tests() {
    if [[ "$QUICK_TESTS" == "true" ]]; then
        log "Skipping integration tests (quick mode)"
        return 0
    fi
    
    log "Running integration tests..."
    
    local tests_passed=0
    local tests_failed=0
    
    # Test 1: End-to-end pipeline with test data
    info "Testing quality control step..."
    if "$ROOT_DIR/scripts/quality_control.sh" \
        --input-dir data/raw \
        --output-dir results/quality_control \
        --threads 1 > logs/integration_qc.log 2>&1; then
        echo "  ✓ Quality control completed"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Quality control failed"
        warning "    See logs/integration_qc.log"
        tests_failed=$((tests_failed + 1))
    fi

    # Test 2: Data cleaning
    info "Testing data cleaning step..."
    if timeout 120 "$ROOT_DIR/scripts/data_cleaning.sh" \
        --input-dir data/raw \
        --output-dir data/processed \
        --sample-id test_sample \
        --threads 1 \
        --min-length 50 > logs/integration_clean.log 2>&1; then
        echo "  ✓ Data cleaning completed"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Data cleaning failed or timed out"
        warning "    Log output:"
        cat logs/integration_clean.log 2>/dev/null | tail -30 || true
        tests_failed=$((tests_failed + 1))
    fi

    # Test 3: Check output files exist
    info "Testing output file generation..."
    local output_files=(
        "results/quality_control/test_sample_R1_fastqc.html"
        "results/quality_control/test_sample_R2_fastqc.html"
        "data/processed/test_sample_clean_R1.fq.gz"
        "data/processed/test_sample_clean_R2.fq.gz"
    )
    
    local files_found=0
    for file in "${output_files[@]}"; do
        if [[ -f "$file" ]]; then
            info "  ✓ Found: $file"
            files_found=$((files_found + 1))
        else
            warning "  ✗ Missing: $file"
        fi
    done
    
    if [[ $files_found -eq ${#output_files[@]} ]]; then
        echo "  ✓ All expected output files generated"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Some output files missing ($files_found/${#output_files[@]})"
        tests_failed=$((tests_failed + 1))
    fi
    
    log "Integration tests completed: $tests_passed passed, $tests_failed failed"
    return $tests_failed
}

# Cleanup function
cleanup_test_environment() {
    if [[ "$CLEANUP" == "true" ]]; then
        log "Cleaning up test environment..."
        cd /tmp
        rm -rf "$TEMP_TEST_DIR"
        info "Test directory removed: $TEMP_TEST_DIR"
    else
        log "Test environment preserved at: $TEMP_TEST_DIR"
    fi
}

# Main test execution
main() {
    log "Starting WGS Pipeline Test Suite"
    log "================================"
    
    local total_failures=0
    
    # Setup
    setup_test_environment

    # Integration data is only required when integration tests will actually run
    if [[ "$UNIT_ONLY" != "true" && "$QUICK_TESTS" != "true" ]]; then
        generate_test_data
    fi
    
    # Run tests
    if [[ "$INTEGRATION_ONLY" != "true" ]]; then
        local unit_failures=0
        set +e
        run_unit_tests
        unit_failures=$?
        set -e
        total_failures=$((total_failures + unit_failures))
    fi
    
    if [[ "$UNIT_ONLY" != "true" ]]; then
        local integration_failures=0
        set +e
        run_integration_tests
        integration_failures=$?
        set -e
        total_failures=$((total_failures + integration_failures))
    fi
    
    # Results
    echo
    if [[ $total_failures -eq 0 ]]; then
        log "🎉 All tests passed!"
        echo "The WGS pipeline is ready for use."
    else
        error "❌ $total_failures test(s) failed"
        echo "Please check the pipeline configuration and dependencies."
    fi
    
    # Cleanup
    cleanup_test_environment
    
    exit $total_failures
}

# Handle interrupts
trap 'error "Tests interrupted"; cleanup_test_environment; exit 1' INT TERM

# Run tests
main "$@"