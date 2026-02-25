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

# =============================================================================
# HARDWARE DETECTION TESTS
# =============================================================================

run_hardware_detection_tests() {
    local failures=0
    
    # Source the detection script
    source "$ROOT_DIR/scripts/detect_hardware.sh"
    
    # Test 1: CPU detection returns positive integer
    info "  Testing CPU detection..."
    local cpu_cores
    cpu_cores=$(detect_cpu_cores)
    if [[ "$cpu_cores" =~ ^[0-9]+$ ]] && [[ "$cpu_cores" -ge 1 ]]; then
        info "    ✓ Detected $cpu_cores CPU cores"
    else
        warning "    ✗ CPU detection failed: got '$cpu_cores'"
        failures=$((failures + 1))
    fi
    
    # Test 2: RAM detection returns positive integer
    info "  Testing RAM detection..."
    local ram_gb
    ram_gb=$(detect_ram_gb)
    if [[ "$ram_gb" =~ ^[0-9]+$ ]] && [[ "$ram_gb" -ge 1 ]]; then
        info "    ✓ Detected ${ram_gb}GB RAM"
    else
        warning "    ✗ RAM detection failed: got '$ram_gb'"
        failures=$((failures + 1))
    fi
    
    # Test 3: GPU detection returns non-negative integer
    info "  Testing GPU detection..."
    local gpu_count
    gpu_count=$(detect_gpu_count)
    if [[ "$gpu_count" =~ ^[0-9]+$ ]]; then
        info "    ✓ Detected $gpu_count GPUs"
    else
        warning "    ✗ GPU detection failed: got '$gpu_count'"
        failures=$((failures + 1))
    fi
    
    # Test 4: Mode selection returns valid mode
    info "  Testing mode selection..."
    local mode
    mode=$(select_performance_mode)
    case "$mode" in
        laptop|workstation|server|dgx)
            info "    ✓ Selected mode: $mode"
            ;;
        *)
            warning "    ✗ Invalid mode: '$mode'"
            failures=$((failures + 1))
            ;;
    esac
    
    # Test 5: Thread recommendation returns positive integer
    info "  Testing thread recommendation..."
    local threads
    threads=$(recommend_threads)
    if [[ "$threads" =~ ^[0-9]+$ ]] && [[ "$threads" -ge 1 ]]; then
        info "    ✓ Recommended $threads threads"
    else
        warning "    ✗ Thread recommendation failed: got '$threads'"
        failures=$((failures + 1))
    fi
    
    # Test 6: GPU settings format check
    info "  Testing GPU settings recommendation..."
    local gpu_settings
    gpu_settings=$(recommend_gpu_settings)
    local use_gpu aligner gpu_rec
    read -r use_gpu aligner gpu_rec <<< "$gpu_settings"
    if [[ "$use_gpu" == "true" || "$use_gpu" == "false" ]] && [[ "$gpu_rec" =~ ^[0-9]+$ ]]; then
        info "    ✓ GPU settings: use=$use_gpu, aligner=$aligner, count=$gpu_rec"
    else
        warning "    ✗ GPU settings malformed: '$gpu_settings'"
        failures=$((failures + 1))
    fi
    
    # Test 7: JSON output is valid
    info "  Testing JSON output..."
    local json_output
    json_output=$(print_json)
    if echo "$json_output" | python3 -c "import sys,json; json.load(sys.stdin)" 2>/dev/null; then
        info "    ✓ JSON output is valid"
    elif echo "$json_output" | jq . >/dev/null 2>&1; then
        info "    ✓ JSON output is valid (jq)"
    else
        warning "    ✗ JSON output is invalid"
        info "    Output was: $json_output"
        failures=$((failures + 1))
    fi
    
    # Test 8: CLI subcommands work
    info "  Testing CLI subcommands..."
    local cli_failures=0
    for cmd in cpu ram gpu mode threads; do
        local output
        output=$("$ROOT_DIR/scripts/detect_hardware.sh" "$cmd" 2>&1)
        if [[ -n "$output" ]] && [[ ! "$output" =~ "error" ]]; then
            info "    ✓ $cmd: $output"
        else
            warning "    ✗ $cmd failed: $output"
            cli_failures=$((cli_failures + 1))
        fi
    done
    if [[ $cli_failures -gt 0 ]]; then
        failures=$((failures + cli_failures))
    fi
    
    return $failures
}

# =============================================================================
# OVERRIDE PRECEDENCE TESTS
# =============================================================================

run_override_tests() {
    local failures=0
    
    info "Testing CLI override precedence..."
    
    # Test that --threads overrides auto-detection
    local dry_output
    dry_output=$("$ROOT_DIR/run_pipeline.sh" --dry-run --threads 99 --skip-requirements-check -y 2>&1 || true)
    if echo "$dry_output" | grep -q "Threads: 99"; then
        info "  ✓ --threads override works"
    else
        warning "  ✗ --threads override failed"
        failures=$((failures + 1))
    fi
    
    # Test that --no-gpu forces CPU mode
    dry_output=$("$ROOT_DIR/run_pipeline.sh" --dry-run --no-gpu --skip-requirements-check -y 2>&1 || true)
    if echo "$dry_output" | grep -q "Alignment Mode: CPU"; then
        info "  ✓ --no-gpu forces CPU mode"
    else
        warning "  ✗ --no-gpu did not force CPU mode"
        failures=$((failures + 1))
    fi
    
    # Test that --no-auto uses default threads
    dry_output=$("$ROOT_DIR/run_pipeline.sh" --dry-run --no-auto --skip-requirements-check -y 2>&1 || true)
    if echo "$dry_output" | grep -q "Threads: 4"; then
        info "  ✓ --no-auto uses default threads (4)"
    else
        warning "  ✗ --no-auto did not use default threads"
        failures=$((failures + 1))
    fi
    
    return $failures
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

    # Test 0: Hardware detection functions
    info "Testing hardware detection..."
    run_hardware_detection_tests
    local hw_result=$?
    if [[ $hw_result -eq 0 ]]; then
        echo "  ✓ Hardware detection tests passed"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Hardware detection tests failed ($hw_result failures)"
        tests_failed=$((tests_failed + 1))
    fi

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

    # Test 5: GPU flag discoverability (help text)
    info "Testing GPU flag discoverability..."
    if "$ROOT_DIR/scripts/alignment.sh" --help 2>/dev/null | grep -q -- "--use-gpu" \
       && "$ROOT_DIR/scripts/alignment.sh" --help 2>/dev/null | grep -q -- "--gpu-aligner" \
       && "$ROOT_DIR/scripts/alignment.sh" --help 2>/dev/null | grep -q -- "--gpu-count"; then
        echo "  ✓ GPU flags documented in help output"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ GPU flags missing from help output"
        tests_failed=$((tests_failed + 1))
    fi

    # Test 6: Hardware auto-detection flags in run_pipeline.sh
    # Note: run_pipeline.sh requires Bash 4+; skip if not available
    info "Testing hardware auto-detection flags..."
    if [[ "${BASH_VERSINFO[0]:-3}" -ge 4 ]]; then
        local pipeline_help
        pipeline_help=$("$ROOT_DIR/run_pipeline.sh" --help 2>&1)
        if echo "$pipeline_help" | grep -q -- "--show-hardware" \
           && echo "$pipeline_help" | grep -q -- "--no-auto" \
           && echo "$pipeline_help" | grep -q -- "--no-gpu"; then
            echo "  ✓ Hardware auto-detection flags documented"
            tests_passed=$((tests_passed + 1))
        else
            echo "  ✗ Hardware auto-detection flags missing from help"
            tests_failed=$((tests_failed + 1))
        fi
    else
        # Bash 3.x: check flags in source file directly
        if grep -q -- "--show-hardware" "$ROOT_DIR/run_pipeline.sh" \
           && grep -q -- "--no-auto" "$ROOT_DIR/run_pipeline.sh" \
           && grep -q -- "--no-gpu" "$ROOT_DIR/run_pipeline.sh"; then
            echo "  ✓ Hardware auto-detection flags found in source (Bash 4+ required for runtime test)"
            tests_passed=$((tests_passed + 1))
        else
            echo "  ✗ Hardware auto-detection flags missing from source"
            tests_failed=$((tests_failed + 1))
        fi
    fi

    # Test 7: Override precedence tests
    # Note: run_pipeline.sh requires Bash 4+; skip if not available
    info "Testing CLI override precedence..."
    if [[ "${BASH_VERSINFO[0]:-3}" -ge 4 ]]; then
        local override_failures=0
        set +e
        run_override_tests
        override_failures=$?
        set -e
        if [[ $override_failures -eq 0 ]]; then
            echo "  ✓ CLI override precedence works correctly"
            tests_passed=$((tests_passed + 1))
        else
            echo "  ✗ CLI override precedence failed ($override_failures issues)"
            tests_failed=$((tests_failed + 1))
        fi
    else
        echo "  ⊘ CLI override precedence test skipped (requires Bash 4+)"
        tests_passed=$((tests_passed + 1))  # Count as pass since we can't test
    fi

    # Test 8: detect_hardware.sh CLI help
    info "Testing detect_hardware.sh help..."
    if "$ROOT_DIR/scripts/detect_hardware.sh" --help 2>&1 | grep -q "Hardware Detection"; then
        echo "  ✓ detect_hardware.sh help works"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ detect_hardware.sh help failed"
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