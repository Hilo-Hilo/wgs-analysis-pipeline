#!/bin/bash
#
# Unit tests for smoke test harness components
#
# These tests validate the data generators and harness logic
# without running the full pipeline. Suitable for CI environments.
#
# Usage:
#   ./tests/smoke/test_smoke_components.sh
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
TEST_TMP="/tmp/smoke-component-test-$$"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

TESTS_PASSED=0
TESTS_FAILED=0

setup() {
    mkdir -p "$TEST_TMP"
}

teardown() {
    rm -rf "$TEST_TMP"
}

pass() {
    echo -e "  ${GREEN}✓${NC} $1"
    TESTS_PASSED=$((TESTS_PASSED + 1))
}

fail() {
    echo -e "  ${RED}✗${NC} $1"
    TESTS_FAILED=$((TESTS_FAILED + 1))
}

warn() {
    echo -e "  ${YELLOW}!${NC} $1"
}

# =============================================================================
# Test Cases
# =============================================================================

test_reference_generator_exists() {
    echo "Testing: reference generator script exists and is executable"
    
    if [[ -x "$SCRIPT_DIR/generate_reference.py" ]]; then
        pass "generate_reference.py exists and is executable"
    elif [[ -f "$SCRIPT_DIR/generate_reference.py" ]]; then
        warn "generate_reference.py exists but is not executable"
        pass "generate_reference.py exists"
    else
        fail "generate_reference.py not found"
    fi
}

test_reference_generator_help() {
    echo "Testing: reference generator --help"
    
    if python3 "$SCRIPT_DIR/generate_reference.py" --help > /dev/null 2>&1; then
        pass "--help works"
    else
        fail "--help failed"
    fi
}

test_reference_generator_output() {
    echo "Testing: reference generator produces valid output"
    
    local out_dir="$TEST_TMP/ref_test"
    
    if python3 "$SCRIPT_DIR/generate_reference.py" \
        --output-dir "$out_dir" \
        --name test_ref > /dev/null 2>&1; then
        
        # Check FASTA exists and has content
        if [[ -f "$out_dir/test_ref.fa" ]]; then
            local size=$(stat -f%z "$out_dir/test_ref.fa" 2>/dev/null || stat -c%s "$out_dir/test_ref.fa")
            if [[ "$size" -gt 1000 ]]; then
                pass "FASTA generated (${size} bytes)"
            else
                fail "FASTA too small (${size} bytes)"
            fi
        else
            fail "FASTA not created"
        fi
        
        # Check manifest
        if [[ -f "$out_dir/test_ref.manifest" ]]; then
            pass "Manifest generated"
        else
            fail "Manifest not created"
        fi
        
        # Check truth file
        if [[ -f "$out_dir/test_ref_variants.truth" ]]; then
            pass "Truth file generated"
        else
            fail "Truth file not created"
        fi
    else
        fail "Reference generation failed"
    fi
}

test_reference_determinism() {
    echo "Testing: reference generator is deterministic"
    
    local out1="$TEST_TMP/ref_det1"
    local out2="$TEST_TMP/ref_det2"
    
    python3 "$SCRIPT_DIR/generate_reference.py" \
        --output-dir "$out1" --name det --seed 42 > /dev/null 2>&1
    python3 "$SCRIPT_DIR/generate_reference.py" \
        --output-dir "$out2" --name det --seed 42 > /dev/null 2>&1
    
    if diff -q "$out1/det.fa" "$out2/det.fa" > /dev/null 2>&1; then
        pass "Same seed produces identical output"
    else
        fail "Same seed produces different output"
    fi
    
    python3 "$SCRIPT_DIR/generate_reference.py" \
        --output-dir "$out2" --name det --seed 43 > /dev/null 2>&1
    
    if ! diff -q "$out1/det.fa" "$out2/det.fa" > /dev/null 2>&1; then
        pass "Different seeds produce different output"
    else
        fail "Different seeds produce same output"
    fi
}

test_fastq_generator_exists() {
    echo "Testing: FASTQ generator script exists"
    
    if [[ -f "$SCRIPT_DIR/generate_fastq.py" ]]; then
        pass "generate_fastq.py exists"
    else
        fail "generate_fastq.py not found"
    fi
}

test_fastq_generator_help() {
    echo "Testing: FASTQ generator --help"
    
    if python3 "$SCRIPT_DIR/generate_fastq.py" --help > /dev/null 2>&1; then
        pass "--help works"
    else
        fail "--help failed"
    fi
}

test_fastq_generator_output() {
    echo "Testing: FASTQ generator produces valid paired output"
    
    # First generate a reference
    local ref_dir="$TEST_TMP/fastq_ref"
    local out_dir="$TEST_TMP/fastq_out"
    
    python3 "$SCRIPT_DIR/generate_reference.py" \
        --output-dir "$ref_dir" --name ref > /dev/null 2>&1
    
    if python3 "$SCRIPT_DIR/generate_fastq.py" \
        --reference "$ref_dir/ref.fa" \
        --output-dir "$out_dir" \
        --profile good \
        --sample-name test \
        --num-reads 100 > /dev/null 2>&1; then
        
        # Check R1 exists
        if [[ -f "$out_dir/test_good_R1.fastq.gz" ]]; then
            pass "R1 FASTQ generated"
        else
            fail "R1 FASTQ not created"
        fi
        
        # Check R2 exists
        if [[ -f "$out_dir/test_good_R2.fastq.gz" ]]; then
            pass "R2 FASTQ generated"
        else
            fail "R2 FASTQ not created"
        fi
        
        # Check both have content (use gzip -dc for cross-platform compatibility)
        local r1_lines=$(gzip -dc "$out_dir/test_good_R1.fastq.gz" | wc -l)
        local r2_lines=$(gzip -dc "$out_dir/test_good_R2.fastq.gz" | wc -l)
        
        if [[ "$r1_lines" -eq "$r2_lines" && "$r1_lines" -eq 400 ]]; then
            pass "R1/R2 have matching read counts (100 reads × 4 lines)"
        else
            fail "R1/R2 line count mismatch or incorrect: R1=$r1_lines, R2=$r2_lines (expected 400)"
        fi
    else
        fail "FASTQ generation failed"
    fi
}

test_fastq_quality_profiles() {
    echo "Testing: FASTQ generator quality profiles"
    
    local ref_dir="$TEST_TMP/profile_ref"
    local out_dir="$TEST_TMP/profile_out"
    
    python3 "$SCRIPT_DIR/generate_reference.py" \
        --output-dir "$ref_dir" --name ref > /dev/null 2>&1
    
    for profile in good poor adapter mixed; do
        if python3 "$SCRIPT_DIR/generate_fastq.py" \
            --reference "$ref_dir/ref.fa" \
            --output-dir "$out_dir" \
            --profile "$profile" \
            --sample-name "profile_test" \
            --num-reads 50 > /dev/null 2>&1; then
            
            if [[ -f "$out_dir/profile_test_${profile}_R1.fastq.gz" ]]; then
                pass "Profile '$profile' generates valid output"
            else
                fail "Profile '$profile' did not create expected file"
            fi
        else
            fail "Profile '$profile' generation failed"
        fi
    done
}

test_smoke_harness_exists() {
    echo "Testing: smoke harness script exists"
    
    if [[ -f "$SCRIPT_DIR/run_smoke.sh" ]]; then
        pass "run_smoke.sh exists"
    else
        fail "run_smoke.sh not found"
    fi
}

test_smoke_harness_help() {
    echo "Testing: smoke harness --help"
    
    if bash "$SCRIPT_DIR/run_smoke.sh" --help > /dev/null 2>&1; then
        pass "--help works"
    else
        fail "--help failed"
    fi
}

test_smoke_harness_syntax() {
    echo "Testing: smoke harness shell syntax"
    
    if bash -n "$SCRIPT_DIR/run_smoke.sh" 2>/dev/null; then
        pass "Shell syntax is valid"
    else
        fail "Shell syntax errors detected"
    fi
}

test_makefile_smoke_targets() {
    echo "Testing: Makefile smoke targets exist"
    
    cd "$ROOT_DIR"
    
    local targets=("smoke" "smoke-mock" "smoke-real" "smoke-poor" "smoke-all")
    
    for target in "${targets[@]}"; do
        if make -n "$target" > /dev/null 2>&1; then
            pass "Makefile target '$target' exists"
        else
            fail "Makefile target '$target' missing"
        fi
    done
}

# =============================================================================
# Main
# =============================================================================

main() {
    echo "=================================================="
    echo "Smoke Test Component Unit Tests"
    echo "=================================================="
    echo
    
    trap teardown EXIT
    setup
    
    test_reference_generator_exists
    test_reference_generator_help
    test_reference_generator_output
    test_reference_determinism
    echo
    
    test_fastq_generator_exists
    test_fastq_generator_help
    test_fastq_generator_output
    test_fastq_quality_profiles
    echo
    
    test_smoke_harness_exists
    test_smoke_harness_help
    test_smoke_harness_syntax
    echo
    
    test_makefile_smoke_targets
    echo
    
    echo "=================================================="
    echo "Results: $TESTS_PASSED passed, $TESTS_FAILED failed"
    echo "=================================================="
    
    if [[ "$TESTS_FAILED" -gt 0 ]]; then
        exit 1
    fi
    exit 0
}

main "$@"
