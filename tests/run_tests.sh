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

# Create deterministic mock tools for edge-case unit tests
create_mock_tools() {
    local mock_bin="$TEMP_TEST_DIR/mock_bin"
    mkdir -p "$mock_bin"

    cat > "$mock_bin/fastqc" << 'EOF'
#!/bin/bash
set -e

outdir=""
input=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir)
            outdir="$2"
            shift 2
            ;;
        --threads)
            shift 2
            ;;
        --extract|--quiet)
            shift
            ;;
        *)
            input="$1"
            shift
            ;;
    esac
done

if [[ "${FASTQC_MOCK_MODE:-ok}" == "bad_exit" ]]; then
    exit 1
fi

if [[ -z "$outdir" || -z "$input" ]]; then
    exit 1
fi

base=$(basename "$input")
base="${base%.fq.gz}"
base="${base%.fastq.gz}"

if [[ "${FASTQC_MOCK_MODE:-ok}" == "missing_output" ]]; then
    exit 0
fi

mkdir -p "$outdir/${base}_fastqc"

if [[ "${FASTQC_MOCK_MODE:-ok}" == "low_quality" ]]; then
    cat > "$outdir/${base}_fastqc/fastqc_data.txt" << 'DATA'
Total Sequences	10000
Sequences flagged as poor quality	8500
Sequence length	100
%GC	45
>>Per base sequence quality	FAIL
>>Per sequence quality scores	FAIL
>>Per base sequence content	FAIL
>>Per sequence GC content	FAIL
>>Adapter Content	FAIL
>>Kmer Content	FAIL
DATA
else
    cat > "$outdir/${base}_fastqc/fastqc_data.txt" << 'DATA'
Total Sequences	10000
Sequences flagged as poor quality	10
Sequence length	100
%GC	45
>>Per base sequence quality	PASS
>>Per sequence quality scores	PASS
>>Per base sequence content	PASS
>>Per sequence GC content	PASS
>>Adapter Content	WARN
>>Kmer Content	PASS
DATA
fi

cat > "$outdir/${base}_fastqc.html" << HTML
<html><body>mock fastqc report</body></html>
HTML
EOF

    cat > "$mock_bin/fastp" << 'EOF'
#!/bin/bash
set -e

input_r1=""
input_r2=""
out_r1=""
out_r2=""
html=""
json=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i)
            input_r1="$2"
            shift 2
            ;;
        -I)
            input_r2="$2"
            shift 2
            ;;
        -o)
            out_r1="$2"
            shift 2
            ;;
        -O)
            out_r2="$2"
            shift 2
            ;;
        --html)
            html="$2"
            shift 2
            ;;
        --json)
            json="$2"
            shift 2
            ;;
        --thread|--qualified_quality_phred|--unqualified_percent_limit|--length_required|--cut_tail_window_size|--cut_tail_mean_quality)
            shift 2
            ;;
        --detect_adapter_for_pe|--correction|--cut_tail|--overrepresentation_analysis)
            shift
            ;;
        *)
            shift
            ;;
    esac
done

mode="${FASTP_MOCK_MODE:-normal}"

if [[ "$mode" == "bad_exit" ]]; then
    exit 1
fi

if [[ "$mode" != "missing_output" ]]; then
    mkdir -p "$(dirname "$out_r1")" "$(dirname "$out_r2")" "$(dirname "$html")" "$(dirname "$json")"
    cp "$input_r1" "$out_r1"
    cp "$input_r2" "$out_r2"
    cat > "$html" << HTML
<html><body>mock fastp report</body></html>
HTML
fi

before=10000
after=8500
q30_before=0.85
q30_after=0.90

if [[ "$mode" == "overtrim" ]]; then
    after=200
    q30_before=0.25
    q30_after=0.60
elif [[ "$mode" == "zero_reads" ]]; then
    after=0
    q30_before=0.20
    q30_after=0.10
elif [[ "$mode" == "low_quality_ok" ]]; then
    after=7000
    q30_before=0.10
    q30_after=0.55
elif [[ "$mode" == "missing_output" ]]; then
    # Simulate tool claiming success but producing no files
    exit 0
fi

cat > "$json" << JSON
{
  "summary": {
    "before_filtering": {
      "total_reads": $before,
      "q30_rate": $q30_before
    },
    "after_filtering": {
      "total_reads": $after,
      "q30_rate": $q30_after
    }
  }
}
JSON

exit 0
EOF

    chmod +x "$mock_bin/fastqc" "$mock_bin/fastp"
    echo "$mock_bin"
}

write_mock_fastq() {
    local output="$1"
    local reads="${2:-200}"
    local length="${3:-100}"

    python3 - "$output" "$reads" "$length" << 'EOF'
import gzip
import random
import sys

out = sys.argv[1]
reads = int(sys.argv[2])
length = int(sys.argv[3])

bases = "ACGT"
qual = "I" * length

with gzip.open(out, "wt") as f:
    for i in range(reads):
        seq = "".join(random.choice(bases) for _ in range(length))
        f.write(f"@read_{i}\n{seq}\n+\n{qual}\n")
EOF
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

    # Test 6: Dependency validator script exists and has correct syntax
    info "Testing dependency validator syntax..."
    if [[ -f "$ROOT_DIR/scripts/validate_deps.sh" ]]; then
        if bash -n "$ROOT_DIR/scripts/validate_deps.sh" 2>/dev/null; then
            echo "  ✓ Dependency validator syntax OK"
            tests_passed=$((tests_passed + 1))
        else
            echo "  ✗ Dependency validator has syntax errors"
            tests_failed=$((tests_failed + 1))
        fi
    else
        echo "  ✗ Dependency validator script not found"
        tests_failed=$((tests_failed + 1))
    fi

    # Test 7: Dependency validator help and version flags
    info "Testing dependency validator CLI..."
    if "$ROOT_DIR/scripts/validate_deps.sh" --help > /dev/null 2>&1; then
        echo "  ✓ Dependency validator --help works"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Dependency validator --help failed"
        tests_failed=$((tests_failed + 1))
    fi

    if "$ROOT_DIR/scripts/validate_deps.sh" --version > /dev/null 2>&1; then
        echo "  ✓ Dependency validator --version works"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Dependency validator --version failed"
        tests_failed=$((tests_failed + 1))
    fi

    # Test 8: Dependency validator rejects invalid arguments
    info "Testing dependency validator argument handling..."
    if "$ROOT_DIR/scripts/validate_deps.sh" --invalid-flag > /dev/null 2>&1; then
        echo "  ✗ Dependency validator accepted invalid flag"
        tests_failed=$((tests_failed + 1))
    else
        echo "  ✓ Dependency validator rejects invalid flags"
        tests_passed=$((tests_passed + 1))
    fi

    # Test 9: Dependency validator CI mode (should not fail for missing optional tools)
    info "Testing dependency validator CI mode..."
    if "$ROOT_DIR/scripts/validate_deps.sh" --ci-mode > /dev/null 2>&1; then
        echo "  ✓ Dependency validator CI mode works"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Dependency validator CI mode failed"
        tests_failed=$((tests_failed + 1))
    fi

    # Test 10: DEPENDENCIES.md exists
    info "Testing DEPENDENCIES.md existence..."
    if [[ -f "$ROOT_DIR/DEPENDENCIES.md" ]]; then
        echo "  ✓ DEPENDENCIES.md exists"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ DEPENDENCIES.md not found"
        tests_failed=$((tests_failed + 1))
    fi

    # Test 6: QC detects broken pairing
    info "Testing QC broken-pairing guardrail..."
    local mock_bin
    mock_bin=$(create_mock_tools)
    local original_path="$PATH"
    PATH="$mock_bin:$PATH"

    local qc_pair_dir="$TEMP_TEST_DIR/edge_qc_pairing"
    mkdir -p "$qc_pair_dir/raw" "$qc_pair_dir/out" "$qc_pair_dir/logs"
    write_mock_fastq "$qc_pair_dir/raw/sampleA_R1.fastq.gz" 400 100

    FASTQC_MOCK_MODE=ok "$ROOT_DIR/scripts/quality_control.sh" \
        --input-dir "$qc_pair_dir/raw" \
        --output-dir "$qc_pair_dir/out" \
        --log-dir "$qc_pair_dir/logs" \
        --threads 1 > "$qc_pair_dir/qc_pairing.log" 2>&1
    local rc=$?
    if [[ "$rc" -eq 12 ]]; then
        echo "  ✓ QC rejects broken FASTQ pairing (exit 12)"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ QC broken-pairing guardrail failed (exit $rc, expected 12)"
        tests_failed=$((tests_failed + 1))
    fi

    # Test 7: QC detects tiny FASTQ files
    info "Testing QC tiny-file guardrail..."
    local qc_tiny_dir="$TEMP_TEST_DIR/edge_qc_tiny"
    mkdir -p "$qc_tiny_dir/raw" "$qc_tiny_dir/out" "$qc_tiny_dir/logs"
    printf "tiny" | gzip -c > "$qc_tiny_dir/raw/sampleB_R1.fastq.gz"
    printf "tiny" | gzip -c > "$qc_tiny_dir/raw/sampleB_R2.fastq.gz"

    FASTQC_MOCK_MODE=ok "$ROOT_DIR/scripts/quality_control.sh" \
        --input-dir "$qc_tiny_dir/raw" \
        --output-dir "$qc_tiny_dir/out" \
        --log-dir "$qc_tiny_dir/logs" \
        --threads 1 > "$qc_tiny_dir/qc_tiny.log" 2>&1
    rc=$?
    if [[ "$rc" -eq 13 ]]; then
        echo "  ✓ QC rejects tiny/truncated FASTQ files (exit 13)"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ QC tiny-file guardrail failed (exit $rc, expected 13)"
        tests_failed=$((tests_failed + 1))
    fi

    # Test 8: QC checks for missing FastQC outputs
    info "Testing QC missing-output guardrail..."
    local qc_missing_dir="$TEMP_TEST_DIR/edge_qc_missing"
    mkdir -p "$qc_missing_dir/raw" "$qc_missing_dir/out" "$qc_missing_dir/logs"
    write_mock_fastq "$qc_missing_dir/raw/sampleC_R1.fastq.gz" 400 100
    write_mock_fastq "$qc_missing_dir/raw/sampleC_R2.fastq.gz" 400 100

    FASTQC_MOCK_MODE=missing_output "$ROOT_DIR/scripts/quality_control.sh" \
        --input-dir "$qc_missing_dir/raw" \
        --output-dir "$qc_missing_dir/out" \
        --log-dir "$qc_missing_dir/logs" \
        --threads 1 > "$qc_missing_dir/qc_missing.log" 2>&1
    rc=$?
    if [[ "$rc" -eq 15 ]]; then
        echo "  ✓ QC detects missing FastQC artifacts (exit 15)"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ QC missing-output guardrail failed (exit $rc, expected 15)"
        tests_failed=$((tests_failed + 1))
    fi

    # Test 9: Data cleaning detects broken pairing
    info "Testing cleaning broken-pairing guardrail..."
    local clean_pair_dir="$TEMP_TEST_DIR/edge_clean_pairing"
    mkdir -p "$clean_pair_dir/raw" "$clean_pair_dir/out"
    write_mock_fastq "$clean_pair_dir/raw/sampleD_R1.fastq.gz" 400 100

    FASTP_MOCK_MODE=normal "$ROOT_DIR/scripts/data_cleaning.sh" \
        --input-dir "$clean_pair_dir/raw" \
        --output-dir "$clean_pair_dir/out" \
        --sample-id sampleD \
        --threads 1 > "$clean_pair_dir/clean_pairing.log" 2>&1
    rc=$?
    if [[ "$rc" -eq 12 ]]; then
        echo "  ✓ Data cleaning rejects broken FASTQ pairing (exit 12)"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Data cleaning pairing guardrail failed (exit $rc, expected 12)"
        tests_failed=$((tests_failed + 1))
    fi

    # Test 10: Data cleaning catches over-aggressive trimming
    info "Testing cleaning over-trimming guardrail..."
    local clean_overtrim_dir="$TEMP_TEST_DIR/edge_clean_overtrim"
    mkdir -p "$clean_overtrim_dir/raw" "$clean_overtrim_dir/out"
    write_mock_fastq "$clean_overtrim_dir/raw/sampleE_R1.fastq.gz" 400 100
    write_mock_fastq "$clean_overtrim_dir/raw/sampleE_R2.fastq.gz" 400 100

    FASTP_MOCK_MODE=overtrim "$ROOT_DIR/scripts/data_cleaning.sh" \
        --input-dir "$clean_overtrim_dir/raw" \
        --output-dir "$clean_overtrim_dir/out" \
        --sample-id sampleE \
        --threads 1 > "$clean_overtrim_dir/clean_overtrim.log" 2>&1
    rc=$?
    if [[ "$rc" -eq 16 ]]; then
        echo "  ✓ Data cleaning detects over-aggressive trimming (exit 16)"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Data cleaning over-trimming guardrail failed (exit $rc, expected 16)"
        tests_failed=$((tests_failed + 1))
    fi

    # Test 11: Data cleaning catches missing output artifacts
    info "Testing cleaning missing-output guardrail..."
    local clean_missing_dir="$TEMP_TEST_DIR/edge_clean_missing"
    mkdir -p "$clean_missing_dir/raw" "$clean_missing_dir/out"
    write_mock_fastq "$clean_missing_dir/raw/sampleF_R1.fastq.gz" 400 100
    write_mock_fastq "$clean_missing_dir/raw/sampleF_R2.fastq.gz" 400 100

    FASTP_MOCK_MODE=missing_output "$ROOT_DIR/scripts/data_cleaning.sh" \
        --input-dir "$clean_missing_dir/raw" \
        --output-dir "$clean_missing_dir/out" \
        --sample-id sampleF \
        --threads 1 > "$clean_missing_dir/clean_missing.log" 2>&1
    rc=$?
    if [[ "$rc" -eq 15 ]]; then
        echo "  ✓ Data cleaning detects missing fastp artifacts (exit 15)"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Data cleaning missing-output guardrail failed (exit $rc, expected 15)"
        tests_failed=$((tests_failed + 1))
    fi

    # Test 12: Data cleaning succeeds on normal-quality mocked data
    info "Testing cleaning normal-path stability..."
    local clean_ok_dir="$TEMP_TEST_DIR/edge_clean_ok"
    mkdir -p "$clean_ok_dir/raw" "$clean_ok_dir/out"
    write_mock_fastq "$clean_ok_dir/raw/sampleG_R1.fastq.gz" 400 100
    write_mock_fastq "$clean_ok_dir/raw/sampleG_R2.fastq.gz" 400 100

    FASTP_MOCK_MODE=normal "$ROOT_DIR/scripts/data_cleaning.sh" \
        --input-dir "$clean_ok_dir/raw" \
        --output-dir "$clean_ok_dir/out" \
        --sample-id sampleG \
        --threads 1 > "$clean_ok_dir/clean_ok.log" 2>&1
    rc=$?
    if [[ "$rc" -eq 0 && -f "$clean_ok_dir/out/sampleG_clean_R1.fq.gz" && -f "$clean_ok_dir/out/sampleG_clean_R2.fq.gz" ]]; then
        echo "  ✓ Data cleaning normal path remains stable"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Data cleaning normal path failed (exit $rc)"
        tests_failed=$((tests_failed + 1))
    fi

    PATH="$original_path"

    # ========================================
    # VEP Annotation Unit Tests
    # ========================================
    
    log "Running VEP annotation unit tests..."
    
    # Create annotation mock tools
    local ann_mock_bin="$TEMP_TEST_DIR/ann_mock_bin"
    mkdir -p "$ann_mock_bin"
    
    # Mock bcftools that works
    cat > "$ann_mock_bin/bcftools" << 'BCFTOOLS_MOCK'
#!/bin/bash
case "$1" in
    view)
        if [[ "$2" == "-H" ]]; then
            # Output variant lines from mocked VCF
            if [[ "${VCF_MOCK_EMPTY:-false}" == "true" ]]; then
                exit 0
            fi
            echo "chr1	100	.	A	G	30	PASS	."
            echo "chr1	200	.	C	T	40	PASS	."
        elif [[ "$2" == "-h" ]]; then
            echo "##fileformat=VCFv4.2"
            echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"
        fi
        ;;
    index)
        touch "${@: -1}.csi"
        ;;
    query)
        echo "chr1	100	A	G	GENE1|missense|MODERATE"
        ;;
    --version)
        echo "bcftools 1.17"
        ;;
esac
exit 0
BCFTOOLS_MOCK
    chmod +x "$ann_mock_bin/bcftools"
    
    # Mock tabix
    cat > "$ann_mock_bin/tabix" << 'TABIX_MOCK'
#!/bin/bash
touch "$3.tbi" 2>/dev/null || touch "${@: -1}.tbi"
exit 0
TABIX_MOCK
    chmod +x "$ann_mock_bin/tabix"
    
    # Test annotation-1: Script rejects invalid arguments
    info "Testing annotation argument rejection..."
    if "$ROOT_DIR/scripts/vep_annotation.sh" --invalid-flag > /dev/null 2>&1; then
        echo "  ✗ vep_annotation.sh accepted invalid flag"
        tests_failed=$((tests_failed + 1))
    else
        echo "  ✓ vep_annotation.sh rejects invalid flags"
        tests_passed=$((tests_passed + 1))
    fi
    
    # Test annotation-2: Help includes exit codes
    info "Testing annotation help includes exit codes..."
    if "$ROOT_DIR/scripts/vep_annotation.sh" --help 2>/dev/null | grep -q "EXIT CODES"; then
        echo "  ✓ vep_annotation.sh help documents exit codes"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ vep_annotation.sh help missing exit codes"
        tests_failed=$((tests_failed + 1))
    fi
    
    # Test annotation-3: Missing VEP binary detected (exit 10)
    info "Testing annotation detects missing VEP..."
    local ann_novep_dir="$TEMP_TEST_DIR/ann_novep"
    mkdir -p "$ann_novep_dir/results/variants" "$ann_novep_dir/config" "$ann_novep_dir/logs"
    cp "$ROOT_DIR/config/default.conf" "$ann_novep_dir/config/"
    
    # Create a minimal valid VCF
    cat > "$ann_novep_dir/results/variants/test_filtered.vcf" << 'MOCK_VCF'
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	G	30	PASS	DP=20
MOCK_VCF
    gzip -f "$ann_novep_dir/results/variants/test_filtered.vcf"
    
    # Use PATH without vep to test missing binary detection
    # Note: prepend mock path to keep system binaries available
    (cd "$ann_novep_dir" && PATH="$ann_mock_bin:$PATH" "$ROOT_DIR/scripts/vep_annotation.sh" \
        --input "results/variants/test_filtered.vcf.gz" \
        --output-dir "results/annotation" \
        --sample-id test \
        > "novep.log" 2>&1)
    rc=$?
    if [[ "$rc" -eq 10 ]]; then
        echo "  ✓ Annotation detects missing VEP binary (exit 10)"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Annotation missing-VEP detection failed (exit $rc, expected 10)"
        tests_failed=$((tests_failed + 1))
    fi
    
    # Test annotation-4: Missing bcftools detected (exit 11)
    info "Testing annotation detects missing bcftools..."
    local ann_nobcf_dir="$TEMP_TEST_DIR/ann_nobcf"
    mkdir -p "$ann_nobcf_dir/results/variants" "$ann_nobcf_dir/config" "$ann_nobcf_dir/logs"
    cp "$ROOT_DIR/config/default.conf" "$ann_nobcf_dir/config/"
    cp "$ann_novep_dir/results/variants/test_filtered.vcf.gz" "$ann_nobcf_dir/results/variants/"
    
    # Create minimal mock directory with only VEP - no bcftools
    mkdir -p "$ann_mock_bin/nobcf_clean"
    cat > "$ann_mock_bin/nobcf_clean/vep" << 'VEP_MOCK'
#!/bin/bash
echo "ensembl-vep : 110"
exit 0
VEP_MOCK
    chmod +x "$ann_mock_bin/nobcf_clean/vep"
    
    # DO NOT create bcftools mock - we want `command -v bcftools` to fail
    # Use a PATH that excludes where bcftools might be installed
    (cd "$ann_nobcf_dir" && PATH="$ann_mock_bin/nobcf_clean:/usr/bin:/bin" "$ROOT_DIR/scripts/vep_annotation.sh" \
        --input "results/variants/test_filtered.vcf.gz" \
        --output-dir "results/annotation" \
        --sample-id test \
        > "nobcf.log" 2>&1)
    rc=$?
    if [[ "$rc" -eq 11 ]]; then
        echo "  ✓ Annotation detects missing bcftools (exit 11)"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Annotation missing-bcftools detection failed (exit $rc, expected 11)"
        tests_failed=$((tests_failed + 1))
    fi
    
    # Test annotation-5: Missing input VCF detected (exit 12)
    info "Testing annotation detects missing input VCF..."
    local ann_novcf_dir="$TEMP_TEST_DIR/ann_novcf"
    mkdir -p "$ann_novcf_dir/results/variants" "$ann_novcf_dir/config" "$ann_novcf_dir/logs"
    cp "$ROOT_DIR/config/default.conf" "$ann_novcf_dir/config/"
    
    # Create mocks for both vep and bcftools
    mkdir -p "$ann_mock_bin/full"
    cat > "$ann_mock_bin/full/vep" << 'VEP_MOCK'
#!/bin/bash
if [[ "$1" == "--version" || "$1" == "--help" ]]; then
    echo "ensembl-vep : 110"
    exit 0
fi
exit 0
VEP_MOCK
    chmod +x "$ann_mock_bin/full/vep"
    cp "$ann_mock_bin/bcftools" "$ann_mock_bin/full/"
    cp "$ann_mock_bin/tabix" "$ann_mock_bin/full/"
    
    (cd "$ann_novcf_dir" && PATH="$ann_mock_bin/full:$PATH" "$ROOT_DIR/scripts/vep_annotation.sh" \
        --input "results/variants/nonexistent.vcf.gz" \
        --output-dir "results/annotation" \
        --sample-id test \
        > "novcf.log" 2>&1)
    rc=$?
    if [[ "$rc" -eq 12 ]]; then
        echo "  ✓ Annotation detects missing input VCF (exit 12)"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Annotation missing-VCF detection failed (exit $rc, expected 12)"
        tests_failed=$((tests_failed + 1))
    fi
    
    # Test annotation-6: Malformed VCF detected (exit 13)
    info "Testing annotation detects malformed VCF..."
    local ann_badvcf_dir="$TEMP_TEST_DIR/ann_badvcf"
    mkdir -p "$ann_badvcf_dir/results/variants" "$ann_badvcf_dir/config" "$ann_badvcf_dir/logs"
    cp "$ROOT_DIR/config/default.conf" "$ann_badvcf_dir/config/"
    
    # Create a malformed VCF (not starting with ##fileformat)
    echo "this is not a vcf file" | gzip > "$ann_badvcf_dir/results/variants/bad_filtered.vcf.gz"
    
    (cd "$ann_badvcf_dir" && PATH="$ann_mock_bin/full:$PATH" "$ROOT_DIR/scripts/vep_annotation.sh" \
        --input "results/variants/bad_filtered.vcf.gz" \
        --output-dir "results/annotation" \
        --sample-id test \
        > "badvcf.log" 2>&1)
    rc=$?
    if [[ "$rc" -eq 13 ]]; then
        echo "  ✓ Annotation detects malformed VCF (exit 13)"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Annotation malformed-VCF detection failed (exit $rc, expected 13)"
        tests_failed=$((tests_failed + 1))
    fi
    
    # Test annotation-7: Empty VCF detected (exit 14)
    info "Testing annotation detects empty VCF..."
    local ann_empty_dir="$TEMP_TEST_DIR/ann_empty"
    mkdir -p "$ann_empty_dir/results/variants" "$ann_empty_dir/config" "$ann_empty_dir/logs"
    cp "$ROOT_DIR/config/default.conf" "$ann_empty_dir/config/"
    
    # Create VCF with header but no variants
    cat > "$ann_empty_dir/results/variants/empty_filtered.vcf" << 'EMPTY_VCF'
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
EMPTY_VCF
    gzip -f "$ann_empty_dir/results/variants/empty_filtered.vcf"
    
    # Use a bcftools that reports empty
    cat > "$ann_mock_bin/full/bcftools_empty" << 'BCF_EMPTY'
#!/bin/bash
case "$1" in
    view)
        if [[ "$2" == "-H" ]]; then
            # No output = empty VCF
            exit 0
        elif [[ "$2" == "-h" ]]; then
            echo "##fileformat=VCFv4.2"
            echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"
        fi
        ;;
    --version)
        echo "bcftools 1.17"
        ;;
esac
exit 0
BCF_EMPTY
    chmod +x "$ann_mock_bin/full/bcftools_empty"
    cp "$ann_mock_bin/full/bcftools" "$ann_mock_bin/full/bcftools_orig"
    cp "$ann_mock_bin/full/bcftools_empty" "$ann_mock_bin/full/bcftools"
    
    (cd "$ann_empty_dir" && PATH="$ann_mock_bin/full:$PATH" "$ROOT_DIR/scripts/vep_annotation.sh" \
        --input "results/variants/empty_filtered.vcf.gz" \
        --output-dir "results/annotation" \
        --sample-id test \
        > "empty.log" 2>&1)
    rc=$?
    if [[ "$rc" -eq 14 ]]; then
        echo "  ✓ Annotation detects empty VCF (exit 14)"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Annotation empty-VCF detection failed (exit $rc, expected 14)"
        tests_failed=$((tests_failed + 1))
    fi
    
    # Restore bcftools
    cp "$ann_mock_bin/full/bcftools_orig" "$ann_mock_bin/full/bcftools"
    
    # Test annotation-8: --skip-empty-check allows empty VCF
    info "Testing annotation --skip-empty-check flag..."
    # Ensure empty bcftools mock is in place
    cp "$ann_mock_bin/full/bcftools_empty" "$ann_mock_bin/full/bcftools"
    (cd "$ann_empty_dir" && PATH="$ann_mock_bin/full:$PATH" "$ROOT_DIR/scripts/vep_annotation.sh" \
        --input "results/variants/empty_filtered.vcf.gz" \
        --output-dir "results/annotation2" \
        --sample-id test \
        --skip-empty-check \
        --dry-run \
        > "skip_empty.log" 2>&1)
    rc=$?
    cp "$ann_mock_bin/full/bcftools_orig" "$ann_mock_bin/full/bcftools"
    if [[ "$rc" -eq 0 ]]; then
        echo "  ✓ Annotation --skip-empty-check allows empty VCF"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Annotation --skip-empty-check failed (exit $rc, expected 0)"
        tests_failed=$((tests_failed + 1))
    fi
    
    # Test annotation-9: --validate-only mode works
    info "Testing annotation --validate-only mode..."
    local ann_validate_dir="$TEMP_TEST_DIR/ann_validate"
    mkdir -p "$ann_validate_dir/results/variants" "$ann_validate_dir/config" "$ann_validate_dir/logs"
    cp "$ROOT_DIR/config/default.conf" "$ann_validate_dir/config/"
    cp "$ann_novep_dir/results/variants/test_filtered.vcf.gz" "$ann_validate_dir/results/variants/"
    
    (cd "$ann_validate_dir" && PATH="$ann_mock_bin/full:$PATH" "$ROOT_DIR/scripts/vep_annotation.sh" \
        --input "results/variants/test_filtered.vcf.gz" \
        --output-dir "results/annotation" \
        --sample-id test \
        --validate-only \
        > "validate.log" 2>&1)
    rc=$?
    if [[ "$rc" -eq 0 ]]; then
        echo "  ✓ Annotation --validate-only mode works"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Annotation --validate-only failed (exit $rc, expected 0)"
        tests_failed=$((tests_failed + 1))
    fi
    
    # Test annotation-10: --dry-run mode works
    info "Testing annotation --dry-run mode..."
    (cd "$ann_validate_dir" && PATH="$ann_mock_bin/full:$PATH" "$ROOT_DIR/scripts/vep_annotation.sh" \
        --input "results/variants/test_filtered.vcf.gz" \
        --output-dir "results/annotation" \
        --sample-id test \
        --dry-run \
        > "dryrun.log" 2>&1)
    rc=$?
    if [[ "$rc" -eq 0 ]] && grep -q "Would run VEP" "$ann_validate_dir/dryrun.log"; then
        echo "  ✓ Annotation --dry-run mode works"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Annotation --dry-run failed (exit $rc)"
        tests_failed=$((tests_failed + 1))
    fi
    
    # Test annotation-11: Version flag works
    info "Testing annotation --version flag..."
    if "$ROOT_DIR/scripts/vep_annotation.sh" --version 2>/dev/null | grep -q "2.1"; then
        echo "  ✓ Annotation --version works"
        tests_passed=$((tests_passed + 1))
    else
        echo "  ✗ Annotation --version failed"
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
