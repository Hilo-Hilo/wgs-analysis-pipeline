#!/bin/bash
#
# WGS Pipeline End-to-End Smoke Test Harness
#
# Runs a complete pipeline execution on synthetic data:
#   1. Generate tiny reference genome + indexes
#   2. Generate synthetic FASTQ pairs (good + poor quality)
#   3. Run QC → Cleaning → Alignment → Variant Calling
#   4. Validate outputs and report pass/fail
#
# This harness is designed to:
#   - Run quickly (<5 minutes with mocked tools, <30 min real)
#   - Be fully deterministic (fixed seeds)
#   - Support both mocked and real tool execution
#   - Provide clear pass/fail signals for CI
#
# Usage:
#   ./tests/smoke/run_smoke.sh [OPTIONS]
#
# Options:
#   --mock              Use mock tools (fast, no bioinformatics deps required)
#   --real              Use real tools (requires bwa, samtools, bcftools, fastp, fastqc)
#   --profile PROFILE   Quality profile: good, poor, adapter, mixed (default: good)
#   --keep              Keep output artifacts after completion
#   --verbose           Enable verbose output
#   --help              Show this help message
#

set -e
set -o pipefail

# =============================================================================
# Configuration
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
SMOKE_DATA_DIR="$SCRIPT_DIR/data"
SMOKE_OUTPUT_DIR="/tmp/wgs-smoke-$$"

# Defaults
USE_MOCK=true
QUALITY_PROFILE="good"
KEEP_ARTIFACTS=false
VERBOSE=false

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

# Exit codes
EX_OK=0
EX_USAGE=2
EX_SETUP_FAIL=10
EX_QC_FAIL=11
EX_CLEAN_FAIL=12
EX_ALIGN_FAIL=13
EX_VARIANT_FAIL=14
EX_VALIDATE_FAIL=15

# Timers
STEP_START=0
TOTAL_START=$(date +%s)

# =============================================================================
# Functions
# =============================================================================

show_help() {
    cat << EOF
WGS Pipeline End-to-End Smoke Test Harness

USAGE:
    $0 [OPTIONS]

OPTIONS:
    --mock              Use mock tools (fast, no dependencies required)
    --real              Use real tools (requires bioinformatics stack)
    --profile PROFILE   Quality profile: good, poor, adapter, mixed (default: good)
    --keep              Keep output artifacts after completion
    --verbose           Enable verbose output
    --help              Show this help message

EXAMPLES:
    # Quick smoke test with mocks
    $0 --mock

    # Full real pipeline test
    $0 --real --verbose

    # Test edge case (poor quality data)
    $0 --real --profile poor --keep

EXIT CODES:
    0   All smoke tests passed
    2   Invalid usage
    10  Setup/data generation failed
    11  Quality control step failed
    12  Data cleaning step failed
    13  Alignment step failed
    14  Variant calling step failed
    15  Output validation failed
EOF
}

log() {
    echo -e "${GREEN}[SMOKE]${NC} $1"
}

error() {
    echo -e "${RED}[SMOKE ERROR]${NC} $1" >&2
}

warning() {
    echo -e "${YELLOW}[SMOKE WARN]${NC} $1"
}

info() {
    if [[ "$VERBOSE" == "true" ]]; then
        echo -e "${BLUE}[SMOKE INFO]${NC} $1"
    fi
}

step() {
    STEP_START=$(date +%s)
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${CYAN}▶ $1${NC}"
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
}

step_done() {
    local elapsed=$(($(date +%s) - STEP_START))
    echo -e "${GREEN}✓ Step completed in ${elapsed}s${NC}"
    echo
}

step_fail() {
    local elapsed=$(($(date +%s) - STEP_START))
    echo -e "${RED}✗ Step failed after ${elapsed}s${NC}"
    echo
}

cleanup() {
    if [[ "$KEEP_ARTIFACTS" != "true" && -d "$SMOKE_OUTPUT_DIR" ]]; then
        info "Cleaning up: $SMOKE_OUTPUT_DIR"
        rm -rf "$SMOKE_OUTPUT_DIR"
    fi
}

# =============================================================================
# Argument Parsing
# =============================================================================

while [[ $# -gt 0 ]]; do
    case $1 in
        --mock)
            USE_MOCK=true
            shift
            ;;
        --real)
            USE_MOCK=false
            shift
            ;;
        --profile)
            QUALITY_PROFILE="$2"
            shift 2
            ;;
        --keep)
            KEEP_ARTIFACTS=true
            shift
            ;;
        --verbose|-v)
            VERBOSE=true
            shift
            ;;
        --help|-h)
            show_help
            exit $EX_OK
            ;;
        *)
            error "Unknown option: $1"
            show_help
            exit $EX_USAGE
            ;;
    esac
done

# Validate profile
case "$QUALITY_PROFILE" in
    good|poor|adapter|mixed)
        ;;
    *)
        error "Invalid profile: $QUALITY_PROFILE"
        error "Valid profiles: good, poor, adapter, mixed"
        exit $EX_USAGE
        ;;
esac

# =============================================================================
# Setup
# =============================================================================

trap cleanup EXIT

log "WGS Pipeline E2E Smoke Test"
log "==========================="
log "Mode: $([ "$USE_MOCK" == "true" ] && echo 'MOCK (fast)' || echo 'REAL (full)')"
log "Profile: $QUALITY_PROFILE"
log "Output: $SMOKE_OUTPUT_DIR"
log "Keep artifacts: $KEEP_ARTIFACTS"
echo

# Create output structure
mkdir -p "$SMOKE_OUTPUT_DIR"/{data/raw,data/processed,data/reference,results,logs}

# =============================================================================
# Step 1: Generate Reference Genome
# =============================================================================

step "Step 1: Generate synthetic reference genome"

REF_DIR="$SMOKE_OUTPUT_DIR/data/reference"
REF_FASTA="$REF_DIR/smoke_ref.fa"

if python3 "$SCRIPT_DIR/generate_reference.py" \
    --output-dir "$REF_DIR" \
    --name smoke_ref \
    $([ "$VERBOSE" == "true" ] && echo "--verbose"); then
    
    # Create indexes - either real or mock depending on mode
    if [[ "$USE_MOCK" == "true" ]]; then
        info "Creating mock reference indexes..."
        # Create mock BWA index files
        touch "$REF_FASTA.bwt" "$REF_FASTA.pac" "$REF_FASTA.ann" "$REF_FASTA.amb" "$REF_FASTA.sa"
        # Create mock samtools index
        awk '/^>/ {contig=substr($1,2); next} {len+=length($0)} END {}' "$REF_FASTA" > /dev/null
        # Generate proper faidx format (contig name, length, offset, bases per line, bytes per line)
        awk '/^>/ {if (contig) print contig"\t"len"\t"offset"\t80\t81"; contig=substr($1,2); len=0; offset=NR*81} 
             !/^>/ {len+=length($0)} 
             END {if (contig) print contig"\t"len"\t"1"\t80\t81"}' "$REF_FASTA" > "$REF_FASTA.fai" 2>/dev/null || touch "$REF_FASTA.fai"
    else
        # Real mode - use actual tools
        if command -v bwa &>/dev/null; then
            info "Indexing reference with BWA..."
            bwa index "$REF_FASTA" 2>/dev/null || warning "BWA indexing failed"
        fi
        
        if command -v samtools &>/dev/null; then
            info "Indexing reference with samtools..."
            samtools faidx "$REF_FASTA" 2>/dev/null || warning "samtools faidx failed"
        fi
    fi
    
    step_done
else
    step_fail
    error "Failed to generate reference genome"
    exit $EX_SETUP_FAIL
fi

# =============================================================================
# Step 2: Generate FASTQ Reads
# =============================================================================

step "Step 2: Generate synthetic FASTQ reads"

RAW_DIR="$SMOKE_OUTPUT_DIR/data/raw"

if python3 "$SCRIPT_DIR/generate_fastq.py" \
    --reference "$REF_FASTA" \
    --output-dir "$RAW_DIR" \
    --profile "$QUALITY_PROFILE" \
    --sample-name "smoke" \
    --num-reads 2000 \
    $([ "$VERBOSE" == "true" ] && echo "--verbose"); then
    
    step_done
else
    step_fail
    error "Failed to generate FASTQ files"
    exit $EX_SETUP_FAIL
fi

# Setup mock tools if requested
MOCK_BIN=""
if [[ "$USE_MOCK" == "true" ]]; then
    step "Step 2b: Setting up mock tools"
    MOCK_BIN="$SMOKE_OUTPUT_DIR/mock_bin"
    mkdir -p "$MOCK_BIN"
    
    # Mock fastqc
    cat > "$MOCK_BIN/fastqc" << 'MOCK'
#!/bin/bash
outdir="."
while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir) outdir="$2"; shift 2 ;;
        --threads|--extract|--quiet) shift ;;
        *) input="$1"; shift ;;
    esac
done
base=$(basename "$input" | sed 's/\.f.*q\.gz$//')
mkdir -p "$outdir/${base}_fastqc"
cat > "$outdir/${base}_fastqc/fastqc_data.txt" << DATA
Total Sequences	2000
Sequences flagged as poor quality	10
Sequence length	150
%GC	45
>>Per base sequence quality	PASS
>>Per sequence quality scores	PASS
>>Per base sequence content	PASS
>>Per sequence GC content	PASS
>>Adapter Content	PASS
DATA
touch "$outdir/${base}_fastqc.html"
MOCK
    
    # Mock fastp
    cat > "$MOCK_BIN/fastp" << 'MOCK'
#!/bin/bash
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) r1_in="$2"; shift 2 ;;
        -I) r2_in="$2"; shift 2 ;;
        -o) r1_out="$2"; shift 2 ;;
        -O) r2_out="$2"; shift 2 ;;
        --html) html="$2"; shift 2 ;;
        --json) json="$2"; shift 2 ;;
        *) shift ;;
    esac
done
cp "$r1_in" "$r1_out"
cp "$r2_in" "$r2_out"
mkdir -p "$(dirname "$html")" "$(dirname "$json")"
echo "<html><body>Mock fastp report</body></html>" > "$html"
cat > "$json" << JSON
{"summary":{"before_filtering":{"total_reads":2000,"q30_rate":0.85},"after_filtering":{"total_reads":1800,"q30_rate":0.92}}}
JSON
MOCK
    
    # Mock bwa
    cat > "$MOCK_BIN/bwa" << 'MOCK'
#!/bin/bash
case "$1" in
    index) exit 0 ;;
    mem)
        # Generate minimal SAM output
        shift
        while [[ $# -gt 0 && "$1" == -* ]]; do shift; [[ -n "$1" && "$1" != -* ]] && shift; done
        ref="$1"; r1="$2"; r2="$3"
        echo "@HD	VN:1.6	SO:unsorted"
        echo "@SQ	SN:chr1_mini	LN:40000"
        echo "@SQ	SN:chr2_mini	LN:30000"
        echo "@SQ	SN:chr3_mini	LN:20000"
        # Generate mock aligned reads
        for i in $(seq 1 100); do
            pos=$((i * 100))
            echo "read_${i}	99	chr1_mini	${pos}	60	150M	=	$((pos+200))	350	$(printf 'A%.0s' {1..150})	$(printf 'I%.0s' {1..150})"
        done
        ;;
esac
MOCK
    
    # Mock samtools
    cat > "$MOCK_BIN/samtools" << 'MOCK'
#!/bin/bash
case "$1" in
    faidx) exit 0 ;;
    sort)
        shift
        while [[ $# -gt 0 ]]; do
            case "$1" in
                -o) out="$2"; shift 2 ;;
                -@|-m) shift 2 ;;
                -) cat > /dev/null; break ;;
                *) shift ;;
            esac
        done
        # Create minimal BAM-like output (actually just a marker file)
        echo "MOCK_BAM_HEADER" | gzip > "$out"
        ;;
    index)
        touch "${2}.bai"
        ;;
    flagstat)
        echo "200 + 0 in total (QC-passed reads + QC-failed reads)"
        echo "200 + 0 mapped (100.00% : N/A)"
        ;;
    --version)
        echo "samtools 1.20 (mock)"
        ;;
esac
MOCK
    
    # Mock bcftools
    cat > "$MOCK_BIN/bcftools" << 'MOCK'
#!/bin/bash
case "$1" in
    mpileup)
        # Generate mock pileup
        echo "##fileformat=VCFv4.2"
        echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"
        echo "chr1_mini	1000	.	A	G	30	.	DP=20"
        echo "chr1_mini	5000	.	T	C	35	.	DP=25"
        ;;
    call)
        # Pass through
        cat
        ;;
    filter)
        shift
        while [[ $# -gt 0 ]]; do
            case "$1" in
                -O) shift 2 ;;
                -o) out="$2"; shift 2 ;;
                -e) shift 2 ;;
                *) input="$1"; shift ;;
            esac
        done
        if [[ -n "$input" ]]; then
            zcat "$input" 2>/dev/null || cat "$input"
        else
            cat
        fi | gzip > "$out"
        ;;
    index)
        touch "${2}.csi"
        ;;
    view)
        shift
        while [[ "$1" == -* ]]; do shift; done
        zcat "$1" 2>/dev/null || cat "$1"
        ;;
    stats)
        echo "# bcftools stats (mock)"
        echo "SN	0	number of samples:	1"
        echo "SN	0	number of SNPs:	2"
        ;;
esac
MOCK
    
    chmod +x "$MOCK_BIN"/*
    export PATH="$MOCK_BIN:$PATH"
    
    step_done
fi

# =============================================================================
# Step 3: Run Quality Control
# =============================================================================

step "Step 3: Quality Control (FastQC)"

QC_OUTPUT="$SMOKE_OUTPUT_DIR/results/quality_control"
mkdir -p "$QC_OUTPUT"

SAMPLE_ID="smoke_${QUALITY_PROFILE}"

if "$ROOT_DIR/scripts/quality_control.sh" \
    --input-dir "$RAW_DIR" \
    --output-dir "$QC_OUTPUT" \
    --log-dir "$SMOKE_OUTPUT_DIR/logs" \
    --threads 2 \
    $([ "$VERBOSE" == "true" ] && echo "--verbose"); then
    
    step_done
else
    step_fail
    
    # For poor quality profile, some warnings are expected
    if [[ "$QUALITY_PROFILE" == "poor" ]]; then
        warning "QC warnings expected for poor quality profile"
    else
        error "Quality control step failed"
        exit $EX_QC_FAIL
    fi
fi

# =============================================================================
# Step 4: Run Data Cleaning
# =============================================================================

step "Step 4: Data Cleaning (fastp)"

CLEAN_OUTPUT="$SMOKE_OUTPUT_DIR/data/processed"

if "$ROOT_DIR/scripts/data_cleaning.sh" \
    --input-dir "$RAW_DIR" \
    --output-dir "$CLEAN_OUTPUT" \
    --sample-id "$SAMPLE_ID" \
    --threads 2 \
    $([ "$VERBOSE" == "true" ] && echo "--verbose"); then
    
    step_done
else
    rc=$?
    step_fail
    
    # Over-trimming expected for very poor quality
    if [[ "$QUALITY_PROFILE" == "poor" && "$rc" -eq 16 ]]; then
        warning "Over-trimming expected for poor quality profile (exit 16)"
        warning "Continuing smoke test to validate error handling..."
    else
        error "Data cleaning step failed (exit $rc)"
        exit $EX_CLEAN_FAIL
    fi
fi

# =============================================================================
# Step 5: Run Alignment (if cleaned files exist)
# =============================================================================

CLEANED_R1="$CLEAN_OUTPUT/${SAMPLE_ID}_clean_R1.fq.gz"
CLEANED_R2="$CLEAN_OUTPUT/${SAMPLE_ID}_clean_R2.fq.gz"

if [[ -f "$CLEANED_R1" && -f "$CLEANED_R2" ]]; then
    step "Step 5: Alignment (BWA/samtools)"
    
    ALIGN_OUTPUT="$SMOKE_OUTPUT_DIR/results/alignment"
    mkdir -p "$ALIGN_OUTPUT"
    
    # Update config to use our reference
    export GRCH38_REFERENCE="$REF_FASTA"
    export PROCESSED_DATA_DIR="$CLEAN_OUTPUT"
    
    if "$ROOT_DIR/scripts/alignment.sh" \
        --input-dir "$CLEAN_OUTPUT" \
        --output-dir "$ALIGN_OUTPUT" \
        --reference "$REF_FASTA" \
        --sample-id "$SAMPLE_ID" \
        --threads 2 \
        $([ "$VERBOSE" == "true" ] && echo "--verbose"); then
        
        step_done
    else
        step_fail
        error "Alignment step failed"
        exit $EX_ALIGN_FAIL
    fi
    
    # ==========================================================================
    # Step 6: Run Variant Calling
    # ==========================================================================
    
    BAM_FILE="$ALIGN_OUTPUT/${SAMPLE_ID}_aligned_sorted.bam"
    
    if [[ -f "$BAM_FILE" ]]; then
        step "Step 6: Variant Calling (bcftools)"
        
        VARIANT_OUTPUT="$SMOKE_OUTPUT_DIR/results/variants"
        mkdir -p "$VARIANT_OUTPUT"
        
        if "$ROOT_DIR/scripts/variant_calling.sh" \
            --input "$BAM_FILE" \
            --output-dir "$VARIANT_OUTPUT" \
            --reference "$REF_FASTA" \
            --sample-id "$SAMPLE_ID" \
            --threads 2 \
            $([ "$VERBOSE" == "true" ] && echo "--verbose"); then
            
            step_done
        else
            step_fail
            error "Variant calling step failed"
            exit $EX_VARIANT_FAIL
        fi
    else
        warning "Skipping variant calling: BAM file not found"
    fi
else
    warning "Skipping alignment + variant calling: cleaned FASTQs not found"
    warning "(This is expected for --profile poor due to over-trimming)"
fi

# =============================================================================
# Step 7: Validate Outputs
# =============================================================================

step "Step 7: Validate outputs"

VALIDATION_PASSED=true
VALIDATION_WARNINGS=0

validate_file() {
    local path="$1"
    local desc="$2"
    local required="${3:-true}"
    
    if [[ -f "$path" ]]; then
        local size=$(stat -f%z "$path" 2>/dev/null || stat -c%s "$path" 2>/dev/null || echo 0)
        if [[ "$size" -gt 0 ]]; then
            info "✓ $desc: $path ($size bytes)"
            return 0
        else
            if [[ "$required" == "true" ]]; then
                error "✗ $desc: file exists but is empty: $path"
                VALIDATION_PASSED=false
            else
                warning "! $desc: file exists but is empty: $path"
                VALIDATION_WARNINGS=$((VALIDATION_WARNINGS + 1))
            fi
            return 1
        fi
    else
        if [[ "$required" == "true" ]]; then
            error "✗ $desc: file not found: $path"
            VALIDATION_PASSED=false
        else
            info "- $desc: skipped (not required)"
        fi
        return 1
    fi
}

# Validate reference generation
validate_file "$REF_FASTA" "Reference FASTA"
validate_file "$REF_DIR/smoke_ref.manifest" "Reference manifest"
validate_file "$REF_DIR/smoke_ref_variants.truth" "Variant truth file"

# Validate FASTQ generation
validate_file "$RAW_DIR/smoke_${QUALITY_PROFILE}_R1.fastq.gz" "Raw FASTQ R1"
validate_file "$RAW_DIR/smoke_${QUALITY_PROFILE}_R2.fastq.gz" "Raw FASTQ R2"

# Validate QC outputs
validate_file "$QC_OUTPUT/quality_control_summary.txt" "QC summary" "false"

# Validate cleaning outputs (may not exist for poor quality)
if [[ "$QUALITY_PROFILE" != "poor" ]]; then
    validate_file "$CLEANED_R1" "Cleaned FASTQ R1"
    validate_file "$CLEANED_R2" "Cleaned FASTQ R2"
    
    # Validate alignment outputs
    validate_file "$ALIGN_OUTPUT/${SAMPLE_ID}_aligned_sorted.bam" "Aligned BAM" "false"
    
    # Validate variant outputs
    validate_file "$VARIANT_OUTPUT/${SAMPLE_ID}_raw.vcf.gz" "Raw VCF" "false"
else
    info "- Skipping downstream validation for poor quality profile"
fi

if [[ "$VALIDATION_PASSED" == "true" ]]; then
    step_done
else
    step_fail
    exit $EX_VALIDATE_FAIL
fi

# =============================================================================
# Summary
# =============================================================================

TOTAL_TIME=$(($(date +%s) - TOTAL_START))

echo
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${CYAN}                   SMOKE TEST SUMMARY${NC}"
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo
echo -e "Mode:            $([ "$USE_MOCK" == "true" ] && echo 'Mock tools' || echo 'Real tools')"
echo -e "Profile:         $QUALITY_PROFILE"
echo -e "Total time:      ${TOTAL_TIME}s"
echo -e "Warnings:        $VALIDATION_WARNINGS"
echo

if [[ "$KEEP_ARTIFACTS" == "true" ]]; then
    echo -e "Artifacts:       $SMOKE_OUTPUT_DIR"
    echo
    echo "Output structure:"
    find "$SMOKE_OUTPUT_DIR" -type f -exec ls -lh {} \; 2>/dev/null | head -20 || true
    echo
fi

if [[ "$VALIDATION_PASSED" == "true" ]]; then
    echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${GREEN}            ✓ ALL SMOKE TESTS PASSED ✓${NC}"
    echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    exit $EX_OK
else
    echo -e "${RED}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${RED}            ✗ SMOKE TESTS FAILED ✗${NC}"
    echo -e "${RED}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    exit $EX_VALIDATE_FAIL
fi
