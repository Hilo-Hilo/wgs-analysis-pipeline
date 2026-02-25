#!/bin/bash

# Hardware Detection Helper for WGS Pipeline
# Auto-detects CPU cores, RAM, and GPU availability for optimal resource allocation
# Source this script or call with subcommands: detect_hardware.sh [cpu|ram|gpu|mode|all|json]

set -e

# Script info
SCRIPT_VERSION="1.0.0"
SCRIPT_NAME="detect_hardware.sh"

# =============================================================================
# HARDWARE DETECTION FUNCTIONS
# =============================================================================

# Detect available CPU cores (physical or logical depending on platform)
detect_cpu_cores() {
    local cores=0
    
    if [[ -f /proc/cpuinfo ]]; then
        # Linux: count processor entries
        cores=$(grep -c ^processor /proc/cpuinfo 2>/dev/null || echo "0")
    elif command -v sysctl &>/dev/null; then
        # macOS/BSD: use sysctl
        cores=$(sysctl -n hw.ncpu 2>/dev/null || echo "0")
    elif command -v nproc &>/dev/null; then
        # Fallback: nproc (coreutils)
        cores=$(nproc 2>/dev/null || echo "0")
    fi
    
    # Sanity check: at least 1 core
    if [[ "$cores" -lt 1 ]]; then
        cores=1
    fi
    
    echo "$cores"
}

# Detect total RAM in GB (rounded down)
detect_ram_gb() {
    local ram_kb=0
    local ram_gb=0
    
    if [[ -f /proc/meminfo ]]; then
        # Linux: read MemTotal from /proc/meminfo
        ram_kb=$(grep MemTotal /proc/meminfo 2>/dev/null | awk '{print $2}' || echo "0")
        ram_gb=$((ram_kb / 1024 / 1024))
    elif command -v sysctl &>/dev/null; then
        # macOS: use sysctl hw.memsize (bytes)
        local ram_bytes
        ram_bytes=$(sysctl -n hw.memsize 2>/dev/null || echo "0")
        ram_gb=$((ram_bytes / 1024 / 1024 / 1024))
    fi
    
    # Sanity check: at least 1 GB
    if [[ "$ram_gb" -lt 1 ]]; then
        ram_gb=1
    fi
    
    echo "$ram_gb"
}

# Detect available GPUs via nvidia-smi
# Returns: count of NVIDIA GPUs (0 if none or nvidia-smi unavailable)
detect_gpu_count() {
    local gpu_count=0
    
    if command -v nvidia-smi &>/dev/null; then
        # Query GPU count (works in containers with proper mounts)
        gpu_count=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | wc -l | tr -d ' ')
        
        # Validate it's a number
        if ! [[ "$gpu_count" =~ ^[0-9]+$ ]]; then
            gpu_count=0
        fi
    fi
    
    echo "$gpu_count"
}

# Get GPU model names (comma-separated)
detect_gpu_models() {
    local models=""
    
    if command -v nvidia-smi &>/dev/null; then
        models=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | paste -sd',' - || echo "")
    fi
    
    echo "$models"
}

# Get total GPU memory in GB (sum across all GPUs)
detect_gpu_memory_gb() {
    local total_mb=0
    
    if command -v nvidia-smi &>/dev/null; then
        total_mb=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits 2>/dev/null | awk '{s+=$1} END {print s}' || echo "0")
    fi
    
    # Convert to GB
    local total_gb=$((total_mb / 1024))
    echo "$total_gb"
}

# =============================================================================
# PERFORMANCE MODE SELECTION
# =============================================================================

# Select optimal performance mode based on detected hardware
# Modes: laptop, workstation, server, dgx
# Returns mode name and recommended settings
select_performance_mode() {
    local cores
    local ram_gb
    local gpu_count
    
    cores=$(detect_cpu_cores)
    ram_gb=$(detect_ram_gb)
    gpu_count=$(detect_gpu_count)
    
    local mode="laptop"  # Conservative default
    
    # DGX-class: 8+ GPUs or 4+ GPUs with 128GB+ RAM
    if [[ "$gpu_count" -ge 8 ]] || { [[ "$gpu_count" -ge 4 ]] && [[ "$ram_gb" -ge 128 ]]; }; then
        mode="dgx"
    # Server: 64GB+ RAM and 16+ cores (with or without GPU)
    elif [[ "$ram_gb" -ge 64 ]] && [[ "$cores" -ge 16 ]]; then
        mode="server"
    # Workstation: 32GB+ RAM and 8+ cores, or 1-3 GPUs
    elif { [[ "$ram_gb" -ge 32 ]] && [[ "$cores" -ge 8 ]]; } || { [[ "$gpu_count" -ge 1 ]] && [[ "$gpu_count" -le 3 ]]; }; then
        mode="workstation"
    # Laptop: Everything else (conservative defaults)
    else
        mode="laptop"
    fi
    
    echo "$mode"
}

# Get recommended thread count for detected hardware
recommend_threads() {
    local cores
    local ram_gb
    local mode
    
    cores=$(detect_cpu_cores)
    ram_gb=$(detect_ram_gb)
    mode=$(select_performance_mode)
    
    local threads=4  # Safe default
    
    case "$mode" in
        dgx)
            # DGX: Use most cores but leave some for system
            threads=$((cores > 8 ? cores - 4 : cores))
            ;;
        server)
            # Server: Use up to 75% of cores
            threads=$((cores * 3 / 4))
            # Cap based on RAM (4GB per thread for WGS)
            local ram_threads=$((ram_gb / 4))
            if [[ "$ram_threads" -lt "$threads" ]]; then
                threads=$ram_threads
            fi
            ;;
        workstation)
            # Workstation: Use up to 50% of cores
            threads=$((cores / 2))
            # Cap based on RAM (4GB per thread)
            local ram_threads=$((ram_gb / 4))
            if [[ "$ram_threads" -lt "$threads" ]]; then
                threads=$ram_threads
            fi
            # At least 4 threads if possible
            if [[ "$threads" -lt 4 ]] && [[ "$cores" -ge 4 ]]; then
                threads=4
            fi
            ;;
        laptop)
            # Laptop: Conservative (2-4 threads)
            threads=2
            if [[ "$cores" -ge 4 ]] && [[ "$ram_gb" -ge 16 ]]; then
                threads=4
            fi
            ;;
    esac
    
    # Minimum 1 thread
    if [[ "$threads" -lt 1 ]]; then
        threads=1
    fi
    
    echo "$threads"
}

# Recommend GPU usage based on hardware
recommend_gpu_settings() {
    local gpu_count
    local mode
    
    gpu_count=$(detect_gpu_count)
    mode=$(select_performance_mode)
    
    local use_gpu="false"
    local gpu_aligner="parabricks"
    local recommended_gpus=0
    
    if [[ "$gpu_count" -ge 1 ]]; then
        # Check if Parabricks is available
        if command -v pbrun &>/dev/null; then
            use_gpu="true"
            
            case "$mode" in
                dgx)
                    # DGX: Use up to 4 GPUs for alignment
                    recommended_gpus=$((gpu_count > 4 ? 4 : gpu_count))
                    ;;
                server|workstation)
                    # Server/Workstation: Use available GPUs (1-2)
                    recommended_gpus=$((gpu_count > 2 ? 2 : gpu_count))
                    ;;
                *)
                    recommended_gpus=1
                    ;;
            esac
        fi
    fi
    
    echo "$use_gpu $gpu_aligner $recommended_gpus"
}

# =============================================================================
# OUTPUT FORMATTING
# =============================================================================

# Print all hardware info as JSON
print_json() {
    local cores
    local ram_gb
    local gpu_count
    local gpu_models
    local gpu_memory
    local mode
    local threads
    local gpu_settings
    
    cores=$(detect_cpu_cores)
    ram_gb=$(detect_ram_gb)
    gpu_count=$(detect_gpu_count)
    gpu_models=$(detect_gpu_models)
    gpu_memory=$(detect_gpu_memory_gb)
    mode=$(select_performance_mode)
    threads=$(recommend_threads)
    
    read -r use_gpu gpu_aligner recommended_gpus <<< "$(recommend_gpu_settings)"
    
    cat << EOF
{
  "hardware": {
    "cpu_cores": $cores,
    "ram_gb": $ram_gb,
    "gpu_count": $gpu_count,
    "gpu_models": "$gpu_models",
    "gpu_memory_gb": $gpu_memory
  },
  "recommendations": {
    "mode": "$mode",
    "threads": $threads,
    "use_gpu": $use_gpu,
    "gpu_aligner": "$gpu_aligner",
    "gpu_count": $recommended_gpus
  }
}
EOF
}

# Print human-readable summary
print_summary() {
    local cores
    local ram_gb
    local gpu_count
    local gpu_models
    local mode
    local threads
    
    cores=$(detect_cpu_cores)
    ram_gb=$(detect_ram_gb)
    gpu_count=$(detect_gpu_count)
    gpu_models=$(detect_gpu_models)
    mode=$(select_performance_mode)
    threads=$(recommend_threads)
    
    read -r use_gpu gpu_aligner recommended_gpus <<< "$(recommend_gpu_settings)"
    
    echo "=== Hardware Detection Summary ==="
    echo "CPU Cores:    $cores"
    echo "RAM:          ${ram_gb}GB"
    echo "GPUs:         $gpu_count"
    if [[ -n "$gpu_models" ]]; then
        echo "GPU Models:   $gpu_models"
    fi
    echo ""
    echo "=== Recommendations ==="
    echo "Mode:         $mode"
    echo "Threads:      $threads"
    echo "Use GPU:      $use_gpu"
    if [[ "$use_gpu" == "true" ]]; then
        echo "GPU Aligner:  $gpu_aligner"
        echo "GPU Count:    $recommended_gpus"
    fi
}

# =============================================================================
# CLI INTERFACE
# =============================================================================

show_help() {
    cat << EOF
$SCRIPT_NAME v$SCRIPT_VERSION - Hardware Detection for WGS Pipeline

DESCRIPTION:
    Detects available hardware (CPU, RAM, GPU) and recommends optimal
    performance settings for the WGS analysis pipeline.

USAGE:
    $0 [COMMAND]

COMMANDS:
    cpu         Print detected CPU core count
    ram         Print detected RAM in GB
    gpu         Print detected GPU count
    gpu-models  Print GPU model names
    gpu-memory  Print total GPU memory in GB
    mode        Print recommended performance mode
    threads     Print recommended thread count
    gpu-settings Print recommended GPU settings (use_gpu aligner count)
    all         Print human-readable summary (default)
    json        Print all info as JSON

SOURCING:
    Source this script to use detection functions in other scripts:
    
    source scripts/detect_hardware.sh
    THREADS=\$(recommend_threads)
    MODE=\$(select_performance_mode)

PERFORMANCE MODES:
    laptop      Conservative (8-16GB RAM, 2-4 cores)
    workstation Balanced (32GB+ RAM, 8+ cores, or 1-3 GPUs)
    server      High-performance (64GB+ RAM, 16+ cores)
    dgx         Maximum (8+ GPUs or 4+ GPUs with 128GB+ RAM)

EXAMPLES:
    # Get recommended thread count
    $0 threads

    # Get full summary as JSON
    $0 json

    # Use in scripts
    source scripts/detect_hardware.sh
    if [[ \$(detect_gpu_count) -ge 1 ]]; then
        echo "GPU available!"
    fi

EOF
}

# Main CLI handler
main() {
    local cmd="${1:-all}"
    
    case "$cmd" in
        -h|--help|help)
            show_help
            ;;
        --version)
            echo "$SCRIPT_NAME v$SCRIPT_VERSION"
            ;;
        cpu)
            detect_cpu_cores
            ;;
        ram)
            detect_ram_gb
            ;;
        gpu)
            detect_gpu_count
            ;;
        gpu-models)
            detect_gpu_models
            ;;
        gpu-memory)
            detect_gpu_memory_gb
            ;;
        mode)
            select_performance_mode
            ;;
        threads)
            recommend_threads
            ;;
        gpu-settings)
            recommend_gpu_settings
            ;;
        json)
            print_json
            ;;
        all|summary)
            print_summary
            ;;
        *)
            echo "Unknown command: $cmd" >&2
            echo "Use --help for usage information" >&2
            exit 1
            ;;
    esac
}

# Run CLI if executed directly (not sourced)
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
