#!/bin/bash

# Docker Helper Script for WGS Pipeline
# Simplifies building, running, and managing the WGS pipeline container

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Configuration
IMAGE_NAME="wgs-pipeline"
CONTAINER_NAME="wgs-analysis"
DEFAULT_DATA_DIR="$(pwd)/data"
DEFAULT_RESULTS_DIR="$(pwd)/results"

# Logging functions
log() {
    echo -e "${GREEN}[DOCKER]${NC} $1"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

warning() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

# Help function
show_help() {
    cat << EOF
WGS Pipeline Docker Helper

USAGE:
    $0 <command> [options]

COMMANDS:
    build           Build the Docker image
    run             Run the container interactively
    exec            Execute command in running container
    stop            Stop the running container
    clean           Remove container and image
    status          Show container status
    logs            Show container logs
    test            Run tests in container
    shell           Open shell in running container
    
BUILD OPTIONS:
    --no-cache      Build without using cache
    --pull          Pull latest base image
    
RUN OPTIONS:
    --data-dir      Data directory to mount (default: ./data)
    --results-dir   Results directory to mount (default: ./results)
    --memory        Memory limit (default: 16g)
    --cpus          CPU limit (default: 4)
    --detach        Run in background
    
EXAMPLES:
    # Build the image
    $0 build

    # Run interactively  
    $0 run

    # Run with custom data directory
    $0 run --data-dir /path/to/my/data

    # Execute a command in the container
    $0 exec "check_requirements.sh"

    # Run tests
    $0 test

    # Open shell in running container
    $0 shell

    # Clean up everything
    $0 clean
EOF
}

# Build Docker image
build_image() {
    local build_args=()
    local no_cache=false
    local pull=false
    
    # Parse build options
    while [[ $# -gt 0 ]]; do
        case $1 in
            --no-cache)
                no_cache=true
                shift
                ;;
            --pull)
                pull=true
                shift
                ;;
            *)
                break
                ;;
        esac
    done
    
    log "Building WGS Pipeline Docker image..."
    
    # Add build arguments
    build_args+=(--build-arg BUILD_DATE="$(date -u +'%Y-%m-%dT%H:%M:%SZ')")
    build_args+=(--build-arg BUILD_VERSION="latest")
    
    if command -v git &> /dev/null && git rev-parse --git-dir &> /dev/null; then
        build_args+=(--build-arg VCS_REF="$(git rev-parse --short HEAD)")
    fi
    
    # Add optional flags
    if [[ "$no_cache" == "true" ]]; then
        build_args+=(--no-cache)
    fi
    
    if [[ "$pull" == "true" ]]; then
        build_args+=(--pull)
    fi
    
    # Build the image
    docker build -t "$IMAGE_NAME:latest" "${build_args[@]}" .
    
    log "✅ Docker image built successfully: $IMAGE_NAME:latest"
}

# Run container
run_container() {
    local data_dir="$DEFAULT_DATA_DIR"
    local results_dir="$DEFAULT_RESULTS_DIR"
    local memory="16g"
    local cpus="4"
    local detach=false
    local run_args=()
    
    # Parse run options
    while [[ $# -gt 0 ]]; do
        case $1 in
            --data-dir)
                data_dir="$2"
                shift 2
                ;;
            --results-dir)
                results_dir="$2"
                shift 2
                ;;
            --memory)
                memory="$2"
                shift 2
                ;;
            --cpus)
                cpus="$2"
                shift 2
                ;;
            --detach)
                detach=true
                shift
                ;;
            *)
                break
                ;;
        esac
    done
    
    # Stop existing container if running
    if docker ps -q -f name="$CONTAINER_NAME" | grep -q .; then
        warning "Stopping existing container..."
        docker stop "$CONTAINER_NAME" > /dev/null
    fi
    
    # Remove existing container if exists
    if docker ps -aq -f name="$CONTAINER_NAME" | grep -q .; then
        info "Removing existing container..."
        docker rm "$CONTAINER_NAME" > /dev/null
    fi
    
    # Ensure directories exist
    mkdir -p "$data_dir" "$results_dir"
    
    # Prepare docker run arguments
    run_args+=(--name "$CONTAINER_NAME")
    run_args+=(--memory "$memory")
    run_args+=(--cpus "$cpus")
    run_args+=(-v "$data_dir:/opt/wgs-pipeline/data")
    run_args+=(-v "$results_dir:/opt/wgs-pipeline/results")
    run_args+=(-v "$(pwd)/logs:/opt/wgs-pipeline/logs")
    
    if [[ "$detach" == "true" ]]; then
        run_args+=(-d)
        log "Starting container in background..."
    else
        run_args+=(-it)
        log "Starting interactive container..."
    fi
    
    # Run the container
    docker run "${run_args[@]}" "$IMAGE_NAME:latest"
    
    if [[ "$detach" == "true" ]]; then
        log "✅ Container started in background: $CONTAINER_NAME"
        log "Use '$0 shell' to open a shell or '$0 exec <command>' to run commands"
    fi
}

# Execute command in container
exec_command() {
    if ! docker ps -q -f name="$CONTAINER_NAME" | grep -q .; then
        error "Container '$CONTAINER_NAME' is not running"
        exit 1
    fi
    
    local command="$*"
    if [[ -z "$command" ]]; then
        command="/bin/bash"
    fi
    
    log "Executing: $command"
    docker exec -it "$CONTAINER_NAME" bash -c "source activate wgs_analysis && $command"
}

# Open shell in container
open_shell() {
    if ! docker ps -q -f name="$CONTAINER_NAME" | grep -q .; then
        error "Container '$CONTAINER_NAME' is not running"
        exit 1
    fi
    
    log "Opening shell in container..."
    docker exec -it "$CONTAINER_NAME" bash -c "source activate wgs_analysis && /bin/bash"
}

# Stop container
stop_container() {
    if docker ps -q -f name="$CONTAINER_NAME" | grep -q .; then
        log "Stopping container..."
        docker stop "$CONTAINER_NAME"
        log "✅ Container stopped"
    else
        warning "Container is not running"
    fi
}

# Show container status
show_status() {
    log "Container status:"
    
    if docker ps -q -f name="$CONTAINER_NAME" | grep -q .; then
        echo -e "${GREEN}Running${NC}"
        docker ps -f name="$CONTAINER_NAME" --format "table {{.Names}}\t{{.Status}}\t{{.Ports}}"
    elif docker ps -aq -f name="$CONTAINER_NAME" | grep -q .; then
        echo -e "${YELLOW}Stopped${NC}"
        docker ps -a -f name="$CONTAINER_NAME" --format "table {{.Names}}\t{{.Status}}"
    else
        echo -e "${RED}Not found${NC}"
    fi
    
    echo
    log "Image status:"
    if docker images -q "$IMAGE_NAME" | grep -q .; then
        docker images "$IMAGE_NAME" --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}\t{{.CreatedSince}}"
    else
        echo "Image not found"
    fi
}

# Show logs
show_logs() {
    if docker ps -aq -f name="$CONTAINER_NAME" | grep -q .; then
        docker logs "$CONTAINER_NAME" "$@"
    else
        error "Container '$CONTAINER_NAME' not found"
        exit 1
    fi
}

# Run tests
run_tests() {
    log "Running tests in container..."
    
    if ! docker ps -q -f name="$CONTAINER_NAME" | grep -q .; then
        log "Starting container for testing..."
        run_container --detach
        sleep 5
    fi
    
    exec_command "cd /opt/wgs-pipeline && tests/run_tests.sh --verbose"
}

# Clean up
clean_up() {
    log "Cleaning up Docker resources..."
    
    # Stop and remove container
    if docker ps -aq -f name="$CONTAINER_NAME" | grep -q .; then
        docker stop "$CONTAINER_NAME" 2>/dev/null || true
        docker rm "$CONTAINER_NAME"
        log "✅ Container removed"
    fi
    
    # Remove image
    if docker images -q "$IMAGE_NAME" | grep -q .; then
        docker rmi "$IMAGE_NAME:latest"
        log "✅ Image removed"
    fi
    
    # Clean up dangling images
    if docker images -q -f dangling=true | grep -q .; then
        docker image prune -f
        log "✅ Dangling images cleaned"
    fi
}

# Main command dispatcher
main() {
    if [[ $# -eq 0 ]]; then
        show_help
        exit 0
    fi
    
    local command="$1"
    shift
    
    case "$command" in
        build)
            build_image "$@"
            ;;
        run)
            run_container "$@"
            ;;
        exec)
            exec_command "$@"
            ;;
        shell)
            open_shell
            ;;
        stop)
            stop_container
            ;;
        status)
            show_status
            ;;
        logs)
            show_logs "$@"
            ;;
        test)
            run_tests
            ;;
        clean)
            clean_up
            ;;
        help|--help|-h)
            show_help
            ;;
        *)
            error "Unknown command: $command"
            echo "Use '$0 help' for usage information"
            exit 1
            ;;
    esac
}

# Run main function
main "$@"