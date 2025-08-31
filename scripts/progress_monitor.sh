#!/bin/bash

# Progress Monitoring System for WGS Pipeline
# Provides real-time progress tracking, resource monitoring, and notifications

# Colors and formatting
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
PURPLE='\033[0;35m'
NC='\033[0m'
BOLD='\033[1m'

# Progress tracking variables
declare -A STEP_STATUS
declare -A STEP_START_TIME
declare -A STEP_DURATION
PIPELINE_START_TIME=""
CURRENT_STEP=""
TOTAL_STEPS=0
COMPLETED_STEPS=0

# Resource monitoring
MONITOR_INTERVAL=30  # seconds
MONITOR_PID=""
LOG_FILE=""
RESOURCE_LOG=""

# Notification settings
ENABLE_NOTIFICATIONS=false
EMAIL_ADDRESS=""
SLACK_WEBHOOK=""

# Initialize progress monitoring
init_progress_monitor() {
    local pipeline_name="$1"
    local total_steps="$2"
    local log_dir="${3:-logs}"
    
    PIPELINE_NAME="$pipeline_name"
    TOTAL_STEPS="$total_steps"
    PIPELINE_START_TIME=$(date +%s)
    
    # Create log files
    mkdir -p "$log_dir"
    LOG_FILE="$log_dir/${pipeline_name}_progress.log"
    RESOURCE_LOG="$log_dir/${pipeline_name}_resources.log"
    
    # Initialize log files
    echo "# WGS Pipeline Progress Log: $pipeline_name" > "$LOG_FILE"
    echo "# Started: $(date)" >> "$LOG_FILE"
    echo "# Total steps: $total_steps" >> "$LOG_FILE"
    echo "" >> "$LOG_FILE"
    
    # Resource log header
    echo "# Resource Usage Log: $pipeline_name" > "$RESOURCE_LOG"
    echo "# Timestamp,CPU%,Memory_MB,Memory%,Disk_Available_GB,Load_Avg" >> "$RESOURCE_LOG"
    
    # Start resource monitoring
    start_resource_monitoring
    
    # Display initial progress
    show_pipeline_header
    log_message "INFO" "Progress monitoring initialized for $pipeline_name"
}

# Show pipeline header
show_pipeline_header() {
    local width=80
    local title="WGS PIPELINE: $PIPELINE_NAME"
    local padding=$(( (width - ${#title}) / 2 ))
    
    echo
    printf "${BOLD}${BLUE}%*s${NC}\n" $width | tr ' ' '='
    printf "${BOLD}${BLUE}%*s%s%*s${NC}\n" $padding "" "$title" $padding ""
    printf "${BOLD}${BLUE}%*s${NC}\n" $width | tr ' ' '='
    echo
    
    printf "${CYAN}Started:${NC} $(date)\n"
    printf "${CYAN}Total Steps:${NC} $TOTAL_STEPS\n"
    printf "${CYAN}Log File:${NC} $LOG_FILE\n"
    echo
}

# Start a new step
start_step() {
    local step_name="$1"
    local step_description="$2"
    
    CURRENT_STEP="$step_name"
    STEP_START_TIME["$step_name"]=$(date +%s)
    STEP_STATUS["$step_name"]="RUNNING"
    
    # Log step start
    log_message "STEP_START" "$step_name: $step_description"
    
    # Display progress
    show_step_start "$step_name" "$step_description"
}

# Complete a step
complete_step() {
    local step_name="$1"
    local success="${2:-true}"
    
    local end_time=$(date +%s)
    local start_time="${STEP_START_TIME[$step_name]}"
    local duration=$((end_time - start_time))
    
    STEP_DURATION["$step_name"]=$duration
    
    if [[ "$success" == "true" ]]; then
        STEP_STATUS["$step_name"]="COMPLETED"
        ((COMPLETED_STEPS++))
        log_message "STEP_COMPLETE" "$step_name completed in $(format_duration $duration)"
        show_step_complete "$step_name" "$duration"
    else
        STEP_STATUS["$step_name"]="FAILED"
        log_message "STEP_FAILED" "$step_name failed after $(format_duration $duration)"
        show_step_failed "$step_name" "$duration"
    fi
    
    # Show overall progress
    show_overall_progress
}

# Display step start
show_step_start() {
    local step_name="$1"
    local description="$2"
    
    local step_num=$((COMPLETED_STEPS + 1))
    echo
    printf "${BOLD}${YELLOW}[STEP $step_num/$TOTAL_STEPS]${NC} ${BOLD}$step_name${NC}\n"
    printf "${CYAN}Description:${NC} $description\n"
    printf "${CYAN}Started:${NC} $(date '+%Y-%m-%d %H:%M:%S')\n"
    
    # Show progress bar
    show_progress_bar $COMPLETED_STEPS $TOTAL_STEPS "Steps"
    echo
}

# Display step completion
show_step_complete() {
    local step_name="$1"
    local duration="$2"
    
    printf "${GREEN}✓${NC} ${BOLD}$step_name${NC} completed in ${GREEN}$(format_duration $duration)${NC}\n"
}

# Display step failure
show_step_failed() {
    local step_name="$1"
    local duration="$2"
    
    printf "${RED}✗${NC} ${BOLD}$step_name${NC} failed after ${RED}$(format_duration $duration)${NC}\n"
}

# Show overall progress
show_overall_progress() {
    echo
    printf "${BOLD}${PURPLE}Overall Progress:${NC}\n"
    show_progress_bar $COMPLETED_STEPS $TOTAL_STEPS "Pipeline"
    
    local elapsed=$(($(date +%s) - PIPELINE_START_TIME))
    printf "${CYAN}Elapsed Time:${NC} $(format_duration $elapsed)\n"
    
    if [[ $COMPLETED_STEPS -gt 0 && $COMPLETED_STEPS -lt $TOTAL_STEPS ]]; then
        local avg_step_time=$((elapsed / COMPLETED_STEPS))
        local remaining_time=$(((TOTAL_STEPS - COMPLETED_STEPS) * avg_step_time))
        printf "${CYAN}Estimated Remaining:${NC} $(format_duration $remaining_time)\n"
        printf "${CYAN}Estimated Completion:${NC} $(date -d "+${remaining_time} seconds" '+%Y-%m-%d %H:%M:%S')\n"
    fi
    echo
}

# Show progress bar
show_progress_bar() {
    local current="$1"
    local total="$2"
    local label="$3"
    local width=50
    
    local percentage=$((current * 100 / total))
    local filled=$((current * width / total))
    local empty=$((width - filled))
    
    # Choose color based on progress
    local color="$RED"
    if [[ $percentage -ge 75 ]]; then
        color="$GREEN"
    elif [[ $percentage -ge 50 ]]; then
        color="$YELLOW"
    elif [[ $percentage -ge 25 ]]; then
        color="$CYAN"
    fi
    
    printf "${BOLD}$label Progress:${NC} ["
    printf "${color}%*s${NC}" $filled | tr ' ' '█'
    printf "%*s" $empty | tr ' ' '░'
    printf "] ${BOLD}$current/$total${NC} (${color}$percentage%%${NC})\n"
}

# Start resource monitoring
start_resource_monitoring() {
    if [[ "$MONITOR_PID" != "" ]]; then
        return  # Already monitoring
    fi
    
    {
        while true; do
            local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
            local cpu_usage=$(top -bn1 | grep "Cpu(s)" | awk '{print $2}' | sed 's/%us,//' || echo "0")
            local memory_info=$(free -m | awk 'NR==2{printf "%.1f,%.1f", $3,$2}' || echo "0,0")
            local disk_available=$(df -BG . | awk 'NR==2 {print $4}' | sed 's/G//' || echo "0")
            local load_avg=$(uptime | awk -F'load average:' '{print $2}' | awk '{print $1}' | sed 's/,//' || echo "0")
            
            echo "$timestamp,$cpu_usage,$memory_info,$disk_available,$load_avg" >> "$RESOURCE_LOG"
            sleep $MONITOR_INTERVAL
        done
    } &
    
    MONITOR_PID=$!
    log_message "INFO" "Resource monitoring started (PID: $MONITOR_PID)"
}

# Stop resource monitoring
stop_resource_monitoring() {
    if [[ "$MONITOR_PID" != "" ]]; then
        kill $MONITOR_PID 2>/dev/null || true
        wait $MONITOR_PID 2>/dev/null || true
        MONITOR_PID=""
        log_message "INFO" "Resource monitoring stopped"
    fi
}

# Log message to file
log_message() {
    local level="$1"
    local message="$2"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    echo "[$timestamp] [$level] $message" >> "$LOG_FILE"
}

# Format duration in human readable format
format_duration() {
    local duration="$1"
    
    if [[ $duration -lt 60 ]]; then
        echo "${duration}s"
    elif [[ $duration -lt 3600 ]]; then
        local minutes=$((duration / 60))
        local seconds=$((duration % 60))
        echo "${minutes}m ${seconds}s"
    else
        local hours=$((duration / 3600))
        local minutes=$(((duration % 3600) / 60))
        local seconds=$((duration % 60))
        echo "${hours}h ${minutes}m ${seconds}s"
    fi
}

# Send notification
send_notification() {
    local title="$1"
    local message="$2"
    local status="$3"  # success, failure, info
    
    if [[ "$ENABLE_NOTIFICATIONS" != "true" ]]; then
        return
    fi
    
    # Email notification
    if [[ -n "$EMAIL_ADDRESS" ]] && command -v mail &> /dev/null; then
        echo "$message" | mail -s "WGS Pipeline: $title" "$EMAIL_ADDRESS"
    fi
    
    # Slack notification
    if [[ -n "$SLACK_WEBHOOK" ]] && command -v curl &> /dev/null; then
        local color="good"
        if [[ "$status" == "failure" ]]; then
            color="danger"
        elif [[ "$status" == "info" ]]; then
            color="warning"
        fi
        
        local payload=$(cat << EOF
{
    "attachments": [
        {
            "color": "$color",
            "title": "WGS Pipeline: $title",
            "text": "$message",
            "footer": "WGS Pipeline Monitor",
            "ts": $(date +%s)
        }
    ]
}
EOF
)
        curl -X POST -H 'Content-type: application/json' --data "$payload" "$SLACK_WEBHOOK" 2>/dev/null || true
    fi
    
    # Desktop notification (if available)
    if command -v notify-send &> /dev/null; then
        notify-send "WGS Pipeline: $title" "$message"
    fi
}

# Finish pipeline monitoring
finish_progress_monitor() {
    local success="${1:-true}"
    
    stop_resource_monitoring
    
    local end_time=$(date +%s)
    local total_duration=$((end_time - PIPELINE_START_TIME))
    
    # Final progress display
    echo
    if [[ "$success" == "true" ]]; then
        printf "${BOLD}${GREEN}🎉 PIPELINE COMPLETED SUCCESSFULLY!${NC}\n"
        send_notification "Completed" "Pipeline '$PIPELINE_NAME' completed successfully in $(format_duration $total_duration)" "success"
    else
        printf "${BOLD}${RED}❌ PIPELINE FAILED!${NC}\n"
        send_notification "Failed" "Pipeline '$PIPELINE_NAME' failed after $(format_duration $total_duration)" "failure"
    fi
    
    echo
    printf "${CYAN}Total Duration:${NC} $(format_duration $total_duration)\n"
    printf "${CYAN}Steps Completed:${NC} $COMPLETED_STEPS/$TOTAL_STEPS\n"
    printf "${CYAN}Log File:${NC} $LOG_FILE\n"
    printf "${CYAN}Resource Log:${NC} $RESOURCE_LOG\n"
    
    # Log final status
    if [[ "$success" == "true" ]]; then
        log_message "PIPELINE_COMPLETE" "Pipeline completed successfully in $(format_duration $total_duration)"
    else
        log_message "PIPELINE_FAILED" "Pipeline failed after $(format_duration $total_duration)"
    fi
    
    # Show step summary
    echo
    printf "${BOLD}${PURPLE}Step Summary:${NC}\n"
    printf "${BOLD}%-30s %-12s %-10s${NC}\n" "Step Name" "Status" "Duration"
    printf "%*s\n" 52 | tr ' ' '-'
    
    for step in "${!STEP_STATUS[@]}"; do
        local status="${STEP_STATUS[$step]}"
        local duration="${STEP_DURATION[$step]:-0}"
        local status_color="$YELLOW"
        
        case "$status" in
            "COMPLETED") status_color="$GREEN" ;;
            "FAILED") status_color="$RED" ;;
        esac
        
        printf "%-30s ${status_color}%-12s${NC} %-10s\n" \
            "$step" "$status" "$(format_duration $duration)"
    done
    
    echo
}

# Enable notifications
enable_notifications() {
    ENABLE_NOTIFICATIONS=true
    EMAIL_ADDRESS="$1"
    SLACK_WEBHOOK="$2"
    log_message "INFO" "Notifications enabled"
}

# Resource usage summary
show_resource_summary() {
    if [[ ! -f "$RESOURCE_LOG" ]]; then
        return
    fi
    
    echo
    printf "${BOLD}${PURPLE}Resource Usage Summary:${NC}\n"
    
    # Skip header line and calculate averages
    local avg_cpu=$(tail -n +2 "$RESOURCE_LOG" | awk -F',' '{sum+=$2; count++} END {if(count>0) print sum/count; else print 0}')
    local avg_memory_mb=$(tail -n +2 "$RESOURCE_LOG" | awk -F',' '{sum+=$3; count++} END {if(count>0) print sum/count; else print 0}')
    local min_disk=$(tail -n +2 "$RESOURCE_LOG" | awk -F',' '{min=999999} {if($5<min) min=$5} END {print min}')
    local max_load=$(tail -n +2 "$RESOURCE_LOG" | awk -F',' '{max=0} {if($6>max) max=$6} END {print max}')
    
    printf "${CYAN}Average CPU Usage:${NC} %.1f%%\n" "$avg_cpu"
    printf "${CYAN}Average Memory Usage:${NC} %.0f MB\n" "$avg_memory_mb"
    printf "${CYAN}Minimum Disk Available:${NC} %.0f GB\n" "$min_disk"
    printf "${CYAN}Maximum Load Average:${NC} %.2f\n" "$max_load"
}

# Trap handler for cleanup
cleanup_progress_monitor() {
    stop_resource_monitoring
    log_message "INFO" "Progress monitoring interrupted"
}

# Set up trap for cleanup
trap cleanup_progress_monitor INT TERM EXIT

# Export functions for use by other scripts
export -f init_progress_monitor
export -f start_step
export -f complete_step
export -f finish_progress_monitor
export -f enable_notifications
export -f show_resource_summary