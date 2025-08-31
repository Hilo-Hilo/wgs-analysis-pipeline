#!/usr/bin/env python3
"""
Real-time Resource Dashboard for WGS Pipeline
Displays live system resource usage during pipeline execution
"""

import sys
import time
import argparse
import subprocess
import shutil
from datetime import datetime
from pathlib import Path

def clear_screen():
    """Clear terminal screen"""
    subprocess.run(['clear'], check=False)

def get_terminal_size():
    """Get terminal width and height"""
    size = shutil.get_terminal_size(fallback=(80, 24))
    return size.columns, size.lines

def get_system_info():
    """Get current system resource usage"""
    try:
        # CPU usage
        cpu_result = subprocess.run(
            ['top', '-bn1'], 
            capture_output=True, 
            text=True, 
            timeout=5
        )
        cpu_line = [line for line in cpu_result.stdout.split('\n') if 'Cpu(s)' in line]
        cpu_usage = 0.0
        if cpu_line:
            # Parse CPU usage from top output
            parts = cpu_line[0].split(',')
            for part in parts:
                if '%us' in part or 'us' in part:
                    cpu_usage = float(part.strip().split('%')[0].split()[-1])
                    break
    except:
        cpu_usage = 0.0
    
    try:
        # Memory usage
        mem_result = subprocess.run(
            ['free', '-m'], 
            capture_output=True, 
            text=True, 
            timeout=5
        )
        mem_lines = mem_result.stdout.split('\n')
        for line in mem_lines:
            if line.startswith('Mem:'):
                parts = line.split()
                total_mem = int(parts[1])
                used_mem = int(parts[2])
                mem_percent = (used_mem / total_mem) * 100
                break
        else:
            total_mem, used_mem, mem_percent = 0, 0, 0.0
    except:
        total_mem, used_mem, mem_percent = 0, 0, 0.0
    
    try:
        # Disk usage
        disk_result = subprocess.run(
            ['df', '-h', '.'], 
            capture_output=True, 
            text=True, 
            timeout=5
        )
        disk_lines = disk_result.stdout.split('\n')
        if len(disk_lines) > 1:
            parts = disk_lines[1].split()
            disk_total = parts[1]
            disk_used = parts[2]
            disk_avail = parts[3]
            disk_percent = int(parts[4].rstrip('%'))
        else:
            disk_total, disk_used, disk_avail, disk_percent = "0G", "0G", "0G", 0
    except:
        disk_total, disk_used, disk_avail, disk_percent = "0G", "0G", "0G", 0
    
    try:
        # Load average
        load_result = subprocess.run(
            ['uptime'], 
            capture_output=True, 
            text=True, 
            timeout=5
        )
        load_line = load_result.stdout.strip()
        load_avg = load_line.split('load average:')[1].split(',')[0].strip()
    except:
        load_avg = "0.00"
    
    return {
        'cpu_usage': cpu_usage,
        'mem_total': total_mem,
        'mem_used': used_mem,
        'mem_percent': mem_percent,
        'disk_total': disk_total,
        'disk_used': disk_used,
        'disk_avail': disk_avail,
        'disk_percent': disk_percent,
        'load_avg': load_avg,
        'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    }

def draw_progress_bar(percentage, width=40, filled_char='█', empty_char='░'):
    """Draw a progress bar"""
    filled = int(width * percentage / 100)
    empty = width - filled
    
    # Color based on percentage
    if percentage < 50:
        color = '\033[92m'  # Green
    elif percentage < 75:
        color = '\033[93m'  # Yellow
    elif percentage < 90:
        color = '\033[91m'  # Red
    else:
        color = '\033[95m'  # Magenta (critical)
    
    reset = '\033[0m'
    
    bar = f"{color}{filled_char * filled}{reset}{empty_char * empty}"
    return f"[{bar}] {percentage:5.1f}%"

def display_dashboard(system_info, pipeline_log=None):
    """Display the resource dashboard"""
    width, height = get_terminal_size()
    
    # Header
    print("=" * width)
    print(f"🖥️  WGS PIPELINE RESOURCE DASHBOARD".center(width))
    print(f"Last Updated: {system_info['timestamp']}".center(width))
    print("=" * width)
    print()
    
    # System Resources
    print("📊 SYSTEM RESOURCES")
    print("-" * 50)
    
    # CPU
    cpu_bar = draw_progress_bar(system_info['cpu_usage'])
    print(f"CPU Usage:    {cpu_bar}")
    
    # Memory
    mem_bar = draw_progress_bar(system_info['mem_percent'])
    print(f"Memory:       {mem_bar} ({system_info['mem_used']} / {system_info['mem_total']} MB)")
    
    # Disk
    disk_bar = draw_progress_bar(system_info['disk_percent'])
    print(f"Disk Usage:   {disk_bar} ({system_info['disk_used']} / {system_info['disk_total']})")
    
    # Load Average
    print(f"Load Average: {system_info['load_avg']}")
    print()
    
    # Pipeline Status (if log file provided)
    if pipeline_log and Path(pipeline_log).exists():
        print("🔄 PIPELINE STATUS")
        print("-" * 50)
        
        try:
            with open(pipeline_log, 'r') as f:
                lines = f.readlines()
            
            # Get recent log entries
            recent_lines = lines[-10:]  # Last 10 lines
            for line in recent_lines:
                if line.strip():
                    # Color code based on log level
                    if '[STEP_START]' in line:
                        print(f"\033[94m{line.strip()}\033[0m")  # Blue
                    elif '[STEP_COMPLETE]' in line:
                        print(f"\033[92m{line.strip()}\033[0m")  # Green
                    elif '[STEP_FAILED]' in line or '[ERROR]' in line:
                        print(f"\033[91m{line.strip()}\033[0m")  # Red
                    elif '[INFO]' in line:
                        print(f"\033[90m{line.strip()}\033[0m")  # Gray
                    else:
                        print(line.strip())
        except Exception as e:
            print(f"Error reading pipeline log: {e}")
    
    print()
    
    # System Warnings
    warnings = []
    if system_info['cpu_usage'] > 90:
        warnings.append("⚠️  HIGH CPU USAGE (>90%)")
    if system_info['mem_percent'] > 90:
        warnings.append("⚠️  HIGH MEMORY USAGE (>90%)")
    if system_info['disk_percent'] > 90:
        warnings.append("⚠️  LOW DISK SPACE (<10% free)")
    
    if warnings:
        print("🚨 SYSTEM WARNINGS")
        print("-" * 50)
        for warning in warnings:
            print(f"\033[91m{warning}\033[0m")
        print()
    
    # Footer
    print("=" * width)
    print("Press Ctrl+C to exit | Refreshing every 5 seconds".center(width))

def monitor_resources(log_file=None, refresh_rate=5):
    """Main monitoring loop"""
    try:
        while True:
            clear_screen()
            system_info = get_system_info()
            display_dashboard(system_info, log_file)
            time.sleep(refresh_rate)
    except KeyboardInterrupt:
        print("\n\n👋 Resource monitoring stopped.")
        sys.exit(0)
    except Exception as e:
        print(f"\n\n❌ Error in monitoring: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description='Real-time resource dashboard for WGS pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic monitoring
    python3 scripts/resource_dashboard.py
    
    # Monitor with pipeline log
    python3 scripts/resource_dashboard.py --log logs/WGS_Pipeline_MySample_progress.log
    
    # Custom refresh rate
    python3 scripts/resource_dashboard.py --refresh 3
        """
    )
    
    parser.add_argument(
        '--log', '-l',
        help='Pipeline log file to monitor',
        default=None
    )
    
    parser.add_argument(
        '--refresh', '-r',
        type=int,
        help='Refresh rate in seconds (default: 5)',
        default=5
    )
    
    args = parser.parse_args()
    
    print("🚀 Starting WGS Pipeline Resource Dashboard...")
    print("   Use Ctrl+C to exit")
    time.sleep(2)
    
    monitor_resources(args.log, args.refresh)

if __name__ == '__main__':
    main()