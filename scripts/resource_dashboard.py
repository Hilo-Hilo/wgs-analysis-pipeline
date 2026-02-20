#!/usr/bin/env python3
"""
Real-time Resource Dashboard for WGS Pipeline
Displays live system resource usage during pipeline execution.
Uses psutil when available, with portable subprocess fallbacks otherwise.
"""

from __future__ import annotations

import argparse
import os
import platform
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

try:
    import psutil  # type: ignore
except ImportError:  # pragma: no cover - expected in minimal environments
    psutil = None


def clear_screen() -> None:
    """Clear terminal screen (TTY-safe)."""
    if not sys.stdout.isatty() or not os.environ.get("TERM"):
        return
    cmd = "cls" if os.name == "nt" else "clear"
    subprocess.run([cmd], check=False)


def get_terminal_size() -> tuple[int, int]:
    """Get terminal width and height."""
    size = shutil.get_terminal_size(fallback=(80, 24))
    return size.columns, size.lines


def _parse_float(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def get_cpu_usage_fallback() -> float:
    """Best-effort CPU usage without psutil."""
    system = platform.system()

    try:
        if system == "Darwin":
            result = subprocess.run(["top", "-l", "1"], capture_output=True, text=True, timeout=5, check=False)
            for line in result.stdout.splitlines():
                if "CPU usage" in line:
                    # Example: CPU usage: 7.82% user, 10.21% sys, 81.95% idle
                    parts = line.split(":", 1)[1].split(",")
                    user = _parse_float(parts[0].strip().split("%")[0]) if len(parts) > 0 else 0.0
                    sys_cpu = _parse_float(parts[1].strip().split("%")[0]) if len(parts) > 1 else 0.0
                    return max(0.0, min(100.0, user + sys_cpu))

        # Linux fallback
        result = subprocess.run(["top", "-bn1"], capture_output=True, text=True, timeout=5, check=False)
        for line in result.stdout.splitlines():
            if "Cpu(s)" in line:
                # Example: %Cpu(s):  7.9 us,  2.1 sy, ...
                user = 0.0
                sys_cpu = 0.0
                segments = line.split(":", 1)[1].split(",") if ":" in line else []
                for seg in segments:
                    seg = seg.strip()
                    if seg.endswith(" us"):
                        user = _parse_float(seg.split()[0])
                    elif seg.endswith(" sy"):
                        sys_cpu = _parse_float(seg.split()[0])
                return max(0.0, min(100.0, user + sys_cpu))
    except Exception:
        return 0.0

    return 0.0


def get_memory_usage_fallback() -> tuple[int, int, float]:
    """Return (total_mb, used_mb, used_percent) without psutil."""
    system = platform.system()

    try:
        if system == "Darwin":
            vm = subprocess.run(["vm_stat"], capture_output=True, text=True, timeout=5, check=False)
            mem = subprocess.run(["sysctl", "-n", "hw.memsize"], capture_output=True, text=True, timeout=5, check=False)

            vm_out = vm.stdout
            total_bytes = int(mem.stdout.strip() or "0")
            if not vm_out or total_bytes <= 0:
                return 0, 0, 0.0

            page_size = 4096
            for line in vm_out.splitlines():
                if "page size of" in line:
                    # Mach Virtual Memory Statistics: (page size of 4096 bytes)
                    digits = "".join(ch for ch in line if ch.isdigit())
                    if digits:
                        page_size = int(digits)
                    break

            def page_count(prefix: str) -> int:
                for ln in vm_out.splitlines():
                    if ln.startswith(prefix):
                        value = ln.split(":", 1)[1].strip().rstrip(".").replace(".", "")
                        return int(value)
                return 0

            active = page_count("Pages active")
            wired = page_count("Pages wired down")
            compressed = page_count("Pages occupied by compressor")
            used_pages = active + wired + compressed

            used_bytes = used_pages * page_size
            total_mb = total_bytes // (1024 * 1024)
            used_mb = used_bytes // (1024 * 1024)
            used_pct = (used_mb / total_mb * 100.0) if total_mb > 0 else 0.0
            return total_mb, used_mb, used_pct

        # Linux fallback
        result = subprocess.run(["free", "-m"], capture_output=True, text=True, timeout=5, check=False)
        for line in result.stdout.splitlines():
            if line.startswith("Mem:"):
                parts = line.split()
                total_mem = int(parts[1])
                used_mem = int(parts[2])
                mem_percent = (used_mem / total_mem) * 100 if total_mem > 0 else 0.0
                return total_mem, used_mem, mem_percent
    except Exception:
        return 0, 0, 0.0

    return 0, 0, 0.0


def get_load_avg_fallback() -> str:
    try:
        if hasattr(os, "getloadavg"):
            return f"{os.getloadavg()[0]:.2f}"

        result = subprocess.run(["uptime"], capture_output=True, text=True, timeout=5, check=False)
        out = result.stdout.strip()
        if "load average:" in out:
            return out.split("load average:", 1)[1].split(",")[0].strip()
        if "load averages:" in out:
            return out.split("load averages:", 1)[1].split(",")[0].strip()
    except Exception:
        return "0.00"

    return "0.00"


def get_system_info() -> dict[str, object]:
    """Get current system resource usage."""
    if psutil is not None:
        try:
            cpu_usage = float(psutil.cpu_percent(interval=0.2))
            mem = psutil.virtual_memory()
            disk = psutil.disk_usage(str(Path.cwd()))
            load_avg = get_load_avg_fallback()

            return {
                "cpu_usage": cpu_usage,
                "mem_total": int(mem.total / (1024 * 1024)),
                "mem_used": int(mem.used / (1024 * 1024)),
                "mem_percent": float(mem.percent),
                "disk_total": f"{disk.total / (1024**3):.1f}G",
                "disk_used": f"{disk.used / (1024**3):.1f}G",
                "disk_avail": f"{disk.free / (1024**3):.1f}G",
                "disk_percent": float(disk.percent),
                "load_avg": load_avg,
                "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "collector": "psutil",
            }
        except Exception:
            # fall through to subprocess implementation
            pass

    # Fallback mode (no psutil or psutil failed)
    cpu_usage = get_cpu_usage_fallback()
    total_mem, used_mem, mem_percent = get_memory_usage_fallback()

    try:
        du = shutil.disk_usage(Path.cwd())
        disk_total = f"{du.total / (1024**3):.1f}G"
        disk_used = f"{du.used / (1024**3):.1f}G"
        disk_avail = f"{du.free / (1024**3):.1f}G"
        disk_percent = (du.used / du.total * 100.0) if du.total > 0 else 0.0
    except Exception:
        disk_total, disk_used, disk_avail, disk_percent = "0G", "0G", "0G", 0.0

    return {
        "cpu_usage": cpu_usage,
        "mem_total": total_mem,
        "mem_used": used_mem,
        "mem_percent": mem_percent,
        "disk_total": disk_total,
        "disk_used": disk_used,
        "disk_avail": disk_avail,
        "disk_percent": disk_percent,
        "load_avg": get_load_avg_fallback(),
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "collector": "fallback",
    }


def draw_progress_bar(percentage: float, width: int = 40, filled_char: str = "█", empty_char: str = "░") -> str:
    """Draw a progress bar."""
    percentage = max(0.0, min(100.0, float(percentage)))
    filled = int(width * percentage / 100)
    empty = width - filled

    if percentage < 50:
        color = "\033[92m"  # Green
    elif percentage < 75:
        color = "\033[93m"  # Yellow
    elif percentage < 90:
        color = "\033[91m"  # Red
    else:
        color = "\033[95m"  # Magenta (critical)

    reset = "\033[0m"
    bar = f"{color}{filled_char * filled}{reset}{empty_char * empty}"
    return f"[{bar}] {percentage:5.1f}%"


def display_dashboard(system_info: dict[str, object], pipeline_log: str | None = None) -> None:
    """Display the resource dashboard."""
    width, _ = get_terminal_size()

    print("=" * width)
    print("🖥️  WGS PIPELINE RESOURCE DASHBOARD".center(width))
    print(f"Last Updated: {system_info['timestamp']}".center(width))
    print(f"Collector: {system_info.get('collector', 'unknown')}".center(width))
    print("=" * width)
    print()

    print("📊 SYSTEM RESOURCES")
    print("-" * 50)

    cpu_bar = draw_progress_bar(float(system_info["cpu_usage"]))
    print(f"CPU Usage:    {cpu_bar}")

    mem_bar = draw_progress_bar(float(system_info["mem_percent"]))
    print(f"Memory:       {mem_bar} ({system_info['mem_used']} / {system_info['mem_total']} MB)")

    disk_bar = draw_progress_bar(float(system_info["disk_percent"]))
    print(f"Disk Usage:   {disk_bar} ({system_info['disk_used']} / {system_info['disk_total']})")

    print(f"Load Average: {system_info['load_avg']}")
    print()

    if pipeline_log and Path(pipeline_log).exists():
        print("🔄 PIPELINE STATUS")
        print("-" * 50)
        try:
            with open(pipeline_log, "r", encoding="utf-8") as f:
                recent_lines = f.readlines()[-10:]

            for line in recent_lines:
                if not line.strip():
                    continue
                if "[STEP_START]" in line:
                    print(f"\033[94m{line.strip()}\033[0m")
                elif "[STEP_COMPLETE]" in line:
                    print(f"\033[92m{line.strip()}\033[0m")
                elif "[STEP_FAILED]" in line or "[ERROR]" in line:
                    print(f"\033[91m{line.strip()}\033[0m")
                elif "[INFO]" in line:
                    print(f"\033[90m{line.strip()}\033[0m")
                else:
                    print(line.strip())
        except Exception as exc:
            print(f"Error reading pipeline log: {exc}")

    print()

    warnings = []
    if float(system_info["cpu_usage"]) > 90:
        warnings.append("⚠️  HIGH CPU USAGE (>90%)")
    if float(system_info["mem_percent"]) > 90:
        warnings.append("⚠️  HIGH MEMORY USAGE (>90%)")
    if float(system_info["disk_percent"]) > 90:
        warnings.append("⚠️  LOW DISK SPACE (<10% free)")

    if warnings:
        print("🚨 SYSTEM WARNINGS")
        print("-" * 50)
        for warning in warnings:
            print(f"\033[91m{warning}\033[0m")
        print()

    print("=" * width)
    print("Press Ctrl+C to exit | Refreshing every 5 seconds".center(width))


def monitor_resources(log_file: str | None = None, refresh_rate: int = 5) -> None:
    """Main monitoring loop."""
    try:
        while True:
            clear_screen()
            system_info = get_system_info()
            display_dashboard(system_info, log_file)
            time.sleep(refresh_rate)
    except KeyboardInterrupt:
        print("\n\n👋 Resource monitoring stopped.")
        sys.exit(0)
    except Exception as exc:
        print(f"\n\n❌ Error in monitoring: {exc}")
        sys.exit(1)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Real-time resource dashboard for WGS pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic monitoring
    python3 scripts/resource_dashboard.py

    # Monitor with pipeline log
    python3 scripts/resource_dashboard.py --log logs/WGS_Pipeline_MySample_progress.log

    # Custom refresh rate
    python3 scripts/resource_dashboard.py --refresh 3
        """,
    )

    parser.add_argument("--log", "-l", help="Pipeline log file to monitor", default=None)
    parser.add_argument("--refresh", "-r", type=int, help="Refresh rate in seconds (default: 5)", default=5)

    args = parser.parse_args()

    print("🚀 Starting WGS Pipeline Resource Dashboard...")
    print("   Use Ctrl+C to exit")
    if psutil is None:
        print("   psutil not found - using fallback collectors")
    time.sleep(1)

    monitor_resources(args.log, max(1, args.refresh))


if __name__ == "__main__":
    main()
