import json
import os
import psutil
import subprocess
import threading
import time
from collections import defaultdict
from typing import Dict, List, Optional, Tuple, Union

from .utilities import human_readable_size


class MemoryMonitor:
    """
    Monitor memory usage of subprocesses and track statistics.

    This class provides functionality to monitor memory usage of subprocesses
    in real-time, collect statistics, and predict memory requirements based
    on input characteristics like contig length.
    """

    def __init__(self, sampling_interval: float = 0.1):
        """
        Initialize the memory monitor.

        Args:
            sampling_interval: How often to sample memory usage (seconds)
        """
        self.sampling_interval = sampling_interval
        self.memory_stats = defaultdict(list)
        self.prediction_data = defaultdict(list)
        self._monitoring = False
        self._monitor_thread = None

    def monitor_process(
        self, process: subprocess.Popen, process_name: str = "unknown"
    ) -> Dict[str, Union[float, int]]:
        """
        Monitor memory usage of a subprocess in real-time.

        Args:
            process: The subprocess.Popen object to monitor
            process_name: Name identifier for the process

        Returns:
            Dictionary containing memory statistics
        """
        try:
            psutil_process = psutil.Process(process.pid)
            memory_samples = []
            start_time = time.time()

            # Monitor until process completes
            while process.poll() is None:
                try:
                    # Get memory info for the process and all its children
                    memory_info = psutil_process.memory_info()
                    rss_mb = memory_info.rss / (1024 * 1024)  # Convert to MB
                    vms_mb = memory_info.vms / (1024 * 1024)  # Convert to MB

                    # Also check children processes
                    children_rss = 0
                    children_vms = 0
                    try:
                        for child in psutil_process.children(recursive=True):
                            child_memory = child.memory_info()
                            children_rss += child_memory.rss / (1024 * 1024)
                            children_vms += child_memory.vms / (1024 * 1024)
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass

                    total_rss = rss_mb + children_rss
                    total_vms = vms_mb + children_vms

                    memory_samples.append(
                        {
                            "timestamp": time.time() - start_time,
                            "rss_mb": total_rss,
                            "vms_mb": total_vms,
                        }
                    )

                    time.sleep(self.sampling_interval)

                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    # Process might have ended
                    break

            # Calculate statistics
            if memory_samples:
                rss_values = [sample["rss_mb"] for sample in memory_samples]
                vms_values = [sample["vms_mb"] for sample in memory_samples]

                stats = {
                    "process_name": process_name,
                    "duration_seconds": time.time() - start_time,
                    "peak_rss_mb": max(rss_values),
                    "peak_vms_mb": max(vms_values),
                    "avg_rss_mb": sum(rss_values) / len(rss_values),
                    "avg_vms_mb": sum(vms_values) / len(vms_values),
                    "sample_count": len(memory_samples),
                    "samples": memory_samples,
                }

                # Store for future analysis
                self.memory_stats[process_name].append(stats)

                return stats
            else:
                return {
                    "process_name": process_name,
                    "duration_seconds": 0,
                    "peak_rss_mb": 0,
                    "peak_vms_mb": 0,
                    "avg_rss_mb": 0,
                    "avg_vms_mb": 0,
                    "sample_count": 0,
                    "samples": [],
                }

        except Exception as e:
            return {
                "process_name": process_name,
                "error": str(e),
                "peak_rss_mb": 0,
                "peak_vms_mb": 0,
                "avg_rss_mb": 0,
                "avg_vms_mb": 0,
                "sample_count": 0,
                "samples": [],
            }


def get_contig_length(contig_file: str) -> int:
    """
    Get the length of a contig from a FASTA file.

    Args:
        contig_file: Path to the contig FASTA file

    Returns:
        Length of the contig in base pairs
    """
    try:
        with open(contig_file, "r") as f:
            lines = f.readlines()
            sequence_lines = [
                line.strip() for line in lines if not line.startswith(">")
            ]
            return sum(len(line) for line in sequence_lines)
    except Exception:
        return 0


def predict_memory_usage(
    tool_name: str, contig_length: int, prediction_data: Dict = None
) -> Dict[str, float]:
    """
    Predict memory usage for an ab initio tool based on contig length.

    Based on real memory monitoring data from funannotate2 runs (2024-07-14).
    Values for SNAP and GlimmerHMM are from actual measurements on 8 contigs (2.8-4.9 Mbp).
    Values for Augustus and GeneMark are conservative estimates due to Docker measurement
    limitations on Apple Silicon.

    Args:
        tool_name: Name of the ab initio tool (snap, augustus, etc.)
        contig_length: Length of the contig in base pairs
        prediction_data: Historical data for predictions (optional)

    Returns:
        Dictionary with predicted memory usage statistics
    """
    # Base memory usage (MB) - minimum memory required regardless of contig size
    # Updated with real data from memory monitoring
    base_memory = {
        "snap": 234.0,  # From real measurements on 8 contigs (2.8-4.9 Mbp)
        "augustus": 250.0,  # Conservative estimate (Docker measurement issues)
        "glimmerhmm": 356.0,  # From real measurements on 8 contigs (2.8-4.9 Mbp)
        "genemark": 200.0,  # Conservative estimate (Docker measurement issues)
    }

    # Memory scaling per MB of sequence (converted from MB per million base pairs)
    # Updated with real data from memory monitoring
    memory_per_mb = {
        "snap": 149.0,  # From real measurements: 149.26 MB/Mbp
        "augustus": 150.0,  # Conservative estimate (similar to SNAP)
        "glimmerhmm": 164.0,  # From real measurements: 163.61 MB/Mbp
        "genemark": 140.0,  # Conservative estimate (slightly lower than SNAP)
    }

    tool = tool_name.lower()
    if tool not in base_memory:
        tool = "augustus"  # Default to augustus if unknown

    # Convert contig length to MB (assuming 1 byte per base)
    contig_mb = contig_length / (1024 * 1024)

    # Simple linear prediction
    predicted_peak = base_memory[tool] + (memory_per_mb[tool] * contig_mb)

    # Add some safety margin (20%)
    predicted_peak_with_margin = predicted_peak * 1.2

    return {
        "tool": tool_name,
        "contig_length_bp": contig_length,
        "contig_size_mb": contig_mb,
        "predicted_peak_mb": predicted_peak,
        "predicted_peak_with_margin_mb": predicted_peak_with_margin,
        "base_memory_mb": base_memory[tool],
        "memory_per_mb_sequence": memory_per_mb[tool],
    }


def save_memory_stats(stats: Dict, output_file: str):
    """
    Save memory statistics to a JSON file.

    Args:
        stats: Memory statistics dictionary
        output_file: Path to output JSON file
    """
    try:
        with open(output_file, "w") as f:
            json.dump(stats, f, indent=2, default=str)
    except Exception as e:
        print(f"Warning: Could not save memory stats to {output_file}: {e}")


def load_memory_stats(input_file: str) -> Dict:
    """
    Load memory statistics from a JSON file.

    Args:
        input_file: Path to input JSON file

    Returns:
        Memory statistics dictionary
    """
    try:
        with open(input_file, "r") as f:
            return json.load(f)
    except Exception:
        return {}


def format_memory_report(stats: Dict) -> str:
    """
    Format memory statistics into a human-readable report.

    Args:
        stats: Memory statistics dictionary

    Returns:
        Formatted string report
    """
    if "error" in stats:
        return f"Memory monitoring failed: {stats['error']}"

    report = []
    report.append(f"Process: {stats['process_name']}")
    report.append(f"Duration: {stats['duration_seconds']:.2f} seconds")
    report.append(
        f"Peak RSS: {human_readable_size(stats['peak_rss_mb'] * 1024 * 1024)}"
    )
    report.append(
        f"Peak VMS: {human_readable_size(stats['peak_vms_mb'] * 1024 * 1024)}"
    )
    report.append(
        f"Average RSS: {human_readable_size(stats['avg_rss_mb'] * 1024 * 1024)}"
    )
    report.append(f"Samples collected: {stats['sample_count']}")

    return "\n".join(report)


def estimate_total_memory_usage(
    contigs: List[str], tools: List[str], prediction_data: Dict = None
) -> Dict[str, float]:
    """
    Estimate total memory usage for running ab initio predictions on multiple contigs.

    Args:
        contigs: List of contig file paths
        tools: List of ab initio tools to run
        prediction_data: Historical data for predictions (optional)

    Returns:
        Dictionary with total memory estimates
    """
    total_estimates = {
        "total_predicted_peak_mb": 0,
        "total_predicted_with_margin_mb": 0,
        "per_contig_estimates": [],
        "per_tool_totals": defaultdict(float),
    }

    for contig in contigs:
        contig_length = get_contig_length(contig)
        contig_name = os.path.basename(contig)

        contig_estimates = {
            "contig": contig_name,
            "length_bp": contig_length,
            "tools": {},
        }

        for tool in tools:
            prediction = predict_memory_usage(tool, contig_length, prediction_data)
            contig_estimates["tools"][tool] = prediction

            # Add to totals (assuming tools run sequentially per contig)
            total_estimates["per_tool_totals"][tool] += prediction[
                "predicted_peak_with_margin_mb"
            ]

        # For this contig, the peak would be the maximum across tools
        contig_peak = max(
            [
                contig_estimates["tools"][tool]["predicted_peak_with_margin_mb"]
                for tool in tools
            ]
        )
        total_estimates["total_predicted_peak_mb"] += contig_peak

        contig_estimates["contig_peak_mb"] = contig_peak
        total_estimates["per_contig_estimates"].append(contig_estimates)

    total_estimates["total_predicted_with_margin_mb"] = total_estimates[
        "total_predicted_peak_mb"
    ]

    return total_estimates


def suggest_cpu_allocation(
    total_memory_estimate: float, available_memory_gb: float, max_cpus: int
) -> Dict[str, Union[int, float]]:
    """
    Suggest optimal CPU allocation based on memory constraints.

    Args:
        total_memory_estimate: Total estimated memory usage in MB
        available_memory_gb: Available system memory in GB
        max_cpus: Maximum number of CPUs available

    Returns:
        Dictionary with CPU allocation suggestions
    """
    available_memory_mb = available_memory_gb * 1024

    # Leave some memory for the system (20% buffer)
    usable_memory_mb = available_memory_mb * 0.8

    # Calculate how many parallel processes we can run
    if total_memory_estimate > 0:
        max_parallel_by_memory = int(usable_memory_mb / total_memory_estimate)
    else:
        max_parallel_by_memory = max_cpus

    # Don't exceed the CPU limit
    suggested_cpus = min(max_parallel_by_memory, max_cpus)

    # Ensure at least 1 CPU
    suggested_cpus = max(1, suggested_cpus)

    return {
        "suggested_cpus": suggested_cpus,
        "max_cpus": max_cpus,
        "memory_limited": max_parallel_by_memory < max_cpus,
        "available_memory_mb": available_memory_mb,
        "usable_memory_mb": usable_memory_mb,
        "estimated_memory_per_process_mb": total_memory_estimate,
        "memory_utilization_percent": (
            suggested_cpus * total_memory_estimate / usable_memory_mb
        )
        * 100,
    }


def update_prediction_model(
    tool_name: str, contig_length: int, actual_memory_mb: float, stats_file: str = None
):
    """
    Update the prediction model with actual memory usage data.

    Args:
        tool_name: Name of the ab initio tool
        contig_length: Length of the contig in base pairs
        actual_memory_mb: Actual peak memory usage in MB
        stats_file: Optional file to save/load historical data
    """
    # Load existing data if available
    historical_data = {}
    if stats_file and os.path.exists(stats_file):
        historical_data = load_memory_stats(stats_file)

    # Initialize tool data if not present
    if "prediction_data" not in historical_data:
        historical_data["prediction_data"] = {}
    if tool_name not in historical_data["prediction_data"]:
        historical_data["prediction_data"][tool_name] = []

    # Add new data point
    data_point = {
        "contig_length_bp": contig_length,
        "actual_memory_mb": actual_memory_mb,
        "timestamp": time.time(),
    }

    historical_data["prediction_data"][tool_name].append(data_point)

    # Save updated data
    if stats_file:
        save_memory_stats(historical_data, stats_file)

    return historical_data


def get_system_memory_info() -> Dict[str, float]:
    """
    Get current system memory information.

    Returns:
        Dictionary with system memory statistics
    """
    try:
        memory = psutil.virtual_memory()
        return {
            "total_gb": memory.total / (1024**3),
            "available_gb": memory.available / (1024**3),
            "used_gb": memory.used / (1024**3),
            "free_gb": memory.free / (1024**3),
            "percent_used": memory.percent,
            "buffers_gb": getattr(memory, "buffers", 0) / (1024**3),
            "cached_gb": getattr(memory, "cached", 0) / (1024**3),
        }
    except Exception as e:
        return {
            "error": str(e),
            "total_gb": 0,
            "available_gb": 0,
            "used_gb": 0,
            "free_gb": 0,
            "percent_used": 0,
            "buffers_gb": 0,
            "cached_gb": 0,
        }
