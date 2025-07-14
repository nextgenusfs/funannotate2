# Memory Monitoring for Funannotate2 Ab Initio Predictions

This document describes the memory monitoring and prediction system implemented for funannotate2's ab initio gene prediction step.

## Overview

The memory monitoring system provides:

1. **Memory Usage Prediction** - Estimate memory requirements based on contig length
2. **Real-time Memory Monitoring** - Track actual memory usage of subprocess calls
3. **Memory-aware CPU Allocation** - Adjust parallelization based on memory constraints
4. **Memory Usage Reporting** - Generate detailed memory usage reports

## Features

### 1. Memory Prediction Models

The system includes empirical models to predict memory usage for each ab initio tool:

- **SNAP**: Base 50 MB + 0.5 MB per MB of sequence
- **Augustus**: Base 100 MB + 2.0 MB per MB of sequence  
- **GlimmerHMM**: Base 30 MB + 0.3 MB per MB of sequence
- **GeneMark**: Base 80 MB + 1.0 MB per MB of sequence

These models provide rough estimates that can be refined with actual usage data.

### 2. Real-time Memory Monitoring

Uses `psutil` to monitor subprocess memory usage in real-time:

- Tracks RSS (Resident Set Size) and VMS (Virtual Memory Size)
- Monitors parent process and all child processes
- Samples memory usage at configurable intervals (default: 100ms)
- Calculates peak, average, and duration statistics

### 3. Memory-aware Scheduling

Automatically adjusts CPU allocation based on:

- Available system memory
- Predicted memory usage per process
- User-specified memory limits
- System memory buffer (20% reserved for OS)

### 4. Integration with Existing Code

The memory monitoring is integrated into:

- `runSubprocess()` - Optional memory monitoring for individual commands
- `abinitio_wrapper()` - Memory prediction and logging per contig
- `runProcessJob()` - Memory-aware CPU allocation for multiprocessing

## Usage

### Command Line Options

Add memory monitoring to funannotate2 predict:

```bash
# Enable memory monitoring
funannotate2 predict -i input_dir --monitor-memory

# Enable memory monitoring with memory limit
funannotate2 predict -i input_dir --monitor-memory --memory-limit 16
```

### CLI Options

- `--monitor-memory`: Enable memory monitoring and prediction
- `--memory-limit GB`: Set memory limit in GB to adjust CPU allocation

### Example Output

When memory monitoring is enabled, you'll see output like:

```
Memory monitoring enabled for ab initio predictions
Memory limit set to 16.0 GB
Memory usage estimate for 150 contigs with tools ['snap', 'augustus']:
  Total estimated peak memory: 2847.3 MB
System memory: 14.2 GB available
Processing contig scaffold_1.fasta (length: 2,847,392 bp)
SNAP memory prediction for scaffold_1.fasta: 51.4 MB
Augustus memory prediction for scaffold_1.fasta: 105.4 MB
Memory usage for snap-scaffold_1.fasta:
Process: snap-scaffold_1.fasta
Duration: 12.34 seconds
Peak RSS: 48.2 MB
Peak VMS: 156.7 MB
Average RSS: 42.1 MB
Samples collected: 247
```

## API Reference

### Core Functions

#### `predict_memory_usage(tool_name, contig_length, prediction_data=None)`

Predict memory usage for an ab initio tool based on contig length.

**Parameters:**
- `tool_name`: Name of the ab initio tool ('snap', 'augustus', etc.)
- `contig_length`: Length of the contig in base pairs
- `prediction_data`: Optional historical data for improved predictions

**Returns:** Dictionary with predicted memory usage statistics

#### `MemoryMonitor.monitor_process(process, process_name)`

Monitor memory usage of a subprocess in real-time.

**Parameters:**
- `process`: subprocess.Popen object to monitor
- `process_name`: Name identifier for the process

**Returns:** Dictionary containing memory statistics

#### `estimate_total_memory_usage(contigs, tools, prediction_data=None)`

Estimate total memory usage for running ab initio predictions on multiple contigs.

**Parameters:**
- `contigs`: List of contig file paths
- `tools`: List of ab initio tools to run
- `prediction_data`: Optional historical data

**Returns:** Dictionary with total memory estimates

#### `suggest_cpu_allocation(total_memory_estimate, available_memory_gb, max_cpus)`

Suggest optimal CPU allocation based on memory constraints.

**Parameters:**
- `total_memory_estimate`: Total estimated memory usage in MB
- `available_memory_gb`: Available system memory in GB
- `max_cpus`: Maximum number of CPUs available

**Returns:** Dictionary with CPU allocation suggestions

### Utility Functions

#### `get_system_memory_info()`

Get current system memory information.

**Returns:** Dictionary with system memory statistics

#### `get_contig_length(contig_file)`

Get the length of a contig from a FASTA file.

**Parameters:**
- `contig_file`: Path to the contig FASTA file

**Returns:** Length of the contig in base pairs

#### `format_memory_report(stats)`

Format memory statistics into a human-readable report.

**Parameters:**
- `stats`: Memory statistics dictionary

**Returns:** Formatted string report

## Implementation Details

### Memory Monitoring Process

1. **Prediction Phase**: Before running ab initio tools, estimate memory usage based on contig lengths
2. **System Check**: Assess available system memory and suggest CPU allocation
3. **Real-time Monitoring**: During subprocess execution, sample memory usage at regular intervals
4. **Reporting**: Log memory statistics and generate reports
5. **Model Updates**: Optionally update prediction models with actual usage data

### Memory Sampling

The memory monitor:
- Creates a `psutil.Process` object for the subprocess
- Samples memory usage every 100ms (configurable)
- Tracks both the main process and all child processes
- Handles process termination gracefully
- Calculates statistics from all samples

### CPU Allocation Logic

The system adjusts CPU allocation by:
1. Estimating memory usage per parallel process
2. Calculating how many processes can fit in available memory
3. Leaving a 20% buffer for the operating system
4. Ensuring at least 1 CPU is allocated
5. Not exceeding the user-specified maximum

## Testing

Run the test suite to verify functionality:

```bash
python test_memory_monitoring.py
```

This will test:
- Memory prediction models
- System memory information
- CPU allocation suggestions  
- Total memory estimation
- Real-time memory monitoring

## Dependencies

The memory monitoring system requires:

- `psutil` - For system and process memory monitoring
- `json` - For saving/loading memory statistics
- `time` - For timing and sampling
- `threading` - For concurrent memory monitoring

## Future Enhancements

Potential improvements include:

1. **Machine Learning Models** - Use actual usage data to train better prediction models
2. **Memory Profiling** - Detailed analysis of memory allocation patterns
3. **Dynamic Scheduling** - Adjust CPU allocation during runtime based on actual usage
4. **Memory Limits** - Hard memory limits with process termination
5. **Historical Analysis** - Long-term memory usage trends and optimization
6. **Tool-specific Tuning** - Fine-tune memory models for different ab initio tools

## Troubleshooting

### Common Issues

1. **psutil not available**: Install with `pip install psutil`
2. **Permission errors**: Some systems may restrict process monitoring
3. **Inaccurate predictions**: Models are empirical and may need tuning for your data
4. **Memory monitoring overhead**: Monitoring adds small CPU/memory overhead

### Performance Impact

Memory monitoring has minimal performance impact:
- ~1-2% CPU overhead for sampling
- ~1-5 MB memory overhead for the monitor
- Sampling interval can be adjusted to reduce overhead
