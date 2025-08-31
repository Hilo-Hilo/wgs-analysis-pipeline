# WGS Pipeline Configuration Profiles

This directory contains optimized configuration profiles for different system types. Each profile is tuned for specific hardware configurations to maximize performance while maintaining stability.

## Available Profiles

### 🔽 Laptop Profile (`laptop.conf`)
**Target Systems:** 8-16GB RAM, 2-4 cores, limited disk space
- **Optimization Focus:** Resource conservation, stability over speed
- **Thread Usage:** Conservative (2 threads)
- **Memory Limits:** 4GB maximum
- **Quality Settings:** Aggressive filtering to reduce data size
- **Disk Management:** Maximum compression, aggressive cleanup
- **Best For:** Development, testing, personal analysis on laptops

### ⚖️ Workstation Profile (`workstation.conf`)
**Target Systems:** 32-64GB RAM, 8-16 cores, good disk space
- **Optimization Focus:** Balanced performance and resource usage
- **Thread Usage:** Moderate (8 threads)
- **Memory Limits:** 16GB maximum
- **Quality Settings:** Standard quality thresholds
- **Disk Management:** Balanced compression, selective cleanup
- **Best For:** Desktop workstations, research environments, production analysis

### 🚀 Server Profile (`server.conf`)
**Target Systems:** 128GB+ RAM, 32+ cores, fast storage
- **Optimization Focus:** Maximum performance and throughput
- **Thread Usage:** High (32 threads)
- **Memory Limits:** 64GB maximum
- **Quality Settings:** Comprehensive analysis
- **Disk Management:** Minimal compression for speed
- **Best For:** High-performance servers, HPC clusters, high-throughput environments

### ☁️ Cloud Profile (`cloud.conf`)
**Target Systems:** Variable cloud resources (AWS, GCP, Azure)
- **Optimization Focus:** Cost efficiency and reliability
- **Thread Usage:** Cost-optimized (16 threads)
- **Memory Limits:** 32GB maximum
- **Quality Settings:** Efficient processing
- **Disk Management:** Network-optimized, cost-aware cleanup
- **Best For:** Cloud computing environments, spot instances, cost-sensitive analysis

## Quick Start

### 1. List Available Profiles
```bash
scripts/manage_profiles.sh list
```

### 2. View Profile Details
```bash
scripts/manage_profiles.sh show laptop
```

### 3. Set Active Profile
```bash
scripts/manage_profiles.sh set workstation
```

### 4. Run Pipeline with Active Profile
```bash
./run_pipeline.sh --sample-id MySample
```

## Profile Management

### Setting a Profile
When you set a profile, it becomes the default configuration for all pipeline runs:
```bash
scripts/manage_profiles.sh set server
```

### Viewing Current Profile
```bash
scripts/manage_profiles.sh list
```
The current active profile will be indicated in the output.

### Creating Custom Profiles
```bash
scripts/manage_profiles.sh create my_custom_profile
```
This creates a new profile based on the workstation template that you can customize.

### Validating Profiles
```bash
scripts/manage_profiles.sh validate laptop
```

### Comparing Profiles
```bash
scripts/manage_profiles.sh compare laptop server
```

## Configuration Integration

Profiles are automatically loaded by the pipeline scripts through the configuration system:

1. **Automatic Loading:** If an active profile is set, it's loaded automatically
2. **Fallback:** If no profile is active, the default configuration is used
3. **Override:** You can still specify custom configurations when needed

## Key Configuration Areas

### Resource Management
- **Memory Limits:** Tool-specific memory caps to prevent system overload
- **Thread Control:** Optimized thread counts for each system type
- **Disk Management:** Storage optimization and cleanup strategies

### Quality Control
- **Quality Thresholds:** Balanced between thoroughness and resource usage
- **Filtering Aggressiveness:** Adjusted for system capabilities
- **Compression Levels:** Optimized for speed vs. storage trade-offs

### Processing Behavior
- **Checkpointing:** Recovery mechanisms for long-running analyses
- **Intermediate Files:** Retention policies based on available storage
- **Monitoring:** Resource tracking appropriate for each system type

### Error Handling
- **Retry Logic:** Appropriate retry counts for system reliability
- **Validation Levels:** Thoroughness of input validation
- **Recovery Mechanisms:** System-appropriate error recovery

## System Requirements by Profile

| Profile | Min RAM | Min Cores | Min Disk | Typical Runtime* |
|---------|---------|-----------|----------|------------------|
| Laptop | 8GB | 2 | 100GB | 8-12 hours |
| Workstation | 32GB | 8 | 200GB | 3-5 hours |
| Server | 128GB | 32 | 500GB | 1-2 hours |
| Cloud | 64GB | 16 | 200GB | 2-4 hours |

*Runtime estimates for 30x coverage human WGS data

## Best Practices

### Profile Selection
1. **Start Conservative:** Use the laptop profile for testing, even on powerful systems
2. **Match Your Hardware:** Choose the profile that best matches your available resources
3. **Consider Workload:** Account for other processes running on your system
4. **Monitor Resources:** Use the built-in resource monitoring to validate your choice

### Custom Profiles
1. **Base on Existing:** Start with the closest matching profile as a template
2. **Test Thoroughly:** Validate custom profiles with small datasets first
3. **Document Changes:** Keep notes on customizations for future reference
4. **Share Profiles:** Consider contributing useful profiles back to the project

### Troubleshooting
1. **Resource Issues:** Switch to a more conservative profile if experiencing memory/disk problems
2. **Slow Performance:** Verify your system matches the profile requirements
3. **Failed Jobs:** Check logs for resource-related errors and adjust profile accordingly
4. **Validation Errors:** Use the validation command to check profile syntax

## Advanced Features

### Cloud-Specific Options
The cloud profile includes specific optimizations for cloud environments:
- Spot instance handling
- Object storage integration
- Network optimization
- Cost monitoring

### HPC Integration
The server profile includes options for HPC environments:
- SLURM integration (optional)
- NUMA optimization
- Hugepage support
- MPI capabilities (future)

## Support and Troubleshooting

For issues with profiles:
1. Check profile validation: `scripts/manage_profiles.sh validate <profile>`
2. Compare with working profile: `scripts/manage_profiles.sh compare <profile1> <profile2>`
3. Reset to default: `scripts/manage_profiles.sh reset`
4. Create custom profile: `scripts/manage_profiles.sh create <name>`

## Profile Development

To contribute new profiles or improvements:
1. Follow the existing profile structure
2. Test thoroughly on target systems
3. Document system requirements and optimizations
4. Validate configuration syntax
5. Submit via pull request with testing results