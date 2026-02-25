# WGS Analysis Pipeline

A memory-optimized whole-genome sequencing analysis pipeline for variant discovery, designed to run on systems with as little as 16 GB RAM.

<!-- badges -->
<!-- ![CI](https://github.com/Hilo-Hilo/wgs-analysis-pipeline/actions/workflows/test.yml/badge.svg) -->

## Features

- End-to-end WGS workflow: QC, trimming, alignment, variant calling, annotation
- Optimized for 16 GB RAM with automatic resource monitoring
- Hardware profiles for laptops, workstations, servers, and cloud instances
- Checkpoint/resume support for long-running analyses
- Docker and Docker Compose support for reproducible environments
- Optional SQLite sample registry for small family projects (`n <= 10`)

## Prerequisites

- **OS:** Linux or macOS
- **Conda/Mamba** (Miniconda or Mambaforge)
- **RAM:** 16 GB minimum (14 GB peak during alignment)
- **Disk:** 400 GB free (SSD recommended)
- **CPU:** 4+ cores

Core pipeline tools are installed via conda:
FastQC, fastp, BWA, SAMtools, BCFtools.

`vep` is optional (annotation step only) and may be installed separately on
platforms where base environments prioritize toolchain stability.

See [DEPENDENCIES.md](DEPENDENCIES.md) for the full version policy and supported ranges.

## Quick Start

```bash
git clone https://github.com/Hilo-Hilo/wgs-analysis-pipeline.git
cd wgs-analysis-pipeline
bash scripts/check_requirements.sh          # verify system requirements
bash scripts/validate_deps.sh --ci-mode     # validate dependency versions
# place paired FASTQ files in data/raw/
bash run_pipeline.sh --sample-id MySample --input-dir data/raw --output-dir results
```

## DGX GPU Alignment (Optional)

Alignment defaults to CPU BWA. On DGX systems you can enable GPU alignment
for the alignment step using Parabricks:

```bash
bash run_pipeline.sh \
  --sample-id E200032534 \
  --input-dir /workspace/genetic-data/01-raw-reads/wgs_analysis \
  --output-dir /workspace/genetic-data \
  --use-gpu --gpu-aligner parabricks --gpu-count 1 --steps alignment
```

Notes:
- `--use-gpu` only affects the `alignment` step.
- Requires `nvidia-smi` and `pbrun` in PATH.
- If GPU flags are not provided, pipeline uses CPU BWA path.

## Optional Sample Registry (Small Family Projects)

Track sample metadata in a local SQLite DB:

```bash
python3 wgs registry init --db .wgs/family_registry.db
python3 wgs registry add --db .wgs/family_registry.db --sample-id child_01 --status pending
python3 wgs registry list --db .wgs/family_registry.db --json
```

You can also wire this into pipeline runs without changing default behavior:

```bash
bash run_pipeline.sh --sample-id child_01 --input-dir data/raw --output-dir results \
  --registry-db .wgs/family_registry.db
```

See [SAMPLE_REGISTRY.md](SAMPLE_REGISTRY.md) for schema and full workflow.

## Pipeline Stages

| Stage | Script | Description | Est. Time | Peak RAM |
|-------|--------|-------------|-----------|----------|
| Requirements check | `scripts/check_requirements.sh` | Validate tools and system resources | < 1 min | — |
| Quality control | `scripts/quality_control.sh` | FastQC reports on raw reads | 15–30 min | 2 GB |
| Data cleaning | `scripts/data_cleaning.sh` | Adapter/quality trimming with fastp | 30–60 min | 2 GB |
| Alignment | `scripts/alignment.sh` | BWA-MEM alignment, sort, markdup | 6–10 h | 14 GB |
| Variant calling | `scripts/variant_calling.sh` | BCFtools mpileup + call, filtering | 2–4 h | 3 GB |
| Annotation | `scripts/vep_annotation.sh` | Ensembl VEP functional annotation | 1–2 h | 6 GB |

Estimates are for 30× human WGS (~90 GB FASTQ) on a 4-core system.

## Configuration

The pipeline reads settings from `config/default.conf`. Key parameters:

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `THREADS` | 4 | CPU threads |
| `MAX_MEMORY_BWA` | 10 | GB reserved for BWA |
| `PARALLEL_JOBS` | 1 | Concurrent pipeline steps |
| `KEEP_INTERMEDIATE_FILES` | false | Retain intermediate BAMs/VCFs |
| `ENABLE_CHECKPOINTING` | true | Allow resume after interruption |

Hardware-specific presets are provided in `config/profiles/` for reference:

| Profile | File | Target |
|---------|------|--------|
| Laptop | `config/profiles/laptop.conf` | 16 GB, 4 cores |
| Workstation | `config/profiles/workstation.conf` | 32–64 GB, 8+ cores |
| Server | `config/profiles/server.conf` | 128+ GB, 16+ cores |
| Cloud | `config/profiles/cloud.conf` | Cloud VMs (AWS/GCP) |

The pipeline entrypoint currently uses `config/default.conf` + CLI overrides.

## Docker

Build and run the pipeline in a container:

```bash
docker compose up --build
```

Or build manually:

```bash
docker build -t wgs-pipeline .
docker run -v $(pwd)/data:/app/data -v $(pwd)/results:/app/results wgs-pipeline
```

See `Dockerfile` and `docker-compose.yml` for details.

## Testing

Generate synthetic test data and run the test suite:

```bash
python3 tests/generate_sample_data.py
bash tests/run_tests.sh
```

Registry-focused checks are covered by:

```bash
bash tests/run_tests.sh --unit-only   # includes DB init/CRUD unit tests
bash tests/run_tests.sh --quick       # includes registry smoke integration
```

CI runs automatically via GitHub Actions (`.github/workflows/test.yml`).

## Project Structure

```
.
├── run_pipeline.sh              # Main entry point
├── scripts/
│   ├── check_requirements.sh    # System requirements validation
│   ├── validate_deps.sh         # Dependency version validation
│   ├── quality_control.sh       # FastQC
│   ├── data_cleaning.sh         # fastp trimming
│   ├── alignment.sh             # CPU BWA or GPU Parabricks
│   ├── variant_calling.sh       # BCFtools
│   ├── vep_annotation.sh        # Ensembl VEP
│   ├── load_config.sh           # Config loader
│   ├── validate_input.sh        # Input file validation
│   ├── manage_profiles.sh       # Profile management
│   ├── progress_monitor.sh      # Runtime progress tracking
│   ├── resource_dashboard.py    # System resource monitor
│   └── sample_registry.py       # SQLite sample metadata registry
├── config/
│   ├── default.conf             # Default configuration
│   └── profiles/                # Hardware-specific profiles
├── tests/
│   ├── run_tests.sh             # Test suite
│   ├── test_sample_registry.py  # Registry unit tests (init/CRUD)
│   └── generate_sample_data.py  # Synthetic data generator
├── SAMPLE_REGISTRY.md           # Registry schema + family workflow
├── DEPENDENCIES.md              # Version policy documentation
├── GETTING_STARTED.md
├── TROUBLESHOOTING.md
├── Dockerfile
├── docker-compose.yml
└── LICENSE
```

## License

MIT — see [LICENSE](LICENSE).

## Development

We use `pre-commit` to maintain code quality.

1. **Setup**: `make dev-setup`
2. **Lint**: `make lint` (or auto-runs on commit)

## Citation

```bibtex
@software{wgs_analysis_pipeline,
  author  = {Hilo-Hilo},
  title   = {WGS Analysis Pipeline},
  url     = {https://github.com/Hilo-Hilo/wgs-analysis-pipeline},
  year    = {2025}
}
```
