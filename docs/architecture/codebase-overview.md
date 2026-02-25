# Codebase Overview

## Purpose

A memory-aware whole-genome sequencing pipeline for QC → cleaning → alignment → variant calling → annotation, with optional GPU alignment and optional sample registry tracking.

## Top-Level Layout

```text
.
├── run_pipeline.sh          # Main orchestrator
├── wgs                      # Unified CLI wrapper
├── scripts/                 # Stage scripts + utilities
├── config/                  # Default config + hardware profiles
├── tests/                   # Unit/integration/smoke tests
├── docs/                    # Architecture/guides/reference docs
├── .github/                 # CI/release workflows + conda env
├── Dockerfile               # Reproducible execution image
└── docker-compose.yml       # Compose convenience entry
```

## Core Runtime Flow

1. `run_pipeline.sh` parses/validates CLI flags and hardware mode.
2. Stage scripts execute in order:
   - `quality_control.sh`
   - `data_cleaning.sh`
   - `alignment.sh` (CPU: BWA-MEM2, GPU: Parabricks)
   - `variant_calling.sh`
   - `vep_annotation.sh` (optional)
3. Checkpoint/resume metadata and logs are updated.
4. Optional sample registry updates status transitions.

## Key Modules

- `scripts/check_requirements.sh`: preflight environment checks.
- `scripts/validate_deps.sh`: enforces dependency policy from `DEPENDENCIES.md`.
- `scripts/detect_hardware.sh`: CPU/RAM/GPU detection and profile hints.
- `scripts/sample_registry.py`: SQLite metadata tracking utility.

## Testing Surface

- `tests/run_tests.sh --unit-only`: fast correctness checks.
- `tests/run_tests.sh`: includes integration flow.
- `tests/smoke/run_smoke.sh`: E2E smoke, mock or real toolchain.

## CI Notes

- Main workflow: `.github/workflows/test.yml`.
- Matrix unit tests: Ubuntu 22.04 + macOS 14.
- Integration environment: `.github/conda/integration-environment.yml`.
