# Sample Registry (SQLite)

A lightweight, optional metadata registry for small-family WGS projects (`n <= 10`).

## Why this exists

For family projects, a full LIMS is overkill. This registry tracks the essentials:

- sample ID
- FASTQ paths
- output path
- reference genome
- run status
- notes
- timestamps

The pipeline works normally without the registry.

---

## Schema

Table: `sample_registry`

| Column | Type | Notes |
|---|---|---|
| `sample_id` | `TEXT PRIMARY KEY` | canonical sample identifier |
| `fastq_r1_path` | `TEXT` | path to R1 FASTQ |
| `fastq_r2_path` | `TEXT` | path to R2 FASTQ |
| `output_dir` | `TEXT` | result/output directory |
| `reference_genome` | `TEXT` | reference FASTA path |
| `run_status` | `TEXT` | one of: `pending`, `running`, `completed`, `failed`, `on-hold` |
| `notes` | `TEXT` | freeform comments |
| `created_at` | `TEXT` | UTC ISO-8601 timestamp |
| `updated_at` | `TEXT` | UTC ISO-8601 timestamp |
| `last_run_at` | `TEXT` | UTC ISO-8601 timestamp for run-state updates |

Index: `idx_sample_registry_status(run_status)`

---

## Commands

You can call the registry directly:

```bash
python3 scripts/sample_registry.py init --db .wgs/family_registry.db
```

Or through the unified CLI wrapper:

```bash
python3 wgs registry init --db .wgs/family_registry.db
```

### Init DB

```bash
python3 wgs registry init --db .wgs/family_registry.db
```

### Add sample

```bash
python3 wgs registry add \
  --db .wgs/family_registry.db \
  --sample-id child_01 \
  --fastq-r1 data/raw/child_01_R1.fastq.gz \
  --fastq-r2 data/raw/child_01_R2.fastq.gz \
  --reference refs/GRCh38.fa \
  --output-dir results/child_01 \
  --status pending \
  --notes "Trio child"
```

### Update sample

```bash
python3 wgs registry update \
  --db .wgs/family_registry.db \
  --sample-id child_01 \
  --status running \
  --notes "Alignment started" \
  --append-notes
```

### List samples

```bash
python3 wgs registry list --db .wgs/family_registry.db
python3 wgs registry list --db .wgs/family_registry.db --json
python3 wgs registry list --db .wgs/family_registry.db --status failed --json
```

---

## Pipeline hook (optional, non-breaking)

`run_pipeline.sh` now accepts:

```bash
--registry-db PATH
```

When provided, the pipeline will:

1. initialize the DB if needed
2. ensure a sample record exists
3. set status to `running` at start
4. set status to `completed` or `failed` at end
5. append status notes

If registry operations fail, pipeline execution continues (warning-only).

---

## Small-family workflow

```bash
# 1) Initialize once
python3 wgs registry init --db .wgs/family_registry.db

# 2) Register all family members
python3 wgs registry add --db .wgs/family_registry.db --sample-id mom --status pending
python3 wgs registry add --db .wgs/family_registry.db --sample-id dad --status pending
python3 wgs registry add --db .wgs/family_registry.db --sample-id child --status pending

# 3) Run per sample with optional hook
bash run_pipeline.sh \
  --sample-id child \
  --input-dir data/raw \
  --output-dir results/child \
  --registry-db .wgs/family_registry.db

# 4) Review status board
python3 wgs registry list --db .wgs/family_registry.db
```
