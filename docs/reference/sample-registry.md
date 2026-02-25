# Sample Registry (SQLite)

A lightweight optional metadata registry for small-family WGS projects (`n <= 10`).

## Schema

Table: `sample_registry`

| Column | Type | Notes |
|---|---|---|
| `sample_id` | `TEXT PRIMARY KEY` | canonical sample identifier |
| `fastq_r1_path` | `TEXT` | path to R1 FASTQ |
| `fastq_r2_path` | `TEXT` | path to R2 FASTQ |
| `output_dir` | `TEXT` | result/output directory |
| `reference_genome` | `TEXT` | reference FASTA path |
| `run_status` | `TEXT` | `pending`, `running`, `completed`, `failed`, `on-hold` |
| `notes` | `TEXT` | freeform comments |
| `created_at` | `TEXT` | UTC ISO-8601 timestamp |
| `updated_at` | `TEXT` | UTC ISO-8601 timestamp |
| `last_run_at` | `TEXT` | UTC ISO-8601 status-update timestamp |

Index: `idx_sample_registry_status(run_status)`

## CLI Usage

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

## Pipeline Hook

`run_pipeline.sh` supports:

```bash
--registry-db PATH
```

When enabled, pipeline state transitions are written automatically (`running` → `completed`/`failed`) without blocking the main run if registry writes fail.
