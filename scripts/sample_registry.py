#!/usr/bin/env python3
"""SQLite-backed sample registry for small-family WGS projects.

This module intentionally keeps scope small and robust for n<=10 samples.
"""

from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ALLOWED_STATUSES = ("pending", "running", "completed", "failed", "on-hold")

SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS sample_registry (
    sample_id TEXT PRIMARY KEY,
    fastq_r1_path TEXT,
    fastq_r2_path TEXT,
    output_dir TEXT,
    reference_genome TEXT,
    run_status TEXT NOT NULL DEFAULT 'pending',
    notes TEXT,
    created_at TEXT NOT NULL,
    updated_at TEXT NOT NULL,
    last_run_at TEXT
);

CREATE INDEX IF NOT EXISTS idx_sample_registry_status
    ON sample_registry(run_status);
"""


class RegistryError(RuntimeError):
    """Raised for expected, user-facing registry errors."""


@dataclass
class SampleRegistry:
    db_path: Path

    def _connect(self) -> sqlite3.Connection:
        conn = sqlite3.connect(str(self.db_path))
        conn.row_factory = sqlite3.Row
        return conn

    @staticmethod
    def _now_iso() -> str:
        return datetime.now(timezone.utc).replace(microsecond=0).isoformat()

    def init_db(self) -> None:
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        with self._connect() as conn:
            conn.executescript(SCHEMA_SQL)

    def _ensure_initialized(self) -> None:
        if not self.db_path.exists():
            raise RegistryError(
                f"Registry DB not found: {self.db_path}. Run 'init' first."
            )

        with self._connect() as conn:
            table = conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name='sample_registry'"
            ).fetchone()
            if table is None:
                raise RegistryError(
                    f"Registry DB exists but schema is missing: {self.db_path}. "
                    "Run 'init' to create schema."
                )

    @staticmethod
    def _validate_status(status: str) -> None:
        if status not in ALLOWED_STATUSES:
            allowed = ", ".join(ALLOWED_STATUSES)
            raise RegistryError(
                f"Invalid run status '{status}'. Allowed values: {allowed}."
            )

    def add_sample(
        self,
        sample_id: str,
        *,
        fastq_r1_path: str | None,
        fastq_r2_path: str | None,
        output_dir: str | None,
        reference_genome: str | None,
        run_status: str,
        notes: str | None,
        allow_existing: bool = False,
    ) -> None:
        self._ensure_initialized()

        normalized_id = sample_id.strip()
        if not normalized_id:
            raise RegistryError("sample_id must be non-empty")

        self._validate_status(run_status)
        now = self._now_iso()

        needs_update = False
        with self._connect() as conn:
            try:
                conn.execute(
                    """
                    INSERT INTO sample_registry (
                        sample_id,
                        fastq_r1_path,
                        fastq_r2_path,
                        output_dir,
                        reference_genome,
                        run_status,
                        notes,
                        created_at,
                        updated_at,
                        last_run_at
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        normalized_id,
                        fastq_r1_path,
                        fastq_r2_path,
                        output_dir,
                        reference_genome,
                        run_status,
                        notes,
                        now,
                        now,
                        now if run_status in {"running", "completed", "failed"} else None,
                    ),
                )
            except sqlite3.IntegrityError as exc:
                if not allow_existing:
                    raise RegistryError(
                        f"Sample '{normalized_id}' already exists. Use update instead."
                    ) from exc
                needs_update = True

        if needs_update:
            # Keep this deterministic: update mutable fields for existing row.
            self.update_sample(
                normalized_id,
                fastq_r1_path=fastq_r1_path,
                fastq_r2_path=fastq_r2_path,
                output_dir=output_dir,
                reference_genome=reference_genome,
                run_status=run_status,
                notes=notes,
                append_notes=False,
            )

    def update_sample(
        self,
        sample_id: str,
        *,
        fastq_r1_path: str | None = None,
        fastq_r2_path: str | None = None,
        output_dir: str | None = None,
        reference_genome: str | None = None,
        run_status: str | None = None,
        notes: str | None = None,
        append_notes: bool = False,
    ) -> None:
        self._ensure_initialized()

        normalized_id = sample_id.strip()
        if not normalized_id:
            raise RegistryError("sample_id must be non-empty")

        fields: dict[str, Any] = {}
        if fastq_r1_path is not None:
            fields["fastq_r1_path"] = fastq_r1_path
        if fastq_r2_path is not None:
            fields["fastq_r2_path"] = fastq_r2_path
        if output_dir is not None:
            fields["output_dir"] = output_dir
        if reference_genome is not None:
            fields["reference_genome"] = reference_genome
        if run_status is not None:
            self._validate_status(run_status)
            fields["run_status"] = run_status

        with self._connect() as conn:
            current = conn.execute(
                "SELECT notes FROM sample_registry WHERE sample_id = ?", (normalized_id,)
            ).fetchone()
            if current is None:
                raise RegistryError(f"Sample '{normalized_id}' does not exist.")

            if notes is not None:
                if append_notes and current["notes"]:
                    fields["notes"] = f"{current['notes']}\n{notes}"
                else:
                    fields["notes"] = notes

            if not fields:
                raise RegistryError("No update fields provided.")

            fields["updated_at"] = self._now_iso()
            if run_status in {"running", "completed", "failed"}:
                fields["last_run_at"] = self._now_iso()

            set_clause = ", ".join(f"{key} = ?" for key in fields)
            values = list(fields.values()) + [normalized_id]

            conn.execute(
                f"UPDATE sample_registry SET {set_clause} WHERE sample_id = ?",
                values,
            )

    def list_samples(
        self,
        *,
        sample_id: str | None = None,
        status: str | None = None,
        limit: int | None = None,
    ) -> list[dict[str, Any]]:
        self._ensure_initialized()

        query = "SELECT * FROM sample_registry"
        clauses = []
        values: list[Any] = []

        if sample_id:
            clauses.append("sample_id = ?")
            values.append(sample_id)

        if status:
            self._validate_status(status)
            clauses.append("run_status = ?")
            values.append(status)

        if clauses:
            query += " WHERE " + " AND ".join(clauses)

        query += " ORDER BY updated_at DESC, sample_id ASC"

        if limit is not None:
            if limit <= 0:
                raise RegistryError("limit must be > 0")
            query += " LIMIT ?"
            values.append(limit)

        with self._connect() as conn:
            rows = conn.execute(query, values).fetchall()
            return [dict(row) for row in rows]


def _add_db_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--db",
        default=".wgs/sample_registry.db",
        help="Path to SQLite registry database (default: .wgs/sample_registry.db)",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Manage a lightweight SQLite sample registry for WGS runs"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # init
    init_parser = subparsers.add_parser("init", help="Initialize registry database")
    _add_db_argument(init_parser)

    # list
    list_parser = subparsers.add_parser("list", help="List sample records")
    _add_db_argument(list_parser)
    list_parser.add_argument("--sample-id", help="Filter by sample ID")
    list_parser.add_argument("--status", choices=ALLOWED_STATUSES, help="Filter by run status")
    list_parser.add_argument("--limit", type=int, help="Maximum records to return")
    list_parser.add_argument("--json", action="store_true", help="Emit JSON")

    # add
    add_parser = subparsers.add_parser("add", help="Add a sample record")
    _add_db_argument(add_parser)
    add_parser.add_argument("--sample-id", required=True, help="Sample identifier")
    add_parser.add_argument("--fastq-r1", dest="fastq_r1_path", help="R1 FASTQ path")
    add_parser.add_argument("--fastq-r2", dest="fastq_r2_path", help="R2 FASTQ path")
    add_parser.add_argument("--output-dir", help="Output directory path")
    add_parser.add_argument("--reference", dest="reference_genome", help="Reference genome path")
    add_parser.add_argument(
        "--status",
        dest="run_status",
        default="pending",
        choices=ALLOWED_STATUSES,
        help="Run status",
    )
    add_parser.add_argument("--notes", help="Freeform notes")
    add_parser.add_argument(
        "--allow-existing",
        action="store_true",
        help="If sample exists, update mutable fields instead of failing",
    )
    add_parser.add_argument("--quiet", action="store_true", help="Suppress success output")

    # update
    update_parser = subparsers.add_parser("update", help="Update a sample record")
    _add_db_argument(update_parser)
    update_parser.add_argument("--sample-id", required=True, help="Sample identifier")
    update_parser.add_argument("--fastq-r1", dest="fastq_r1_path", help="R1 FASTQ path")
    update_parser.add_argument("--fastq-r2", dest="fastq_r2_path", help="R2 FASTQ path")
    update_parser.add_argument("--output-dir", help="Output directory path")
    update_parser.add_argument("--reference", dest="reference_genome", help="Reference genome path")
    update_parser.add_argument("--status", dest="run_status", choices=ALLOWED_STATUSES, help="Run status")
    update_parser.add_argument("--notes", help="Freeform notes")
    update_parser.add_argument(
        "--append-notes",
        action="store_true",
        help="Append new notes to existing notes (newline-separated)",
    )
    update_parser.add_argument("--quiet", action="store_true", help="Suppress success output")

    return parser


def _print_table(rows: list[dict[str, Any]]) -> None:
    if not rows:
        print("No sample records found.")
        return

    headers = [
        "sample_id",
        "run_status",
        "reference_genome",
        "fastq_r1_path",
        "fastq_r2_path",
        "output_dir",
        "updated_at",
    ]
    print("\t".join(headers))
    for row in rows:
        print("\t".join(str(row.get(h) or "") for h in headers))


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    db_path = Path(args.db).expanduser().resolve()
    registry = SampleRegistry(db_path=db_path)

    try:
        if args.command == "init":
            registry.init_db()
            print(f"Initialized sample registry at: {db_path}")
            return 0

        if args.command == "list":
            rows = registry.list_samples(
                sample_id=args.sample_id,
                status=args.status,
                limit=args.limit,
            )
            if args.json:
                print(json.dumps(rows, indent=2, sort_keys=True))
            else:
                _print_table(rows)
            return 0

        if args.command == "add":
            registry.add_sample(
                args.sample_id,
                fastq_r1_path=args.fastq_r1_path,
                fastq_r2_path=args.fastq_r2_path,
                output_dir=args.output_dir,
                reference_genome=args.reference_genome,
                run_status=args.run_status,
                notes=args.notes,
                allow_existing=args.allow_existing,
            )
            if not args.quiet:
                print(f"Added sample '{args.sample_id}'")
            return 0

        if args.command == "update":
            registry.update_sample(
                args.sample_id,
                fastq_r1_path=args.fastq_r1_path,
                fastq_r2_path=args.fastq_r2_path,
                output_dir=args.output_dir,
                reference_genome=args.reference_genome,
                run_status=args.run_status,
                notes=args.notes,
                append_notes=args.append_notes,
            )
            if not args.quiet:
                print(f"Updated sample '{args.sample_id}'")
            return 0

        parser.print_help()
        return 1
    except RegistryError as exc:
        print(f"Registry error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
