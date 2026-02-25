#!/usr/bin/env python3
"""Unit tests for SQLite sample registry CRUD flows."""

from __future__ import annotations

import json
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class SampleRegistryCliTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmpdir = tempfile.TemporaryDirectory(prefix="wgs_registry_unit_")
        self.db_path = Path(self.tmpdir.name) / "family" / "registry.db"
        self.script_path = (
            Path(__file__).resolve().parents[1] / "scripts" / "sample_registry.py"
        )

    def tearDown(self) -> None:
        self.tmpdir.cleanup()

    def run_cli(self, *args: str, expect_ok: bool = True) -> subprocess.CompletedProcess:
        cmd = [sys.executable, str(self.script_path), *args]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if expect_ok and proc.returncode != 0:
            self.fail(
                "CLI command failed:\n"
                f"cmd={' '.join(cmd)}\n"
                f"stdout={proc.stdout}\n"
                f"stderr={proc.stderr}"
            )
        if not expect_ok and proc.returncode == 0:
            self.fail(
                "CLI command unexpectedly succeeded:\n"
                f"cmd={' '.join(cmd)}\n"
                f"stdout={proc.stdout}\n"
                f"stderr={proc.stderr}"
            )
        return proc

    def test_init_creates_schema(self) -> None:
        self.run_cli("init", "--db", str(self.db_path))
        self.assertTrue(self.db_path.exists(), "DB file should exist after init")

        with sqlite3.connect(self.db_path) as conn:
            columns = conn.execute("PRAGMA table_info(sample_registry)").fetchall()

        column_names = {col[1] for col in columns}
        expected = {
            "sample_id",
            "fastq_r1_path",
            "fastq_r2_path",
            "output_dir",
            "reference_genome",
            "run_status",
            "notes",
            "created_at",
            "updated_at",
            "last_run_at",
        }
        self.assertTrue(expected.issubset(column_names))

    def test_crud_flow_add_list_update(self) -> None:
        self.run_cli("init", "--db", str(self.db_path))

        self.run_cli(
            "add",
            "--db",
            str(self.db_path),
            "--sample-id",
            "child_01",
            "--fastq-r1",
            "data/raw/child_01_R1.fastq.gz",
            "--fastq-r2",
            "data/raw/child_01_R2.fastq.gz",
            "--reference",
            "refs/GRCh38.fa",
            "--status",
            "pending",
            "--notes",
            "Imported trio sample",
        )

        listed = self.run_cli("list", "--db", str(self.db_path), "--json")
        rows = json.loads(listed.stdout)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["sample_id"], "child_01")
        self.assertEqual(rows[0]["run_status"], "pending")

        self.run_cli(
            "update",
            "--db",
            str(self.db_path),
            "--sample-id",
            "child_01",
            "--status",
            "running",
            "--notes",
            "Aligned with BWA",
            "--append-notes",
        )

        listed = self.run_cli(
            "list",
            "--db",
            str(self.db_path),
            "--sample-id",
            "child_01",
            "--json",
        )
        rows = json.loads(listed.stdout)
        self.assertEqual(rows[0]["run_status"], "running")
        self.assertIn("Imported trio sample", rows[0]["notes"])
        self.assertIn("Aligned with BWA", rows[0]["notes"])
        self.assertIsNotNone(rows[0]["updated_at"])

    def test_duplicate_add_requires_allow_existing(self) -> None:
        self.run_cli("init", "--db", str(self.db_path))

        self.run_cli(
            "add",
            "--db",
            str(self.db_path),
            "--sample-id",
            "parent_01",
        )

        proc = self.run_cli(
            "add",
            "--db",
            str(self.db_path),
            "--sample-id",
            "parent_01",
            expect_ok=False,
        )
        self.assertIn("already exists", proc.stderr)

        self.run_cli(
            "add",
            "--db",
            str(self.db_path),
            "--sample-id",
            "parent_01",
            "--status",
            "completed",
            "--allow-existing",
        )

        listed = self.run_cli(
            "list",
            "--db",
            str(self.db_path),
            "--sample-id",
            "parent_01",
            "--json",
        )
        rows = json.loads(listed.stdout)
        self.assertEqual(rows[0]["run_status"], "completed")

    def test_update_nonexistent_sample_fails(self) -> None:
        self.run_cli("init", "--db", str(self.db_path))
        proc = self.run_cli(
            "update",
            "--db",
            str(self.db_path),
            "--sample-id",
            "missing_sample",
            "--status",
            "failed",
            expect_ok=False,
        )
        self.assertIn("does not exist", proc.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
