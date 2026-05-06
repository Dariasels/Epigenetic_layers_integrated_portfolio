#!/usr/bin/env python3
"""
Import GSE102538 H3K27ac peak-count files, store raw counts, and build
condition-aware enhancer and plasticity summaries.

Inputs:
  - data/GSE102538_H3K27ac_EntorhinalCortex.bed.gz
  - data/chip_GSE102538 /*.txt.gz

Outputs in MySQL (`brain_multiomics`):
  - chip_peak_counts: raw peak counts per sample
  - enhancers: adds H3K27ac AD/Control summary columns
  - plasticity_summary: adds gene-level H3K27ac summary columns

Usage:
  DB_PASSWORD=simba python3 STEP3_additionaldata/import_chipseq_condition_counts.py
"""

from __future__ import annotations

import gzip
import os
import re
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

import mysql.connector


DB_CONFIG = {
    "host": "localhost",
    "user": "daria",
    "password": os.getenv("DB_PASSWORD", ""),
    "database": "brain_multiomics",
}

COUNT_DIR = os.path.join("data", "chip_GSE102538 ")
PEAK_BED = os.path.join("data", "GSE102538_H3K27ac_EntorhinalCortex.bed.gz")
DATASET = "GSE102538"


@dataclass
class PeakSummary:
    ad_sum: float = 0.0
    ad_n: int = 0
    control_sum: float = 0.0
    control_n: int = 0

    def add(self, condition: str, value: float) -> None:
        if condition == "Alzheimer":
            self.ad_sum += value
            self.ad_n += 1
        elif condition == "Control":
            self.control_sum += value
            self.control_n += 1

    @property
    def ad_mean(self):
        return self.ad_sum / self.ad_n if self.ad_n else None

    @property
    def control_mean(self):
        return self.control_sum / self.control_n if self.control_n else None

    @property
    def delta(self):
        if self.ad_mean is None or self.control_mean is None:
            return None
        return self.ad_mean - self.control_mean


def connect():
    return mysql.connector.connect(**DB_CONFIG)


def normalize_condition(raw: str) -> str:
    raw = (raw or "").strip().lower()
    if raw.startswith("ad") or "alzheimer" in raw:
        return "Alzheimer"
    if raw.startswith("control"):
        return "Control"
    return "unknown"


def sample_condition_from_filename(filename: str) -> str:
    name = filename.upper()
    if re.search(r"_AD\d+", name):
        return "Alzheimer"
    if re.search(r"_CONTROL\d+", name):
        return "Control"
    return "unknown"


def open_text(path: str):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.endswith(".gz") else open(path, "r", encoding="utf-8", errors="replace")


def load_peak_map(peak_bed: str) -> Dict[str, str]:
    """Map Peak_1-style identifiers to the enhancer peak_name used in MySQL."""
    peak_map: Dict[str, str] = {}
    with open_text(peak_bed) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            cols = line.split()
            if len(cols) < 4:
                continue
            start = int(cols[1])
            end = int(cols[2])
            peak_id = cols[3]
            enhancer_peak_name = f"{DATASET}_H3K27ac_EntorhinalCortex_{start}_{end}"
            peak_map[peak_id] = enhancer_peak_name
    return peak_map


def ensure_schema(conn):
    cursor = conn.cursor()
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS chip_peak_counts (
            count_id     INT AUTO_INCREMENT PRIMARY KEY,
            sample_id    VARCHAR(50) NOT NULL,
            peak_name    VARCHAR(100) NOT NULL,
            count_value  DOUBLE NOT NULL,
            dataset      VARCHAR(50) NOT NULL,
            condition_label VARCHAR(20) NOT NULL,
            created_at   TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            UNIQUE KEY uq_sample_peak (sample_id, peak_name),
            KEY idx_peak_name (peak_name),
            KEY idx_dataset (dataset),
            KEY idx_condition (condition_label),
            CONSTRAINT fk_chip_counts_sample
                FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
                ON DELETE CASCADE
        )
        """
    )

    for sql in [
        "ALTER TABLE enhancers ADD COLUMN h3k27ac_ad_mean DOUBLE NULL",
        "ALTER TABLE enhancers ADD COLUMN h3k27ac_control_mean DOUBLE NULL",
        "ALTER TABLE enhancers ADD COLUMN h3k27ac_delta DOUBLE NULL",
        "ALTER TABLE enhancers ADD COLUMN h3k27ac_ad_samples INT NULL",
        "ALTER TABLE enhancers ADD COLUMN h3k27ac_control_samples INT NULL",
        "ALTER TABLE plasticity_summary ADD COLUMN h3k27ac_ad_mean DOUBLE NULL",
        "ALTER TABLE plasticity_summary ADD COLUMN h3k27ac_control_mean DOUBLE NULL",
        "ALTER TABLE plasticity_summary ADD COLUMN h3k27ac_delta DOUBLE NULL",
        "ALTER TABLE plasticity_summary ADD COLUMN h3k27ac_enhancer_count INT NULL",
    ]:
        try:
            cursor.execute(sql)
        except mysql.connector.Error:
            pass

    conn.commit()
    cursor.close()


def normalize_existing_metadata(conn):
    cursor = conn.cursor()
    for table in ("samples", "metadata_import"):
        cursor.execute(f"UPDATE {table} SET `condition`='Alzheimer' WHERE dataset=%s AND UPPER(`condition`)='AD'", (DATASET,))
        cursor.execute(f"UPDATE {table} SET `condition`='Control' WHERE dataset=%s AND LOWER(`condition`)='control'", (DATASET,))
    conn.commit()
    cursor.close()


def list_count_files() -> List[str]:
    return sorted(
        os.path.join(COUNT_DIR, fn)
        for fn in os.listdir(COUNT_DIR)
        if fn.endswith(".txt.gz")
    )


def import_counts(conn, files: Iterable[str], peak_map: Dict[str, str], write_raw: bool = True):
    cursor = conn.cursor()
    sql = """
        INSERT IGNORE INTO chip_peak_counts
            (sample_id, peak_name, count_value, dataset, condition_label)
        VALUES (%s, %s, %s, %s, %s)
    """
    batch = []
    batch_size = 10000
    total = 0
    summaries: Dict[str, PeakSummary] = defaultdict(PeakSummary)

    for filepath in files:
        filename = os.path.basename(filepath)
        m = re.match(r"^(GSM\d+)", filename)
        if not m:
            continue
        sample_id = m.group(1)
        condition = sample_condition_from_filename(filename)
        if condition == "unknown":
            continue

        with open_text(filepath) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track"):
                    continue
                cols = line.split()
                if len(cols) < 2:
                    continue
                peak_name = cols[0]
                try:
                    count_value = float(cols[1])
                except ValueError:
                    continue
                enhancer_peak_name = peak_map.get(peak_name)
                if not enhancer_peak_name:
                    continue
                if write_raw:
                    batch.append((sample_id, peak_name, count_value, DATASET, condition))
                    if len(batch) >= batch_size:
                        cursor.executemany(sql, batch)
                        conn.commit()
                        total += len(batch)
                        batch.clear()
                summaries[enhancer_peak_name].add(condition, count_value)

    if write_raw and batch:
        cursor.executemany(sql, batch)
        conn.commit()
        total += len(batch)

    cursor.close()
    return total, summaries


def build_enhancer_summary(conn, summaries: Dict[str, PeakSummary]):
    cursor = conn.cursor(dictionary=True)
    cursor.execute(
        """
        SELECT enhancer_id, peak_name
        FROM enhancers
        WHERE dataset = %s
        """,
        (DATASET,),
    )
    enhancer_rows = cursor.fetchall()
    cursor.close()

    enhancer_updates = []
    enhancer_to_summary = {}
    for row in enhancer_rows:
        peak_name = row["peak_name"]
        summary = summaries.get(peak_name)
        if not summary:
            continue
        enhancer_id = row["enhancer_id"]
        enhancer_to_summary[enhancer_id] = summary
        enhancer_updates.append(
            (
                summary.ad_mean,
                summary.control_mean,
                summary.delta,
                summary.ad_n,
                summary.control_n,
                enhancer_id,
            )
        )

    cursor = conn.cursor()
    cursor.executemany(
        """
        UPDATE enhancers
        SET h3k27ac_ad_mean=%s,
            h3k27ac_control_mean=%s,
            h3k27ac_delta=%s,
            h3k27ac_ad_samples=%s,
            h3k27ac_control_samples=%s
        WHERE enhancer_id=%s
        """,
        enhancer_updates,
    )
    conn.commit()
    cursor.close()

    return enhancer_to_summary


def update_plasticity_summary(conn, enhancer_to_summary: Dict[int, PeakSummary]):
    cursor = conn.cursor(dictionary=True)
    cursor.execute(
        """
        SELECT egl.gene_symbol,
               egl.enhancer_id
        FROM enhancer_gene_links egl
        JOIN enhancers e ON e.enhancer_id = egl.enhancer_id
        WHERE e.dataset = %s
        """,
        (DATASET,),
    )
    gene_map: Dict[str, List[int]] = defaultdict(list)
    for row in cursor.fetchall():
        gene_map[row["gene_symbol"]].append(row["enhancer_id"])
    cursor.close()

    rows = []
    for gene_symbol, enhancer_ids in gene_map.items():
        ad_vals = []
        control_vals = []
        for enhancer_id in enhancer_ids:
            summary = enhancer_to_summary.get(enhancer_id)
            if not summary:
                continue
            if summary.ad_mean is not None:
                ad_vals.append(summary.ad_mean)
            if summary.control_mean is not None:
                control_vals.append(summary.control_mean)
        if not ad_vals and not control_vals:
            continue
        ad_mean = sum(ad_vals) / len(ad_vals) if ad_vals else None
        control_mean = sum(control_vals) / len(control_vals) if control_vals else None
        delta = None if ad_mean is None or control_mean is None else ad_mean - control_mean
        rows.append((ad_mean, control_mean, delta, len(set(enhancer_ids)), gene_symbol))

    cursor = conn.cursor()
    cursor.executemany(
        """
        UPDATE plasticity_summary
        SET h3k27ac_ad_mean=%s,
            h3k27ac_control_mean=%s,
            h3k27ac_delta=%s,
            h3k27ac_enhancer_count=%s
        WHERE gene_symbol=%s
        """,
        rows,
    )
    conn.commit()
    cursor.close()


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-raw-import", action="store_true", help="Do not write to chip_peak_counts; rebuild summaries only")
    args = parser.parse_args()

    conn = connect()
    try:
        normalize_existing_metadata(conn)
        ensure_schema(conn)
        peak_map = load_peak_map(PEAK_BED)
        files = list_count_files()
        inserted, summaries = import_counts(conn, files, peak_map, write_raw=not args.skip_raw_import)
        enhancer_to_summary = build_enhancer_summary(conn, summaries)
        update_plasticity_summary(conn, enhancer_to_summary)
        if args.skip_raw_import:
            print("Skipped raw peak-count import")
        else:
            print(f"Imported {inserted:,} raw peak-count rows")
        print(f"Summarised {len(enhancer_to_summary):,} enhancers")
    finally:
        conn.close()


if __name__ == "__main__":
    main()