"""
02b_import_rnaseq_GSE33000_FINAL.py
=====================================
Imports GSE33000 into rna_expression using GSM IDs as sample_ids.

KEY FACTS FROM THE SERIES MATRIX FILE:
  - Platform: GPL4372 (custom Agilent 44K array) — NOT GPL10558 Illumina
  - Data table columns: already GSM IDs (GSM1423780, GSM1423781, ...)
  - Probe IDs: numeric (e.g. 10019475365) in the ID_REF column
  - Sample titles: PFC_1 ... PFC_624 (internal names, NOT what we store)
  - Source names ch1: HBTRC_PF_Pool_1 ... (pool names, ignore)
  - GSM range: GSM1423780 – GSM1424403 (624 samples)

PROBE → GENE MAPPING:
  GPL4372 is a custom array. The series matrix file does NOT have a
  Gene column in the data table. We need to get gene symbols from
  the GPL4372 platform annotation on GEO.

  Download command (run once):
    wget "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL4nnn/GPL4372/soft/GPL4372.annot.gz"
    gunzip GPL4372.annot.gz

  The annotation file has columns including:
    ID          → probe ID (matches ID_REF in data table)
    Gene symbol → HGNC symbol  ← what we want
    Gene title  → full gene name

  If the wget fails (NCBI paths change), use GEOquery in R:
    library(GEOquery); gpl <- getGEO("GPL4372")
    write.csv(Table(gpl)[,c("ID","Gene symbol")], "GPL4372_annotation.csv", row.names=FALSE)

Usage:
  python 02b_import_rnaseq_GSE33000_FINAL.py \
    --file GSE33000_series_matrix.txt \
    --annotation GPL4372.annot

  # Test with 1000 probes first:
  python 02b_import_rnaseq_GSE33000_FINAL.py \
    --file GSE33000_series_matrix.txt \
    --annotation GPL4372.annot \
    --limit 1000
"""

import os
import argparse
import math
import pandas as pd
import mysql.connector
import gzip
import io
import sys
import csv

DB_CONFIG = {
    "host":     "localhost",
    "user":     "daria",
    "password": os.getenv("DB_PASSWORD", ""),
    "database": "brain_multiomics"
}

def connect():
    conn = mysql.connector.connect(**DB_CONFIG)
    print("✓ Connected to database")
    return conn

def open_file(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") \
           if path.endswith(".gz") else \
           open(path, "r", encoding="utf-8", errors="replace")

# ── STEP 1: parse the series matrix data table ────────────────────────────────
def parse_series_matrix(filepath):
    """
    Extract only the data between !series_matrix_table_begin and
    !series_matrix_table_end. Column headers are GSM IDs already.
    Returns a DataFrame: index = probe ID, columns = GSM IDs.
    """
    print(f"  Parsing series matrix: {filepath}")
    lines = []
    in_table = False

    with open_file(filepath) as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                lines.append(line)

    if not lines:
        print("✗ Could not find data table in file.")
        print("  Check that the file is a proper series matrix.")
        sys.exit(1)

    df = pd.read_csv(io.StringIO("".join(lines)), sep="\t",
                     index_col=0, low_memory=False)

    # Strip quotes from column names (GEO sometimes includes them)
    df.columns = df.columns.str.strip('"')
    df.index   = df.index.astype(str).str.strip('"')

    print(f"  ✓ Data table: {df.shape[0]:,} probes × {df.shape[1]} samples")
    print(f"  ✓ Sample ID example: {df.columns[0]}")
    print(f"  ✓ Probe ID example:  {df.index[0]}")
    return df

# ── STEP 2: load probe → gene annotation ─────────────────────────────────────
def load_annotation(filepath):
    """
    Load GPL4372 annotation file. Handles both .annot (GEO format with
    header lines starting with #/!) and plain CSV from GEOquery R export.
    Returns dict: probe_id (str) → gene_symbol (str)
    """
    print(f"\n  Loading annotation: {filepath}")

    # Detect format
    with open_file(filepath) as fh:
        first_lines = [next(fh) for _ in range(5)]

    is_geo_annot = any(l.startswith("!") or l.startswith("#")
                       for l in first_lines)

    if is_geo_annot:
        # GEO .annot format: skip lines starting with ! or #,
        # find !platform_table_begin, read from there
        lines = []
        in_table = False
        with open_file(filepath) as fh:
            for line in fh:
                if line.startswith("!platform_table_begin"):
                    in_table = True
                    continue
                if line.startswith("!platform_table_end"):
                    break
                if in_table:
                    lines.append(line)
        df = pd.read_csv(io.StringIO("".join(lines)), sep="\t",
                         low_memory=False)
    else:
        # Plain CSV (from R GEOquery export)
        df = pd.read_csv(filepath, low_memory=False)

    df.columns = df.columns.str.strip()
    print(f"  Annotation columns: {list(df.columns[:8])}")

    # Find ID column
    col_lower = {c.lower().strip(): c for c in df.columns}
    id_col  = col_lower.get("id", col_lower.get("probe_id"))
    # Find gene symbol column — GPL4372 uses "Gene symbol"
    sym_col = col_lower.get("gene symbol",
              col_lower.get("gene_symbol",
              col_lower.get("symbol")))

    if not id_col:
        print(f"✗ Cannot find probe ID column. Available: {list(df.columns[:10])}")
        sys.exit(1)
    if not sym_col:
        print(f"✗ Cannot find gene symbol column. Available: {list(df.columns[:10])}")
        print("  Try the R download method described in the script header.")
        sys.exit(1)

    print(f"  Using: id='{id_col}', symbol='{sym_col}'")

    probe_map = {}
    no_symbol = 0
    for _, row in df.iterrows():
        pid  = str(row[id_col]).strip()
        gene = str(row.get(sym_col, "")).strip()
        if gene and gene not in ("nan", "---", "", "N/A"):
            # Take first symbol if multiple (e.g. "BDNF /// NTRK2")
            gene = gene.replace("///", ";").split(";")[0].strip()
            probe_map[pid] = gene
        else:
            no_symbol += 1

    print(f"  ✓ {len(probe_map):,} probes mapped to gene symbols")
    print(f"  ⚠ {no_symbol:,} probes have no symbol (controls/unnamed)")
    return probe_map

# ── STEP 3: ensure all GSM samples exist in samples table ────────────────────
def ensure_samples(gsm_ids, conn):
    """
    Insert GSM IDs into samples table. condition='unknown' for now —
    the metadata update script will fill in AD/control from metadata_import.
    """
    cursor = conn.cursor()

    # Check which samples are already in metadata_import so we can
    # set condition correctly if available
    cursor.execute("SHOW COLUMNS FROM metadata_import")
    columns = {row[0].lower() for row in cursor.fetchall()}
    label_col = None
    if "condition" in columns:
        label_col = "condition"
    elif "disease" in columns:
        label_col = "disease"

    meta = {}
    if label_col:
        cursor.execute(
            f"SELECT sample_id, `{label_col}` FROM metadata_import WHERE dataset='GSE33000'"
        )
        for sid, label in cursor.fetchall():
            if label:
                c = str(label).lower()
                if "alzheimer" in c or " ad" in c or "dementia" in c:
                    meta[sid] = "AD"
                elif "control" in c or "normal" in c or "non-demented" in c:
                    meta[sid] = "control"
                else:
                    meta[sid] = label

    sql = ("INSERT IGNORE INTO samples "
           "(sample_id, `condition`, tissue, dataset) VALUES (%s,%s,%s,%s)")
    rows = []
    for gsm in gsm_ids:
        cond = meta.get(gsm, "unknown")
        rows.append((gsm, cond, "prefrontal cortex", "GSE33000"))

    cursor.executemany(sql, rows)
    conn.commit()

    known = sum(1 for r in rows if r[1] != "unknown")
    print(f"  ✓ {len(rows)} samples inserted/updated")
    print(f"  ✓ {known} already had condition labels from metadata_import")
    if known == 0:
        if label_col:
            print("  ⚠ No condition labels found — run metadata import first, then:")
            print("    UPDATE samples s JOIN metadata_import m ON m.sample_id=s.sample_id")
            print(f"    SET s.condition = CASE WHEN LOWER(m.{label_col}) LIKE '%alzheimer%'")
            print(f"    THEN 'AD' WHEN LOWER(m.{label_col}) LIKE '%control%' THEN 'control'")
            print(f"    ELSE m.{label_col} END WHERE s.dataset='GSE33000';")
        else:
            print("  ⚠ metadata_import has no 'condition' or 'disease' column;")
            print("    cannot prefill sample labels from metadata_import.")
    cursor.close()

# ── STEP 4: import expression values ─────────────────────────────────────────
def import_expression(df, probe_map, conn, limit=None):
    """
    For each probe with a known gene symbol, insert one row per sample
    into rna_expression. Uses GSM IDs directly as sample_id.
    """
    cursor = conn.cursor()
    sql = ("INSERT IGNORE INTO rna_expression "
           "(gene_symbol, sample_id, expression_value) VALUES (%s,%s,%s)")

    gsm_ids       = list(df.columns)
    batch         = []
    BATCH_SIZE    = 8000
    total         = 0
    skipped       = 0
    probes_done   = 0

    print(f"\n  Importing {len(df):,} probes × {len(gsm_ids)} samples...")
    print(f"  (Only probes with gene symbol annotations will be inserted)\n")

    for probe_id, row in df.iterrows():
        gene = probe_map.get(str(probe_id).strip())
        if not gene:
            skipped += 1
            continue

        for gsm in gsm_ids:
            try:
                val = float(row[gsm])
            except (ValueError, TypeError):
                continue
            if not math.isfinite(val):
                continue
            batch.append((gene, gsm, val))

        probes_done += 1
        if len(batch) >= BATCH_SIZE:
            cursor.executemany(sql, batch)
            conn.commit()
            total += len(batch)
            batch  = []
            print(f"  Probes: {probes_done:>6,} | Rows inserted: {total:>10,}", end="\r")

        if limit and probes_done >= limit:
            print(f"\n  --limit {limit} reached.")
            break

    if batch:
        cursor.executemany(sql, batch)
        conn.commit()
        total += len(batch)

    cursor.close()
    print(f"\n  ✓ Probes with gene symbol : {probes_done:,}")
    print(f"  ✓ Rows inserted           : {total:,}")
    print(f"  ⚠ Probes skipped (no gene): {skipped:,}")

# ── MAIN ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file",       required=True,
                        help="Path to GSE33000_series_matrix.txt(.gz)")
    parser.add_argument("--annotation", required=True,
                        help="Path to GPL4372 annotation file")
    parser.add_argument("--limit",      type=int, default=None,
                        help="Stop after N probes (for testing)")
    args = parser.parse_args()

    print("\n=== Importing GSE33000 (GPL4372 Agilent array) ===\n")
    if args.limit:
        print(f"  TEST MODE: first {args.limit} probes only\n")

    df        = parse_series_matrix(args.file)
    probe_map = load_annotation(args.annotation)
    conn      = connect()
    ensure_samples(list(df.columns), conn)
    import_expression(df, probe_map, conn, limit=args.limit)
    conn.close()

    print("\n=== Done ===")
    print("\nVerify in MySQL:")
    print("  SELECT COUNT(DISTINCT sample_id) FROM rna_expression;")
    print("  SELECT COUNT(DISTINCT gene_symbol) FROM rna_expression;")
    print("  SELECT gene_symbol, sample_id, expression_value")
    print("    FROM rna_expression WHERE gene_symbol='BDNF' LIMIT 5;")

if __name__ == "__main__":
    main()
