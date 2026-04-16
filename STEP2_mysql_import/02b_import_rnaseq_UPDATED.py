"""
02b_import_rnaseq_UPDATED.py
=============================
Updated importer for GSE33000 — fixes TWO problems:
  1. The series_matrix file has ~80 header lines of metadata BEFORE the actual
     data table. Simply renaming to .csv and opening it reads those header lines
     as columns, which is why you only see 2 columns instead of 1200.
  2. The 'Gene' column contains RSE_######### IDs (Illumina probe transcript IDs),
     NOT gene symbols. We need the GPL10558 annotation file to map them.

PLATFORM: GPL10558 — Illumina HumanHT-12 V4.0 Expression BeadChip
PROBE IDs: Illumina numeric IDs (e.g. 10019475365), NOT RSE_ IDs.
           The RSE_######### values in the 'Gene'/'Transcript' columns are
           Illumina transcript cluster IDs — not standard gene symbols.
           The GPL10558 annotation file maps:
               Illumina probe ID → ILMN_Gene → Symbol → RefSeq → Entrez

HOW TO GET THE ANNOTATION FILE (do this once before running this script):
---------------------------------------------------------------------------
Option A — Download via GEO FTP (simplest):
  wget "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz"
  gunzip GPL10558.annot.gz

Option B — Use GEOquery in R (most reliable):
  library(GEOquery)
  gpl <- getGEO("GPL10558")
  anno <- Table(gpl)
  write.csv(anno[, c("ID", "Symbol", "ILMN_Gene", "Entrez_Gene_ID")],
            "GPL10558_annotation.csv", row.names=FALSE)

  The key columns in GPL10558 are:
    ID          → the numeric probe ID (matches 'reporterID' in your data)
    Symbol      → HGNC gene symbol  ← THIS IS WHAT WE WANT
    ILMN_Gene   → Illumina gene name (similar but not always same as Symbol)
    Entrez_Gene_ID → NCBI Entrez ID

Usage:
  python 02b_import_rnaseq_UPDATED.py \
    --file GSE33000_series_matrix.txt \
    --annotation GPL10558_annotation.csv

Requirements:
  pip install mysql-connector-python pandas
"""

import os
import argparse
import pandas as pd
import mysql.connector
import sys
import gzip
import io

DB_CONFIG = {
    "host":     "localhost",
    "user":     "daria",
    "password": os.getenv("DB_PASSWORD", ""),   # ← CHANGE THIS
    "database": "brain_multiomics"
}

# ── STEP 1 ─────────────────────────────────────────────────────────────────────
def parse_series_matrix(filepath):
    """
    Parse a GEO Series Matrix file correctly.

    The file has TWO sections:
      1. Metadata lines starting with '!' (e.g. !Series_title, !Sample_geo_accession)
         These must be SKIPPED — they cause the "only 2 columns" problem.
      2. The actual data table, delimited by:
           !series_matrix_table_begin
           ... tab-separated expression data ...
           !series_matrix_table_end

    We extract ONLY the data table section.
    """
    print(f"  Parsing series matrix file: {filepath}")

    # Handle both plain .txt and .gz compressed files
    if filepath.endswith(".gz"):
        opener = lambda: gzip.open(filepath, "rt", encoding="utf-8", errors="replace")
    else:
        opener = lambda: open(filepath, "r", encoding="utf-8", errors="replace")

    data_lines = []
    in_table   = False

    with opener() as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                data_lines.append(line)

    if not data_lines:
        print("  ✗ ERROR: Could not find '!series_matrix_table_begin' in file.")
        print("    This means either:")
        print("    a) You opened the supplementary expression file (not the series matrix)")
        print("    b) The file was renamed but is not a proper series matrix file")
        print("    Try downloading GSE33000_series_matrix.txt.gz directly from GEO.")
        sys.exit(1)

    print(f"  ✓ Found data table with {len(data_lines)-1} data rows")  # -1 for header

    # Parse the collected lines as a TSV
    content = "".join(data_lines)
    df = pd.read_csv(io.StringIO(content), sep="\t", low_memory=False)

    # The first column is always ID_REF in a series matrix file
    print(f"  ✓ Data shape: {df.shape[0]} probes × {df.shape[1]} columns")
    print(f"  ✓ First column name: '{df.columns[0]}'")
    print(f"  ✓ Sample column example: '{df.columns[1]}'")
    return df


# ── STEP 2 ─────────────────────────────────────────────────────────────────────
# def load_annotation(annotation_path):
#     """
#     Load the GPL10558 annotation file and build a probe_id → gene_symbol dict.

#     The annotation file (from GEO) has these relevant columns:
#       ID       → numeric Illumina probe ID (e.g. 10019475365)
#       Symbol   → HGNC gene symbol (e.g. BDNF, APP, CAMK2A)
#       ILMN_Gene → Illumina gene name (fallback if Symbol is empty)

#     Some probes map to multiple genes (semicolon-separated) — we take the first.
#     Probes with no gene symbol are kept as None (will be skipped at import).
#     """
#     print(f"\n  Loading probe annotation from: {annotation_path}")

#     if annotation_path.endswith(".gz"):
#         df = pd.read_csv(annotation_path, sep=None, engine="python",
#                          comment="#", compression="gzip")
#     else:
#         df = pd.read_csv(annotation_path, sep=None, engine="python",
#                          comment="#")

#     print(f"  Annotation columns: {list(df.columns)}")

#     # Normalise column names (GEO uses different capitalisation sometimes)
#     df.columns = df.columns.str.strip()
#     col_map = {c.lower(): c for c in df.columns}

#     id_col  = col_map.get("id",     col_map.get("probe_id", None))
#     sym_col = col_map.get("symbol", col_map.get("gene_symbol",
#               col_map.get("gene symbol", None)))
#     ilmn_col = col_map.get("ilmn_gene", col_map.get("gene_name", None))

#     if not id_col:
#         print("  ✗ Could not find probe ID column in annotation file.")
#         print("    Expected a column named 'ID' or 'Probe_Id'")
#         sys.exit(1)

#     print(f"  Using: ID col='{id_col}', Symbol col='{sym_col}', "
#           f"ILMN_Gene col='{ilmn_col}'")

#     probe_to_gene = {}
#     no_symbol     = 0

#     for _, row in df.iterrows():
#         probe_id = str(row[id_col]).strip()

#         # Try Symbol first, fall back to ILMN_Gene
#         gene = None
#         if sym_col and str(row.get(sym_col, "")).strip() not in ("", "nan"):
#             gene = str(row[sym_col]).strip()
#             # Take first gene if multiple (e.g. "BDNF///NTRK2")
#             gene = gene.replace("///", ";").split(";")[0].strip()
#         elif ilmn_col and str(row.get(ilmn_col, "")).strip() not in ("", "nan"):
#             gene = str(row[ilmn_col]).strip()

#         if gene and gene not in ("nan", "---", ""):
#             probe_to_gene[probe_id] = gene
#         else:
#             no_symbol += 1

#     print(f"  ✓ Mapped {len(probe_to_gene):,} probes to gene symbols")
#     print(f"  ⚠ {no_symbol:,} probes have no gene symbol (will be skipped)")
#     return probe_to_gene
def load_annotation(annotation_path):
    print(f"\n  Loading probe annotation from: {annotation_path}")
    
    # 1. Find the exact row where the data table starts
    skip_rows = 0
    with open(annotation_path, 'r', encoding='utf-8', errors='ignore') as f:
        for i, line in enumerate(f):
            # We look for the line that defines the columns
            if line.startswith("ID\t") or line.startswith('"ID"\t'):
                skip_rows = i
                break
    
    # 2. Read the table using Tab separator
    df = pd.read_csv(annotation_path, sep="\t", skiprows=skip_rows, low_memory=False)

    # 3. Clean column names (remove quotes/spaces)
    df.columns = [c.strip().strip('"') for c in df.columns]
    
    # 4. Set the specific columns found in your file
    id_col = 'ID'
    sym_col = 'Gene symbol'  # Matched to your file's lowercase 's'

    print(f"  ✓ Found columns: {id_col} and {sym_col}")

    probe_to_gene = {}
    no_symbol = 0

    # 5. Build the dictionary
    for _, row in df.iterrows():
        probe_id = str(row[id_col]).strip().strip('"')
        gene = str(row.get(sym_col, "")).strip()

        if gene and gene not in ("nan", "---", "", "null", "None"):
            # If multiple genes (e.g. "BDNF///NTRK2"), take the first one
            gene = gene.split("///")[0].split(";")[0].strip()
            probe_to_gene[probe_id] = gene
        else:
            no_symbol += 1

    print(f"  ✓ Mapped {len(probe_to_gene):,} probes to gene symbols")
    print(f"  ⚠ {no_symbol:,} probes have no gene symbol (will be skipped)")
    return probe_to_gene


# ── STEP 3 ─────────────────────────────────────────────────────────────────────
def get_sample_columns(df):
    """
    Return the list of sample (GSM) columns from the expression dataframe.
    In a series matrix file the first col is 'ID_REF', rest are GSM IDs.
    """
    id_col      = df.columns[0]   # always 'ID_REF'
    sample_cols = list(df.columns[1:])
    print(f"  ID column: '{id_col}'")
    print(f"  {len(sample_cols)} sample columns found")
    print(f"  Example sample IDs: {sample_cols[:3]}")
    return id_col, sample_cols


def ensure_samples_exist(sample_cols, conn):
    """Add any missing samples to the samples table."""
    cursor = conn.cursor()
    sql = ("INSERT IGNORE INTO samples (sample_id, `condition`, tissue, dataset) "
           "VALUES (%s, %s, %s, %s)")
    rows = [(s, "unknown", "prefrontal cortex", "GSE33000") for s in sample_cols]
    cursor.executemany(sql, rows)
    conn.commit()
    cursor.close()
    print(f"  ✓ Ensured {len(rows)} samples exist in samples table")


def import_expression(df, id_col, sample_cols, probe_to_gene, conn):
    """
    For each probe with a known gene symbol, insert expression values
    into rna_expression (one row per gene × sample).
    """
    cursor = conn.cursor()
    sql = """
        INSERT IGNORE INTO rna_expression (gene_symbol, sample_id, expression_value)
        VALUES (%s, %s, %s)
    """

    batch        = []
    batch_size   = 5000
    total        = 0
    skipped_probes = 0

    for _, row in df.iterrows():
        probe_id = str(row[id_col]).strip().strip('"')
        gene     = probe_to_gene.get(probe_id)

        if not gene:
            skipped_probes += 1
            continue

        for sample_col in sample_cols:
            val = row[sample_col]
            try:
                val = float(val)
            except (ValueError, TypeError):
                continue
            batch.append((gene, str(sample_col).strip(), val))

        if len(batch) >= batch_size:
            cursor.executemany(sql, batch)
            conn.commit()
            total += len(batch)
            batch = []
            print(f"    ... {total:,} rows inserted", end="\r")

    if batch:
        cursor.executemany(sql, batch)
        conn.commit()
        total += len(batch)

    cursor.close()
    print(f"\n  ✓ rna_expression: {total:,} rows inserted")
    print(f"  ⚠ {skipped_probes:,} probes skipped (no gene symbol in annotation)")


# ── MAIN ────────────────────────────────────────────────────────────────────────
def connect():
    conn = mysql.connector.connect(**DB_CONFIG)
    print("✓ Connected to database")
    return conn


def main():
    parser = argparse.ArgumentParser(
        description="Import GSE33000 (GPL10558 Illumina HT-12 v4) into brain_multiomics"
    )
    parser.add_argument(
        "--file", required=True,
        help="Path to GSE33000_series_matrix.txt (or .txt.gz)"
    )
    parser.add_argument(
        "--annotation", required=True,
        help="Path to GPL10558 annotation CSV (see script header for how to download)"
    )
    args = parser.parse_args()

    print("\n=== Importing GSE33000 microarray expression ===\n")

    # 1. Parse the series matrix — skipping header metadata correctly
    df = parse_series_matrix(args.file)

    # 2. Load probe → gene mapping
    probe_to_gene = load_annotation(args.annotation)

    # 3. Identify sample columns
    id_col, sample_cols = get_sample_columns(df)

    # 4. Connect and import
    conn = connect()
    ensure_samples_exist(sample_cols, conn)
    import_expression(df, id_col, sample_cols, probe_to_gene, conn)
    conn.close()

    print("\n=== Done ===\n")
    print("Quick check — run in MySQL:")
    print("  SELECT gene_symbol, COUNT(*) FROM rna_expression")
    print("  WHERE gene_symbol='BDNF' GROUP BY gene_symbol;")


if __name__ == "__main__":
    main()
