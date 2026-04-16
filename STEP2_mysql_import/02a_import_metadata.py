"""
02a_import_metadata.py
======================
Imports clean_metadata_full.csv into the `samples` AND `metadata_import` tables.

The `samples` table is what all other tables link to (foreign key).
The `metadata_import` table keeps the full original metadata for reference.

Usage:

    python 02a_import_metadata.py --file /home/daria/Documents/UGent/databases/multiomics/chatgpt_version/multiomics_data/STEP1_download-meta-data/clean_metadata_full.csv

Requirements:
    pip install mysql-connector-python pandas
"""

import argparse
import pandas as pd
import mysql.connector
import sys

# ── Database connection settings ──────────────────────────────────────────────
DB_CONFIG = {
    "host":     "localhost",
    "user":     "daria",          
    "password": "simba", 
    "database": "brain_multiomics"
}

# ── Column name mapping ────────────────────────────────────────────────────────
# Maps column names in YOUR csv → what the database expects.
# Edit the left side to match your actual CSV column headers.
COLUMN_MAP = {
    "sample_id": "sample_id",
    "dataset": "dataset",
    "organism_ch1": "organism",
    "source_name_ch1": "tissue",
    "sex": "sex",
    "age": "age",
    "condition": "condition",
}

DB_TO_CSV = {db_col: csv_col for csv_col, db_col in COLUMN_MAP.items()}

def connect():
    """Create and return a database connection."""
    try:
        conn = mysql.connector.connect(**DB_CONFIG)
        print("✓ Connected to MySQL database")
        return conn
    except mysql.connector.Error as e:
        print(f"✗ Could not connect to MySQL: {e}")
        sys.exit(1)

def load_csv(filepath):
    """Load the metadata CSV into a pandas DataFrame."""
    print(f"  Loading {filepath} ...")
    df = pd.read_csv(filepath)
    print(f"  Found {len(df)} rows and columns: {list(df.columns)}")
    return df

def normalise_condition(value):

    if pd.isna(value):
        return "unknown"

    v = str(value).lower().strip()

    if "ad" in v or "alzheimer" in v:
        return "AD"

    if "control" in v:
        return "control"

    if "exclude" in v:
        return "exclude"

    return v

def clean_str(value, default=""):
    if pd.isna(value):
        return default
    return str(value).strip()

def get_mapped_value(row, db_col, default=""):
    csv_col = DB_TO_CSV.get(db_col, db_col)
    return clean_str(row.get(csv_col, default), default=default)

def get_table_columns(conn, table_name):
    cursor = conn.cursor()
    cursor.execute(
        """
        SELECT COLUMN_NAME
        FROM INFORMATION_SCHEMA.COLUMNS
        WHERE TABLE_SCHEMA = %s AND TABLE_NAME = %s
        """,
        (DB_CONFIG["database"], table_name),
    )
    cols = {r[0] for r in cursor.fetchall()}
    cursor.close()
    return cols

def import_metadata(df, conn):
    """Insert rows into metadata_import and samples tables."""
    cursor = conn.cursor()

    inserted_metadata = 0
    inserted_samples  = 0
    skipped           = 0

    meta_existing = get_table_columns(conn, "metadata_import")
    samples_existing = get_table_columns(conn, "samples")

    meta_desired = ["sample_id", "dataset", "organism", "condition", "tissue", "sex", "age"]
    meta_cols = [c for c in meta_desired if c in meta_existing]
    if not meta_cols:
        raise RuntimeError("metadata_import has no expected columns. Check table schema.")

    sample_desired = ["sample_id", "condition", "tissue", "dataset"]
    sample_cols = [c for c in sample_desired if c in samples_existing]
    if not sample_cols:
        raise RuntimeError("samples has no expected columns. Check table schema.")

    for _, row in df.iterrows():
        sample_id = str(row.get("sample_id", "")).strip()
        if not sample_id:
            skipped += 1
            continue

        # ── Insert into metadata_import (full record) ────────────────────────
        meta_sql = f"""
            INSERT IGNORE INTO metadata_import
                ({", ".join(f"`{c}`" for c in meta_cols)})
            VALUES ({", ".join(["%s"] * len(meta_cols))})
        """
        meta_vals_map = {
            "sample_id": sample_id,
            "dataset": get_mapped_value(row, "dataset", ""),
            "organism": get_mapped_value(row, "organism", ""),
            "condition": get_mapped_value(row, "condition", ""),
            "tissue": get_mapped_value(row, "tissue", ""),
            "sex": get_mapped_value(row, "sex", ""),
            "age": get_mapped_value(row, "age", ""),
        }
        meta_vals = tuple(meta_vals_map[c] for c in meta_cols)
        cursor.execute(meta_sql, meta_vals)
        inserted_metadata += cursor.rowcount

        # ── Insert into samples (simplified, linked table) ───────────────────
        condition = normalise_condition(get_mapped_value(row, "condition", ""))
        tissue    = get_mapped_value(row, "tissue", "prefrontal cortex")
        dataset   = get_mapped_value(row, "dataset", "")

        sample_vals_map = {
            "sample_id": sample_id,
            "condition": condition,
            "tissue": tissue,
            "dataset": dataset,
        }
        sample_sql = f"""
            INSERT IGNORE INTO samples ({", ".join(f"`{c}`" for c in sample_cols)})
            VALUES ({", ".join(["%s"] * len(sample_cols))})
        """
        cursor.execute(sample_sql, tuple(sample_vals_map[c] for c in sample_cols))
        inserted_samples += cursor.rowcount

    conn.commit()
    cursor.close()

    print(f"  ✓ metadata_import: {inserted_metadata} rows inserted")
    print(f"  ✓ samples:         {inserted_samples} rows inserted")
    print(f"  ⚠ skipped:         {skipped} rows (missing sample_id)")

def main():
    parser = argparse.ArgumentParser(description="Import metadata CSV into brain_multiomics")
    parser.add_argument("--file", required=True, help="Path to clean_metadata_full.csv")
    args = parser.parse_args()

    print("\n=== Importing metadata ===")
    df   = load_csv(args.file)
    conn = connect()
    import_metadata(df, conn)
    conn.close()
    print("=== Done ===\n")

if __name__ == "__main__":
    main()
