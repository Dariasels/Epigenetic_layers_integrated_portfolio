import pandas as pd
import mysql.connector
import sys
import argparse
import math
import re

DB_CONFIG = {
    "host": "localhost",
    "user": "daria",
    "password": "simba", 
    "database": "brain_multiomics"
}

def parse_and_import(file_path, full_run=False):
    print(f"Reading {file_path}...")
    try:
        df = pd.read_csv(file_path, sep=None, engine='python', header=None)
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    # 1. Determine Sample Names
    first_row = df.iloc[0]
    # Check if first row contains strings (likely headers)
    if isinstance(first_row[2], str) and not str(first_row[2]).replace('.','').isdigit():
        df.columns = df.iloc[0]
        df = df[1:].reset_index(drop=True)
    else:
        # If no headers, create 'Sample_0', 'Sample_1' etc.
        df.columns = ['ID', 'Gene'] + [f'Sample_{i}' for i in range(len(df.columns)-2)]

    gene_col_idx = 1
    sample_cols = df.columns[2:] 

    conn = mysql.connector.connect(**DB_CONFIG)
    cursor = conn.cursor()

    # --- CRITICAL FIX: Add Samples to the 'samples' table first ---
    print(f"Ensuring {len(sample_cols)} samples are registered...")
    sql_samples = "INSERT IGNORE INTO samples (sample_id, dataset, tissue) VALUES (%s, %s, %s)"
    sample_batch = [(str(s), "GSE33000", "prefrontal cortex") for s in sample_cols]
    cursor.executemany(sql_samples, sample_batch)
    conn.commit()
    # --------------------------------------------------------------

    sql_expr = "INSERT IGNORE INTO rna_expression (gene_symbol, sample_id, expression_value) VALUES (%s, %s, %s)"
    
    batch = []
    inserted_count = 0
    skipped_rows = 0

    process_rows = df if full_run else df.head(500)
    print("Starting import...")
    for _, row in process_rows.iterrows():
        raw_gene = str(row.iloc[gene_col_idx]).strip()
        
        # Aggressive filter
        if (raw_gene.upper().startswith("RSE_") or 
            raw_gene.upper().startswith("LOC") or 
            raw_gene.isdigit() or 
            len(raw_gene) < 2 or 
            raw_gene.lower() in ("nan", "null", "", "symbol")):
            skipped_rows += 1
            continue
            
        gene_symbol = re.split(r'[,\s;]+', raw_gene)[0].strip()

        for sample_id in sample_cols:
            val = row[sample_id]
            if pd.isna(val): continue
            try:
                f_val = float(val)
                if math.isfinite(f_val):
                    batch.append((gene_symbol, str(sample_id), f_val))
            except:
                continue

            if len(batch) >= 10000:
                cursor.executemany(sql_expr, batch)
                conn.commit()
                inserted_count += len(batch)
                batch = []
                print(f"  Inserted {inserted_count:,} records...", end="\r")

    if batch:
        cursor.executemany(sql_expr, batch)
        conn.commit()
        inserted_count += len(batch)

    cursor.close()
    conn.close()
    print(f"\nDONE! Total inserted: {inserted_count:,}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", required=True)
    parser.add_argument("--full", action="store_true")
    args = parser.parse_args()
    parse_and_import(args.file, full_run=args.full)