import os
import argparse
import pandas as pd
import mysql.connector
import sys
import io

DB_CONFIG = {
    "host":     "localhost",
    "user":     "daria",
    "password": os.getenv("DB_PASSWORD", ""), 
    "database": "brain_multiomics"
}

def clean_id(val):
    """Force value to string and strip quotes/spaces to ensure matching works."""
    return str(val).strip().strip('"').strip("'")

# ── STEP 1: PARSE THE SERIES MATRIX ──────────────────────────────────────────
def parse_series_matrix(filepath):
    print(f"  Parsing series matrix file: {filepath}")
    data_lines = []
    in_table   = False

    with open(filepath, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                data_lines.append(line)

    if not data_lines:
        print("  ✗ ERROR: Could not find '!series_matrix_table_begin'")
        sys.exit(1)

    content = "".join(data_lines)
    # Read as dtype=str to prevent pandas from mangling long ID numbers
    df = pd.read_csv(io.StringIO(content), sep="\t", low_memory=False, dtype=str)
    return df

# ── STEP 2: LOAD ANNOTATION ──────────────────────────────────────────────────
def load_annotation(annotation_path):
    print(f"\n  Loading probe annotation from: {annotation_path}")
    
    skip_rows = 0
    with open(annotation_path, 'r', encoding='utf-8', errors='ignore') as f:
        for i, line in enumerate(f):
            if line.startswith("ID\t") or line.startswith('"ID"\t'):
                skip_rows = i
                break
    
    df = pd.read_csv(annotation_path, sep="\t", skiprows=skip_rows, low_memory=False, dtype=str)
    df.columns = [c.strip().strip('"') for c in df.columns]
    
    id_col = 'ID'
    sym_options = ['Gene symbol', 'Gene Symbol', 'Symbol']
    sym_col = next((c for c in sym_options if c in df.columns), None)

    if not sym_col:
        print(f"  ✗ Error: Could not find Gene Symbol column. Found: {list(df.columns)}")
        sys.exit(1)

    probe_to_gene = {}
    for _, row in df.iterrows():
        # CLEAN THE ID
        p_id = clean_id(row[id_col])
        gene = str(row.get(sym_col, "")).strip()

        if gene and gene not in ("nan", "---", "", "null", "None"):
            gene = gene.replace("///", ";").split(";")[0].strip()
            probe_to_gene[p_id] = gene

    print(f"  ✓ Mapped {len(probe_to_gene):,} probes to gene symbols")
    return probe_to_gene

# ── STEP 3: IMPORT TO MYSQL ──────────────────────────────────────────────────
def import_expression(df, probe_to_gene):
    conn = mysql.connector.connect(**DB_CONFIG)
    cursor = conn.cursor()
    
    # In a series matrix, column 0 is 'ID_REF', others are GSM IDs
    id_col_name = df.columns[0] 
    sample_cols = list(df.columns[1:])
    
    print(f"  Ensuring {len(sample_cols)} samples exist in database...")
    sql_samp = "INSERT IGNORE INTO samples (sample_id, dataset, tissue) VALUES (%s, %s, %s)"
    # Clean the GSM IDs too just in case they have quotes
    cursor.executemany(sql_samp, [(clean_id(s), "GSE33000", "prefrontal cortex") for s in sample_cols])
    conn.commit()

    sql_expr = "INSERT IGNORE INTO rna_expression (gene_symbol, sample_id, expression_value) VALUES (%s, %s, %s)"
    
    batch = []
    total = 0
    skipped_probes = 0

    print("  Starting data import...")
    for _, row in df.iterrows():
        # CLEAN THE ID from the matrix file to match the annotation map
        matrix_probe_id = clean_id(row[id_col_name])
        gene = probe_to_gene.get(matrix_probe_id)

        if not gene:
            skipped_probes += 1
            continue

        for sample_id in sample_cols:
            val = row[sample_id]
            try:
                # Clean sample_id to match the GSM in the samples table
                gsm_id = clean_id(sample_id)
                batch.append((gene, gsm_id, float(val)))
            except:
                continue

        if len(batch) >= 10000:
            cursor.executemany(sql_expr, batch)
            conn.commit()
            total += len(batch)
            print(f"    ... {total:,} rows inserted", end="\r")
            batch = []

    if batch:
        cursor.executemany(sql_expr, batch)
        conn.commit()
        total += len(batch)

    cursor.close()
    conn.close()
    print(f"\n  ✓ DONE: {total:,} rows inserted! ({skipped_probes:,} probes skipped)")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", required=True)
    parser.add_argument("--annotation", required=True)
    args = parser.parse_args()

    matrix_df = parse_series_matrix(args.file)
    mapping = load_annotation(args.annotation)
    import_expression(matrix_df, mapping)

if __name__ == "__main__":
    main()