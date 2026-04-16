"""
02d_import_methylation_SMART.py
================================
Solves the "table is full" problem by only importing CpGs that are
relevant to your analysis — i.e. CpGs that map to plasticity genes.

WHY 257M ROWS IS TOO MANY:
  485,000 CpGs × 531 samples = 257M rows. Each row is ~80 bytes.
  That's ~20GB just for methylation. MySQL's default InnoDB tablespace
  limit is often 64GB but on a VM with limited disk it fills up fast.
  More importantly: you don't need all 485,000 CpGs — you only need
  the ones near your plasticity genes.

TWO STRATEGIES (pick one):

  Strategy A — plasticity genes only (RECOMMENDED, ~50,000–100,000 rows)
    Filter to CpGs that map to your plasticity_genes list via the
    Illumina 450k annotation file. Only ~200-500 CpGs per gene.
    Run time: minutes. Table size: manageable.

  Strategy B — per-sample averages (medium, ~250,000 rows)
    Instead of storing every individual sample, store the mean beta
    value per CpG per condition (AD vs control). You lose per-sample
    variance but the table stays small.
    Use --strategy average

  Strategy C — full import to a different storage engine
    If you really need all 257M rows, switch the table to MyISAM or
    partition it. See the SQL commands printed by --strategy full.

Usage:
  # Download annotation first (if not done yet):
  wget "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13534/annot/GPL13534.annot.gz"
  gunzip GPL13534.annot.gz

  # Strategy A — only plasticity CpGs (recommended):
  python 02d_import_methylation_SMART.py \
    --file GSE59685_betas.csv \
    --annotation GPL13534.annot \
    --strategy plasticity

  # Strategy B — per-condition averages:
  python 02d_import_methylation_SMART.py \
    --file GSE59685_betas.csv \
    --annotation GPL13534.annot \
    --strategy average

  # Test any strategy with 5000 CpG rows first:
  python 02d_import_methylation_SMART.py ... --limit 5000
"""

import os
import argparse
import mysql.connector
import csv
import gzip
import sys
from collections import defaultdict

DB_CONFIG = {
    "host":     "localhost",
    "user":     "daria",
    "password": os.getenv("DB_PASSWORD", ""),
    "database": "brain_multiomics"
}

def connect():
    return mysql.connector.connect(**DB_CONFIG)

def open_file(filepath):
    if filepath.endswith(".gz"):
        return gzip.open(filepath, "rt", encoding="utf-8", errors="replace")
    return open(filepath, "r", encoding="utf-8", errors="replace")

# ── STEP 1: load plasticity CpG list from annotation ─────────────────────────
def get_plasticity_cpgs(annotation_path, conn):
    """
    Returns a set of cpg_ids (e.g. {'cg00004067', ...}) that map to
    any gene in your plasticity_genes table.
    Also returns a dict: cpg_id → (gene_symbol, relation)
    """
    print("  Loading plasticity genes from database...")
    cursor = conn.cursor()
    cursor.execute("SELECT gene_symbol FROM plasticity_genes")
    plasticity_genes = {row[0].upper() for row in cursor.fetchall()}
    cursor.close()
    print(f"  Found {len(plasticity_genes)} plasticity genes")

    if not plasticity_genes:
        print("  ⚠ plasticity_genes table is empty!")
        print("    Run 04_filter_plasticity_genes.py first.")
        sys.exit(1)

    print(f"  Loading 450k annotation from {annotation_path}...")
    plasticity_cpgs   = set()
    cpg_gene_map      = {}   # cpg_id → (gene_symbol, relation)
    total_probes      = 0
    matched_probes    = 0

    # Detect separator
    with open_file(annotation_path) as fh:
        for line in fh:
            if not line.startswith("#") and line.strip():
                sep = "\t" if "\t" in line else ","
                break

    with open_file(annotation_path) as fh:
        reader = csv.DictReader(
            (l for l in fh if not l.startswith("#")),
            delimiter=sep
        )
        # Normalise column names
        first_row = None
        for row in reader:
            if first_row is None:
                first_row = row
                # Find the right column names
                cols = list(row.keys())
                id_col   = next((c for c in cols if c.strip().lower() in
                                 ("id", "ilmnid", "probe_id", "name")), None)
                gene_col = next((c for c in cols if "refgene_name" in c.lower()
                                 or c.strip().lower() in ("symbol","gene_symbol")), None)
                grp_col  = next((c for c in cols if "refgene_group" in c.lower()), None)
                if not id_col or not gene_col:
                    print(f"  ✗ Cannot find ID/gene columns. Columns: {cols[:8]}")
                    sys.exit(1)
                print(f"  Using columns: id='{id_col}', gene='{gene_col}', "
                      f"region='{grp_col}'")

            total_probes += 1
            cpg_id   = row[id_col].strip()
            gene_str = row[gene_col].strip() if gene_col else ""
            grp_str  = row[grp_col].strip()  if grp_col  else ""

            if not gene_str or gene_str in ("nan", ""):
                continue

            # Probes can map to multiple genes (semicolon-separated)
            genes    = [g.strip().upper() for g in gene_str.split(";") if g.strip()]
            groups   = [g.strip() for g in grp_str.split(";")]  if grp_str else []

            for i, gene in enumerate(genes):
                if gene in plasticity_genes:
                    relation = groups[i] if i < len(groups) else "unknown"
                    plasticity_cpgs.add(cpg_id)
                    # Store first plasticity gene match per CpG
                    if cpg_id not in cpg_gene_map:
                        cpg_gene_map[cpg_id] = (gene.title(), relation)
                    matched_probes += 1
                    break

    print(f"  Scanned {total_probes:,} probes → "
          f"{len(plasticity_cpgs):,} map to plasticity genes")
    return plasticity_cpgs, cpg_gene_map

# ── STEP 2: read header ───────────────────────────────────────────────────────
def read_header(filepath):
    gsm_ids = []
    with open_file(filepath) as fh:
        for i, line in enumerate(fh):
            if i == 5:  # line 6 = new GSMs
                cols    = line.rstrip("\n").split(",")
                gsm_ids = [c.strip() for c in cols[1:] if c.strip().startswith("GSM")]
                break
    print(f"  Found {len(gsm_ids)} GSM sample IDs")
    return gsm_ids

# ── STEP 3a: import plasticity CpGs only ─────────────────────────────────────
def import_plasticity_only(filepath, gsm_ids, plasticity_cpgs,
                           cpg_gene_map, conn, limit=None):
    """
    Stream the file and only INSERT rows where cpg_id is in plasticity_cpgs.
    Estimated rows: ~500 CpGs × 531 samples = ~265,000 rows — tiny.
    """
    cursor   = conn.cursor()
    insert_sql = """
        INSERT IGNORE INTO methylation (cpg_id, beta_value, sample_id, tissue)
        VALUES (%s, %s, %s, %s)
    """
    # Also populate methylation_gene_links while we're at it
    link_sql = """
        INSERT IGNORE INTO methylation_gene_links (cpg_id, gene_symbol, relation)
        VALUES (%s, %s, %s)
    """

    batch        = []
    link_batch   = []
    inserted     = 0
    cpgs_seen    = set()
    cpg_rows     = 0
    BATCH_SIZE   = 3000

    print(f"\n  Importing {len(plasticity_cpgs):,} plasticity CpGs × "
          f"{len(gsm_ids)} samples...")

    with open_file(filepath) as fh:
        reader = csv.reader(fh)
        for line_num, row in enumerate(reader):
            if line_num < 7:
                continue
            if not row or not row[0].startswith("cg"):
                continue

            cpg_id = row[0].strip()
            if cpg_id not in plasticity_cpgs:
                continue

            cpg_rows += 1

            for j, gsm in enumerate(gsm_ids):
                col_idx = j + 1
                if col_idx >= len(row):
                    continue
                try:
                    beta = float(row[col_idx])
                except (ValueError, TypeError):
                    continue
                batch.append((cpg_id, beta, gsm, "prefrontal cortex"))

            # Add gene link if not done yet
            if cpg_id not in cpgs_seen and cpg_id in cpg_gene_map:
                gene, relation = cpg_gene_map[cpg_id]
                link_batch.append((cpg_id, gene, relation))
                cpgs_seen.add(cpg_id)

            if len(batch) >= BATCH_SIZE:
                cursor.executemany(insert_sql, batch)
                if link_batch:
                    cursor.executemany(link_sql, link_batch)
                conn.commit()
                inserted   += len(batch)
                batch       = []
                link_batch  = []

            if cpg_rows % 50 == 0:
                print(f"  Plasticity CpGs found: {cpg_rows:>5,}  |  "
                      f"rows inserted: {inserted:>8,}", end="\r")

            if limit and cpg_rows >= limit:
                break

    if batch:
        cursor.executemany(insert_sql, batch)
        if link_batch:
            cursor.executemany(link_sql, link_batch)
        conn.commit()
        inserted += len(batch)

    cursor.close()
    print(f"\n  ✓ Plasticity CpGs imported : {cpg_rows:,}")
    print(f"  ✓ DB rows inserted         : {inserted:,}")
    print(f"  ✓ Gene links added         : {len(cpgs_seen):,}")

# ── STEP 3b: per-condition averages ──────────────────────────────────────────
def import_averages(filepath, gsm_ids, plasticity_cpgs, conn, limit=None):
    """
    Compute mean beta per CpG per condition (AD / control / unknown)
    and store one row per CpG per condition instead of per sample.
    Requires a methylation_averages table (created here if not exists).
    """
    cursor = conn.cursor()
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS methylation_averages (
            cpg_id      VARCHAR(50),
            condition   VARCHAR(50),
            mean_beta   FLOAT,
            n_samples   INT,
            PRIMARY KEY (cpg_id, `condition`)
        )
    """)
    conn.commit()

    # Get condition per sample
    cursor.execute("SELECT sample_id, `condition` FROM samples WHERE dataset='GSE59685'")
    sample_condition = {row[0]: row[1] for row in cursor.fetchall()}

    print(f"\n  Computing per-condition averages for "
          f"{len(plasticity_cpgs) if plasticity_cpgs else 'ALL'} CpGs...")

    cpg_rows = 0
    inserted = 0
    BATCH_SIZE = 2000
    batch    = []

    with open_file(filepath) as fh:
        reader = csv.reader(fh)
        for line_num, row in enumerate(reader):
            if line_num < 7:
                continue
            if not row or not row[0].startswith("cg"):
                continue

            cpg_id = row[0].strip()
            if plasticity_cpgs and cpg_id not in plasticity_cpgs:
                continue

            # Accumulate by condition
            by_cond = defaultdict(list)
            for j, gsm in enumerate(gsm_ids):
                col_idx = j + 1
                if col_idx >= len(row):
                    continue
                try:
                    beta = float(row[col_idx])
                except (ValueError, TypeError):
                    continue
                cond = sample_condition.get(gsm, "unknown")
                by_cond[cond].append(beta)

            for cond, vals in by_cond.items():
                mean_beta = sum(vals) / len(vals)
                batch.append((cpg_id, cond, mean_beta, len(vals)))

            cpg_rows += 1

            if len(batch) >= BATCH_SIZE:
                cursor.executemany("""
                    INSERT INTO methylation_averages
                        (cpg_id, `condition`, mean_beta, n_samples)
                    VALUES (%s, %s, %s, %s)
                    ON DUPLICATE KEY UPDATE mean_beta=VALUES(mean_beta),
                                            n_samples=VALUES(n_samples)
                """, batch)
                conn.commit()
                inserted += len(batch)
                batch     = []

            if cpg_rows % 200 == 0:
                print(f"  CpGs processed: {cpg_rows:>6,}  |  "
                      f"rows inserted: {inserted:>8,}", end="\r")

            if limit and cpg_rows >= limit:
                break

    if batch:
        cursor.executemany("""
            INSERT INTO methylation_averages
                (cpg_id, `condition`, mean_beta, n_samples)
            VALUES (%s, %s, %s, %s)
            ON DUPLICATE KEY UPDATE mean_beta=VALUES(mean_beta),
                                    n_samples=VALUES(n_samples)
        """, batch)
        conn.commit()
        inserted += len(batch)

    cursor.close()
    print(f"\n  ✓ CpGs processed  : {cpg_rows:,}")
    print(f"  ✓ Rows inserted   : {inserted:,}")
    print("  Table: methylation_averages")

# ── MAIN ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file",       required=True)
    parser.add_argument("--annotation", required=True,
                        help="Illumina 450k annotation file (GPL13534_annotation.csv)")
    parser.add_argument("--strategy",   default="plasticity",
                        choices=["plasticity", "average"],
                        help="plasticity=only plasticity CpGs (default); "
                             "average=per-condition means for plasticity CpGs")
    parser.add_argument("--limit",      type=int, default=None,
                        help="Stop after N CpG rows (for testing)")
    args = parser.parse_args()

    print("\n=== Importing methylation — smart mode ===\n")

    conn    = connect()
    gsm_ids = read_header(args.file)

    # Always filter to plasticity CpGs to keep the table small
    plasticity_cpgs, cpg_gene_map = get_plasticity_cpgs(args.annotation, conn)

    if args.strategy == "plasticity":
        # Truncate existing methylation table first to avoid mixing with
        # the 56M rows already inserted from the failed full import
        print("\n  Clearing existing methylation table...")
        cursor = conn.cursor()
        cursor.execute("TRUNCATE TABLE methylation")
        conn.commit()
        cursor.close()
        print("  ✓ Table cleared\n")

        import_plasticity_only(
            args.file, gsm_ids, plasticity_cpgs,
            cpg_gene_map, conn, limit=args.limit
        )

    elif args.strategy == "average":
        import_averages(
            args.file, gsm_ids, plasticity_cpgs,
            conn, limit=args.limit
        )

    conn.close()
    print("\n=== Done ===\n")
    print("Quick check:")
    print("  SELECT COUNT(*) FROM methylation;")
    print("  SELECT cpg_id, sample_id, beta_value FROM methylation LIMIT 5;")

if __name__ == "__main__":
    main()