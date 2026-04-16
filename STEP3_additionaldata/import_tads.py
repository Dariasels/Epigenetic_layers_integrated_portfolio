"""
import_tads.py
==============
Imports TAD (Topologically Associating Domain) regions from a BED file
into the tads table, then maps plasticity genes and ATAC peaks into TADs.

About TADs:
  A TAD is a 3D chromosomal neighbourhood — typically 0.1 to 3 Mb in size.
  DNA regions INSIDE the same TAD physically interact much more with each
  other than with regions OUTSIDE. This means:
    - An enhancer and a gene in the SAME TAD can regulate each other
    - An enhancer and a gene in DIFFERENT TADs typically cannot
  TAD boundaries (marked by CTCF protein binding) act like insulators.

  Your file: GSE105194_ENCFF306YQN_topologically_associated_domains_hg19.bed
  Source: ENCODE project, hg19
  Format: BED3 (chrom, start, end) — possibly BED6 with score/strand

Why use this even though it's not brain-specific?
  TAD boundaries are largely conserved across human cell types.
  Using ENCODE TADs as a reference is standard in AD epigenomics papers.

Usage:
  python import_tads.py \
    --file ~/Documents/UGent/databases/multiomics/chatgpt_version/\
multiomics_data/extra_data_notthisproj/\
GSE105194_ENCFF306YQN_topologically_associated_domains_hg19.bed

  Then run the mapping:
  python import_tads.py --map-genes
  python import_tads.py --map-atac
"""

import os
import argparse
import mysql.connector
import sys

DB_CONFIG = {
    "host":     "localhost",
    "user":     "daria",
    "password": os.getenv("DB_PASSWORD", ""),
    "database": "brain_multiomics"
}

def connect():
    return mysql.connector.connect(**DB_CONFIG)

# ── IMPORT TAD REGIONS ────────────────────────────────────────────────────────
def import_tads(filepath, conn):
    cursor = conn.cursor()

    sql = """
        INSERT INTO tads (chrom, start_pos, end_pos, tad_size_kb, source)
        VALUES (%s, %s, %s, %s, %s)
    """

    batch    = []
    inserted = 0
    skipped  = 0
    BATCH    = 2000
    SOURCE   = "GSE105194 ENCODE hg19"

    with open(filepath, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue

            cols = line.split("\t")
            if len(cols) < 3:
                continue

            chrom = cols[0].strip()
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom

            chrom_num = chrom.replace("chr", "").split("_")[0]
            if chrom_num not in [str(i) for i in range(1, 23)] + ["X", "Y"]:
                skipped += 1
                continue

            try:
                start = int(cols[1])
                end   = int(cols[2])
            except ValueError:
                skipped += 1
                continue

            size_kb = (end - start) // 1000
            batch.append((chrom, start, end, size_kb, SOURCE))

            if len(batch) >= BATCH:
                cursor.executemany(sql, batch)
                conn.commit()
                inserted += len(batch)
                batch = []
                print(f"  {inserted:,} TADs inserted...", end="\r")

    if batch:
        cursor.executemany(sql, batch)
        conn.commit()
        inserted += len(batch)

    cursor.close()
    print(f"\n  ✓ {inserted:,} TADs inserted  ({skipped} skipped — non-standard chroms)")

    # Print size distribution
    cursor = conn.cursor()
    cursor.execute("""
        SELECT
            SUM(tad_size_kb < 100)   AS tiny_under100kb,
            SUM(tad_size_kb BETWEEN 100 AND 500) AS small_100_500kb,
            SUM(tad_size_kb BETWEEN 500 AND 1000) AS medium_500kb_1mb,
            SUM(tad_size_kb > 1000)  AS large_over1mb
        FROM tads
    """)
    row = cursor.fetchone()
    print(f"\n  TAD size distribution:")
    print(f"    < 100kb   : {row[0]:,}")
    print(f"    100–500kb : {row[1]:,}")
    print(f"    500kb–1Mb : {row[2]:,}")
    print(f"    > 1Mb     : {row[3]:,}")
    cursor.close()

# ── MAP GENES TO TADs ─────────────────────────────────────────────────────────
def map_genes_to_tads(conn):
    """
    For each gene in gene_coordinates, find which TAD it falls in.
    Uses the gene's TSS (start_pos) as the anchor point.
    A gene's TSS must fall WITHIN the TAD boundaries.
    """
    print("\n  Mapping genes to TADs (using TSS position)...")
    cursor = conn.cursor()

    # Clear existing links first
    cursor.execute("TRUNCATE TABLE tad_gene_links")
    conn.commit()

    sql_insert = """
        INSERT INTO tad_gene_links (tad_id, gene_symbol)
        SELECT t.tad_id, gc.gene_symbol
        FROM gene_coordinates gc
        JOIN tads t ON t.chrom = gc.chrom
                   AND gc.start_pos BETWEEN t.start_pos AND t.end_pos
    """

    cursor.execute(sql_insert)
    inserted = cursor.rowcount
    conn.commit()
    cursor.close()

    print(f"  ✓ {inserted:,} gene-TAD links created")

    # How many plasticity genes are covered?
    cursor = conn.cursor()
    cursor.execute("""
        SELECT COUNT(DISTINCT tgl.gene_symbol)
        FROM tad_gene_links tgl
        JOIN plasticity_genes pg ON pg.gene_symbol = tgl.gene_symbol
    """)
    n_plasticity = cursor.fetchone()[0]
    cursor.close()
    print(f"  ✓ {n_plasticity} plasticity genes assigned to a TAD")

# ── MAP ATAC PEAKS TO TADs ────────────────────────────────────────────────────
def map_atac_to_tads(conn):
    """
    For each ATAC peak, find which TAD it overlaps.
    A peak overlaps a TAD if any part of it falls within the TAD.
    """
    print("\n  Mapping ATAC peaks to TADs...")
    cursor = conn.cursor()

    cursor.execute("TRUNCATE TABLE tad_atac_links")
    conn.commit()

    sql_insert = """
        INSERT INTO tad_atac_links (tad_id, peak_id)
        SELECT t.tad_id, ap.peak_id
        FROM atac_peaks ap
        JOIN tads t ON t.chrom = ap.chrom
                   AND ap.chrom_start < t.end_pos
                   AND ap.chrom_end   > t.start_pos
    """

    cursor.execute(sql_insert)
    inserted = cursor.rowcount
    conn.commit()
    cursor.close()
    print(f"  ✓ {inserted:,} ATAC peak-TAD links created")

# ── MAIN ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file",       help="Path to TAD BED file")
    parser.add_argument("--map-genes",  action="store_true",
                        help="Map genes to TADs (run after import)")
    parser.add_argument("--map-atac",   action="store_true",
                        help="Map ATAC peaks to TADs (run after import)")
    args = parser.parse_args()

    conn = connect()
    print("✓ Connected to database")

    if args.file:
        print(f"\n=== Importing TADs from {args.file} ===")
        import_tads(args.file, conn)

    if args.map_genes:
        map_genes_to_tads(conn)

    if args.map_atac:
        map_atac_to_tads(conn)

    if not args.file and not args.map_genes and not args.map_atac:
        print("Specify at least one of: --file  --map-genes  --map-atac")
        print("\nRun order:")
        print("  1. python import_tads.py --file path/to/TAD.bed")
        print("  2. python import_tads.py --map-genes")
        print("  3. python import_tads.py --map-atac")

    conn.close()

if __name__ == "__main__":
    main()
