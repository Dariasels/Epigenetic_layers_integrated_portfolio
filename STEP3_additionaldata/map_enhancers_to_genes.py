"""
map_enhancers_to_genes.py
==========================
Maps H3K27ac ChIP-seq peaks (enhancers) to nearby genes,
populating the enhancer_gene_links table.

Logic is the same as ATAC peak mapping — we use gene_coordinates
to find genes within 2000bp (promoter-proximal) or up to 50kb (distal).

H3K27ac peaks near a gene's TSS = active promoter
H3K27ac peaks far from TSS = active enhancer

Run after import_chipseq_H3K27ac.py has finished.

Usage:
  python map_enhancers_to_genes.py
"""

import os
import mysql.connector

DB_CONFIG = {
    "host":     "localhost",
    "user":     "daria",
    "password": os.getenv("DB_PASSWORD", ""),
    "database": "brain_multiomics"
}

PROMOTER_WINDOW = 2000    # bp upstream of TSS
DISTAL_MAX      = 50000   # bp — anything further is not reliably linked

def connect():
    return mysql.connector.connect(**DB_CONFIG)

def map_enhancers(conn):
    cursor = conn.cursor()

    # Clear existing links
    cursor.execute("TRUNCATE TABLE enhancer_gene_links")
    conn.commit()

    # Load all enhancers
    print("  Loading enhancers...")
    cursor.execute("""
        SELECT enhancer_id, chrom,
               (start_pos + end_pos) / 2 AS midpoint
        FROM enhancers
    """)
    enhancers = cursor.fetchall()
    print(f"  {len(enhancers):,} enhancers loaded")

    # Load gene coordinates (indexed by chrom)
    print("  Loading gene coordinates...")
    cursor.execute("""
        SELECT gene_symbol, chrom, start_pos AS tss, start_pos, end_pos, strand
        FROM gene_coordinates
    """)
    from collections import defaultdict
    by_chrom = defaultdict(list)
    for row in cursor.fetchall():
        by_chrom[row[1]].append(row)
    print(f"  Genes on {len(by_chrom)} chromosomes loaded")

    sql = """
        INSERT IGNORE INTO enhancer_gene_links
            (enhancer_id, gene_symbol, distance_bp, region_type)
        VALUES (%s, %s, %s, %s)
    """

    batch    = []
    BATCH    = 5000
    inserted = 0
    unmapped = 0

    for enhancer_id, chrom, midpoint in enhancers:
        genes = by_chrom.get(chrom, [])
        # Filter to genes within DISTAL_MAX of midpoint
        nearby = [g for g in genes
                  if abs(midpoint - g[2]) <= DISTAL_MAX]

        mapped = False
        for gene_symbol, _, tss, tx_start, tx_end, strand in nearby:
            dist = abs(midpoint - tss)

            if dist <= PROMOTER_WINDOW:
                region = "promoter_proximal"
            elif tx_start <= midpoint <= tx_end:
                region = "intragenic_enhancer"
            else:
                region = "distal_enhancer"

            batch.append((enhancer_id, gene_symbol, int(dist), region))
            mapped = True

            if len(batch) >= BATCH:
                cursor.executemany(sql, batch)
                conn.commit()
                inserted += len(batch)
                batch = []
                print(f"  {inserted:,} enhancer-gene links...", end="\r")

        if not mapped:
            unmapped += 1

    if batch:
        cursor.executemany(sql, batch)
        conn.commit()
        inserted += len(batch)

    cursor.close()
    print(f"\n  ✓ {inserted:,} enhancer-gene links created")
    print(f"  ⚠ {unmapped:,} enhancers not within 50kb of any gene")

    # Check plasticity coverage
    cursor = conn.cursor()
    cursor.execute("""
        SELECT COUNT(DISTINCT egl.gene_symbol)
        FROM enhancer_gene_links egl
        JOIN plasticity_genes pg ON pg.gene_symbol = egl.gene_symbol
    """)
    n = cursor.fetchone()[0]
    cursor.close()
    print(f"  ✓ {n} plasticity genes have nearby H3K27ac peaks")

if __name__ == "__main__":
    conn = connect()
    print("✓ Connected")
    map_enhancers(conn)
    conn.close()
    print("\nDone. Run the integration queries in integration_chip_tad.sql")
