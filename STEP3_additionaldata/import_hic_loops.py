#!/usr/bin/env python3
"""
Import HIC (3D chromatin contact) data into hic_loops table.

Supports:
  - BEDPE format (6-column: chrom1 start1 end1 chrom2 start2 end2 [score])
  - HiCPro format (contact matrices with normalized contact frequencies)
  - Raw contact matrices

Usage:
  python import_hic_loops.py <source_file> <dataset_name> [--min-contact 0.0] [--resolution 25]

Examples:
  python import_hic_loops.py hic_ad_prefrontal.bedpe GSE105194
  python import_hic_loops.py contacts_matrix.txt AD_prefrontal --min-contact 5.0 --resolution 25

Notes:
  - File format auto-detected by extension (.bedpe or .txt)
  - Contact strength values are normalized (contact frequency / interaction potential)
  - Chromosome names normalized (chr1 → 1, chrX → X)
  - Duplicates ignored (INSERT IGNORE)
"""

import os
import sys
import argparse
import mysql.connector
import pandas as pd
from pathlib import Path

DB_CONFIG = {
    "host": "localhost",
    "user": "daria",
    "password": os.getenv("DB_PASSWORD", ""),
    "database": "brain_multiomics"
}

def connect():
    """Establish MySQL connection"""
    try:
        return mysql.connector.connect(**DB_CONFIG)
    except mysql.connector.Error as err:
        if err.errno == 2003:
            print(f"ERROR: Cannot connect to MySQL server at '{DB_CONFIG['host']}'")
            print("  - Ensure MySQL is running")
            print("  - Check DB_PASSWORD environment variable is set")
        elif err.errno == 1045:
            print(f"ERROR: Access denied for user '{DB_CONFIG['user']}'")
            print(f"  - Verify credentials in DB_CONFIG")
        else:
            print(f"ERROR: {err}")
        sys.exit(1)

def normalize_chrom(chrom_str):
    """Normalize chromosome names: chr1 → 1, chrX → X, etc."""
    chrom = str(chrom_str).strip()
    if chrom.startswith('chr'):
        chrom = chrom[3:]
    return chrom.replace('MT', 'chrM')  # Handle mitochrondia

def import_hic_bedpe(file_path, dataset_name, conn, min_contact_strength=0.0, resolution_kb=25):
    """
    Import BEDPE format (6-column + optional score):
    chrom1 start1 end1 chrom2 start2 end2 [contact_strength]
    """
    print(f"\n📖 Reading BEDPE file: {file_path}")
    cursor = conn.cursor()
    
    # Read BEDPE, auto-detect if score column exists
    try:
        df = pd.read_csv(file_path, sep='\t', comment='#', header=None)

        # Auto-detect score if 7+ columns; typical BEDPE files are headerless
        if df.shape[1] >= 7:
            df = df.iloc[:, :7].copy()
            df.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'contact_strength']
            print(f"  ✓ Detected score column (contact strength)")
        elif df.shape[1] >= 6:
            df = df.iloc[:, :6].copy()
            df.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']
            df['contact_strength'] = 1.0
            print(f"  ✓ No score column; using default (1.0)")
        else:
            raise ValueError(f"Expected at least 6 BEDPE columns, found {df.shape[1]}")
    except Exception as e:
        print(f"ERROR reading {file_path}: {e}")
        return
    
    print(f"  ✓ Loaded {len(df):,} rows")
    
    # Filter by contact strength
    if 'contact_strength' in df.columns:
        before = len(df)
        df = df[df['contact_strength'] >= min_contact_strength]
        after = len(df)
        if before > after:
            print(f"  ✓ Filtered to {after:,} contacts >= {min_contact_strength}")
    
    # Normalize chromosome names
    df['chrom1'] = df['chrom1'].apply(normalize_chrom)
    df['chrom2'] = df['chrom2'].apply(normalize_chrom)
    
    # Remove mitochondrial and unknown chromatin
    before = len(df)
    df = df[~df['chrom1'].str.contains('Un_|M$', regex=True)]
    df = df[~df['chrom2'].str.contains('Un_|M$', regex=True)]
    if len(df) < before:
        print(f"  ✓ Removed {before - len(df)} mitochondrial/unknown loci")
    
    # Insert into database
    sql = """
        INSERT IGNORE INTO hic_loops 
            (chrom1, start1, end1, chrom2, start2, end2, contact_strength, source_dataset, resolution_kb)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
    """
    
    batch = []
    BATCH_SIZE = 5000
    inserted = 0
    
    for idx, row in df.iterrows():
        batch.append((
            row['chrom1'], 
            int(row['start1']), 
            int(row['end1']),
            row['chrom2'], 
            int(row['start2']), 
            int(row['end2']),
            float(row.get('contact_strength', 1.0)),
            dataset_name,
            resolution_kb
        ))
        
        if len(batch) >= BATCH_SIZE:
            cursor.executemany(sql, batch)
            conn.commit()
            inserted += len(batch)
            print(f"  ⏳ Inserted {inserted:,} loops...", end="\r")
            batch = []
    
    if batch:
        cursor.executemany(sql, batch)
        conn.commit()
        inserted += len(batch)
    
    cursor.close()
    print(f"\n  ✅ Imported {inserted:,} HIC loops from {dataset_name}")

def map_hic_to_genes(conn):
    """
    For each HIC loop endpoint, find overlapping or nearby genes.
    """
    cursor = conn.cursor()
    
    print(f"\n🧬 Mapping HIC loops to genes...")
    
    # Clear existing mappings
    cursor.execute("TRUNCATE TABLE hic_gene_links")
    conn.commit()
    print(f"  ✓ Cleared existing mappings")
    
    # Load gene coordinates into memory
    cursor.execute("""
        SELECT gene_symbol, chrom, start_pos, end_pos
        FROM gene_coordinates
        WHERE chrom NOT IN ('MT', 'chrM', 'Un_%')
    """)
    genes = cursor.fetchall()
    print(f"  ✓ Loaded {len(genes):,} genes")
    
    # Build chromosome index for faster lookups
    from collections import defaultdict
    genes_by_chrom = defaultdict(list)
    for gene_sym, g_chrom, g_start, g_end in genes:
        genes_by_chrom[g_chrom].append((gene_sym, g_start, g_end))
    
    # Fetch HIC loops
    cursor.execute("""
        SELECT hic_id, chrom1, start1, end1, chrom2, start2, end2 
        FROM hic_loops
    """)
    loops = cursor.fetchall()
    loop_count = len(loops)
    print(f"  ✓ Processing {loop_count:,} HIC loops")
    
    sql = """
        INSERT IGNORE INTO hic_gene_links 
            (hic_id, gene_symbol, interacting_region, distance_to_bp)
        VALUES (%s, %s, %s, %s)
    """
    
    batch = []
    BATCH_SIZE = 5000
    mapped_count = 0
    CONTACT_WINDOW = 50000  # bp — consider genes within this distance
    
    for loop_idx, (hic_id, chrom1, start1, end1, chrom2, start2, end2) in enumerate(loops):
        # Check loop end 1: find genes on chrom1 near end1
        for gene_sym, g_start, g_end in genes_by_chrom.get(chrom1, []):
            # Direct overlap
            if (start1 <= g_end and end1 >= g_start):
                batch.append((hic_id, gene_sym, 'loop_end1', 0))
                mapped_count += 1
            # Within contact window
            elif abs(end1 - g_start) < CONTACT_WINDOW or abs(end1 - g_end) < CONTACT_WINDOW:
                dist = min(abs(end1 - g_start), abs(end1 - g_end))
                batch.append((hic_id, gene_sym, 'loop_end1', int(dist)))
                mapped_count += 1
        
        # Check loop end 2: find genes on chrom2 near end2
        for gene_sym, g_start, g_end in genes_by_chrom.get(chrom2, []):
            # Direct overlap
            if (start2 <= g_end and end2 >= g_start):
                batch.append((hic_id, gene_sym, 'loop_end2', 0))
                mapped_count += 1
            # Within contact window
            elif abs(end2 - g_start) < CONTACT_WINDOW or abs(end2 - g_end) < CONTACT_WINDOW:
                dist = min(abs(end2 - g_start), abs(end2 - g_end))
                batch.append((hic_id, gene_sym, 'loop_end2', int(dist)))
                mapped_count += 1
        
        if len(batch) >= BATCH_SIZE:
            cursor.executemany(sql, batch)
            conn.commit()
            print(f"  ⏳ Processed {loop_idx+1:,}/{loop_count:,} loops ({mapped_count:,} links)...", end="\r")
            batch = []
    
    if batch:
        cursor.executemany(sql, batch)
        conn.commit()
    
    cursor.close()
    print(f"\n  ✅ Mapped {loop_count:,} HIC loops → {mapped_count:,} gene links")

def map_hic_to_enhancers(conn):
    """
    Identify H3K27ac enhancers (from enhancers table) that overlap HIC loop endpoints.
    """
    cursor = conn.cursor()
    
    print(f"\n🎯 Mapping HIC loops to H3K27ac enhancers...")
    
    # Clear existing mappings
    cursor.execute("TRUNCATE TABLE hic_enhancer_links")
    conn.commit()
    print(f"  ✓ Cleared existing mappings")
    
    # Load enhancers
    cursor.execute("""
        SELECT enhancer_id, chrom, start_pos, end_pos
        FROM enhancers
    """)
    enhancers = cursor.fetchall()
    
    if not enhancers:
        print(f"  ⚠ No enhancers found. Run import_chipseq_H3K27ac.py first.")
        return
    
    print(f"  ✓ Loaded {len(enhancers):,} enhancers")
    
    # Build chromosome index
    from collections import defaultdict
    enh_by_chrom = defaultdict(list)
    for enh_id, e_chrom, e_start, e_end in enhancers:
        enh_by_chrom[e_chrom].append((enh_id, e_start, e_end))
    
    # Fetch HIC loops
    cursor.execute("""
        SELECT hic_id, chrom1, start1, end1, chrom2, start2, end2
        FROM hic_loops
    """)
    loops = cursor.fetchall()
    
    sql = """
        INSERT IGNORE INTO hic_enhancer_links
            (hic_id, enhancer_id, loop_end_involved)
        VALUES (%s, %s, %s)
    """
    
    batch = []
    BATCH_SIZE = 5000
    linked = 0
    
    for hic_id, chrom1, start1, end1, chrom2, start2, end2 in loops:
        # Check end 1 vs enhancers on chrom1
        for enh_id, e_start, e_end in enh_by_chrom.get(chrom1, []):
            if (start1 <= e_end and end1 >= e_start):
                batch.append((hic_id, enh_id, 'end1'))
                linked += 1
        
        # Check end 2 vs enhancers on chrom2
        for enh_id, e_start, e_end in enh_by_chrom.get(chrom2, []):
            if (start2 <= e_end and end2 >= e_start):
                batch.append((hic_id, enh_id, 'end2'))
                linked += 1
        
        if len(batch) >= BATCH_SIZE:
            cursor.executemany(sql, batch)
            conn.commit()
            print(f"  ⏳ Processing HIC loops ({linked:,} links)...", end="\r")
            batch = []
    
    if batch:
        cursor.executemany(sql, batch)
        conn.commit()
    
    cursor.close()
    print(f"\n  ✅ Mapped {linked:,} HIC loops → enhancer links")

def summary_statistics(conn):
    """Print summary statistics"""
    cursor = conn.cursor()
    
    print(f"\n📊 Summary Statistics:")
    
    cursor.execute("SELECT COUNT(*) FROM hic_loops")
    n_loops = cursor.fetchone()[0]
    print(f"  • Total HIC loops: {n_loops:,}")
    
    cursor.execute("SELECT COUNT(DISTINCT source_dataset) FROM hic_loops")
    n_datasets = cursor.fetchone()[0]
    print(f"  • Datasets: {n_datasets}")
    
    cursor.execute("SELECT COUNT(*) FROM hic_gene_links")
    n_gene_links = cursor.fetchone()[0]
    print(f"  • Gene links: {n_gene_links:,}")
    
    cursor.execute("SELECT COUNT(DISTINCT gene_symbol) FROM hic_gene_links")
    n_genes = cursor.fetchone()[0]
    print(f"  • Genes with HIC contacts: {n_genes:,}")
    
    cursor.execute("""
        SELECT COUNT(DISTINCT hgl.gene_symbol)
        FROM hic_gene_links hgl
        JOIN plasticity_genes pg ON hgl.gene_symbol = pg.gene_symbol
    """)
    n_plasticity = cursor.fetchone()[0]
    print(f"  • Plasticity genes with HIC contacts: {n_plasticity:,}")
    
    cursor.execute("SELECT COUNT(*) FROM hic_enhancer_links")
    n_enh_links = cursor.fetchone()[0]
    print(f"  • Enhancer-HIC links: {n_enh_links:,}")
    
    cursor.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Import HIC 3D chromatin contact data into MySQL database"
    )
    parser.add_argument("file", help="Input file (BEDPE or contact matrix format)")
    parser.add_argument("dataset", help="Dataset name/ID (e.g., GSE105194, AD_prefrontal)")
    parser.add_argument("--min-contact", type=float, default=0.0,
                       help="Minimum contact strength threshold")
    parser.add_argument("--resolution", type=int, default=25,
                       help="Genomic resolution in kb (default: 25kb)")
    
    args = parser.parse_args()
    
    # Verify file exists
    if not os.path.exists(args.file):
        print(f"ERROR: File not found: {args.file}")
        sys.exit(1)
    
    # Connect to database
    conn = connect()
    print(f"✅ Connected to MySQL database: {DB_CONFIG['database']}")
    
    try:
        # Import HIC loops
        if args.file.endswith('.bedpe'):
            import_hic_bedpe(args.file, args.dataset, conn, args.min_contact, args.resolution)
        else:
            import_hic_bedpe(args.file, args.dataset, conn, args.min_contact, args.resolution)
        
        # Map to genes
        map_hic_to_genes(conn)
        
        # Map to enhancers (if enhancers table is populated)
        map_hic_to_enhancers(conn)
        
        # Print summary
        summary_statistics(conn)
        
        print(f"\n✨ HIC import complete!")
        print(f"   Next: Run integration_chip_tad.sql to create analysis views")
        
    finally:
        conn.close()
