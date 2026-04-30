#!/usr/bin/env python3
"""
Fix Gene Coordinates — Fetch Authoritative hg38 Coordinates from Ensembl
Alternative version using PyMySQL instead of mysql-connector
"""

import sys
import pymysql
import requests
import time
import argparse
from collections import defaultdict

DB_CONFIG = {
    "host": "localhost",
    "user": "daria",
    "password": "simba",
    "database": "brain_multiomics"
}

ENSEMBL_API = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"
BATCH_SIZE = 50
SLEEP_SEC = 0.05

def connect():
    try:
        return pymysql.connect(**DB_CONFIG)
    except pymysql.Error as e:
        print(f"✗ Database connection failed: {e}")
        sys.exit(1)

def fetch_ensembl_coordinates(gene_symbols):
    """Fetch hg38 coordinates from Ensembl REST API"""
    results = {}
    
    for i, symbol in enumerate(gene_symbols):
        try:
            response = requests.get(
                f"{ENSEMBL_API}/{symbol}",
                headers={"Content-Type": "application/json"},
                params={"expand": 0},
                timeout=10
            )
            
            if response.status_code == 200:
                data = response.json()
                results[symbol] = {
                    "chrom": f"chr{data.get('seq_region_name')}",
                    "start": data.get("start"),
                    "end": data.get("end"),
                    "strand": "+" if data.get("strand", 1) == 1 else "-",
                    "biotype": data.get("biotype", "unknown"),
                    "source": "Ensembl_GRCh38"
                }
            elif response.status_code == 404:
                results[symbol] = {"error": "not_found_in_ensembl"}
            else:
                results[symbol] = {"error": f"http_{response.status_code}"}
        
        except requests.exceptions.RequestException as e:
            results[symbol] = {"error": f"request_error: {str(e)[:50]}"}
        
        if (i + 1) % 10 == 0:
            print(f"  Fetched {i + 1:,}/{len(gene_symbols):,} genes...", end="\r")
        
        time.sleep(SLEEP_SEC)
    
    print(f"  Fetched {len(gene_symbols):,} genes from Ensembl API       ")
    return results

def get_current_genes():
    """Get all genes currently in gene_coordinates table"""
    conn = connect()
    cursor = conn.cursor()
    
    cursor.execute("""
        SELECT DISTINCT gene_symbol 
        FROM gene_coordinates 
        ORDER BY gene_symbol
    """)
    
    genes = [row[0] for row in cursor.fetchall()]
    cursor.close()
    conn.close()
    
    return genes

def preview_changes(genes_to_update, ensembl_data):
    """Show what will change"""
    conn = connect()
    cursor = conn.cursor()
    
    print("\n" + "="*70)
    print("PREVIEW: Sample of coordinate changes (first 10 genes)")
    print("="*70)
    
    for gene in genes_to_update[:10]:
        if "error" in ensembl_data.get(gene, {}):
            print(f"\n❌ {gene}")
            print(f"   Error: {ensembl_data[gene]['error']}")
            continue
        
        cursor.execute("""
            SELECT chrom, start_pos, end_pos, strand, genome_build
            FROM gene_coordinates
            WHERE gene_symbol = %s
            LIMIT 1
        """, (gene,))
        current = cursor.fetchone()
        
        new = ensembl_data[gene]
        
        print(f"\n✓ {gene}")
        if current:
            old_chrom, old_start, old_end, old_strand, old_build = current
            print(f"   OLD: {old_chrom}:{old_start:,}-{old_end:,} {old_strand} (build={old_build})")
        else:
            print(f"   OLD: (not in database)")
        
        print(f"   NEW: {new['chrom']}:{new['start']:,}-{new['end']:,} {new['strand']} (hg38)")
        
        if current and (old_chrom != new['chrom'] or old_start != new['start']):
            print(f"   ⚠️  COORDINATES CHANGED")
    
    cursor.close()
    conn.close()
    
    print(f"\n... and {max(0, len(genes_to_update) - 10)} more genes")
    print("="*70)

def apply_updates(genes_to_update, ensembl_data, check_only=False):
    """Update gene_coordinates with Ensembl hg38 coordinates"""
    conn = connect()
    cursor = conn.cursor()
    
    updated = 0
    not_found = 0
    errors = 0
    
    print(f"\n{'='*70}")
    if check_only:
        print("CHECK-ONLY MODE: Would update the following")
    else:
        print("Updating gene_coordinates with Ensembl hg38 coordinates...")
    print(f"{'='*70}\n")
    
    for i, gene in enumerate(genes_to_update):
        data = ensembl_data.get(gene, {})
        
        if "error" in data:
            if data["error"] == "not_found_in_ensembl":
                not_found += 1
                if i < 5:
                    print(f"  ⚠️  {gene}: Not found in Ensembl (may need manual check)")
            else:
                errors += 1
                if i < 5:
                    print(f"  ✗ {gene}: {data['error']}")
            continue
        
        if not check_only:
            try:
                cursor.execute("""
                    UPDATE gene_coordinates
                    SET chrom = %s, start_pos = %s, end_pos = %s, strand = %s,
                        genome_build = 'hg38', gene_biotype = %s
                    WHERE gene_symbol = %s
                """, (
                    data["chrom"],
                    data["start"],
                    data["end"],
                    data["strand"],
                    data.get("biotype", "unknown"),
                    gene
                ))
                updated += 1
            except pymysql.Error as e:
                print(f"  ✗ DB error updating {gene}: {e}")
                errors += 1
        else:
            updated += 1
    
    if not check_only:
        conn.commit()
    
    cursor.close()
    conn.close()
    
    print(f"\n{'='*70}")
    print(f"SUMMARY:")
    print(f"  ✓ Updated coordinates: {updated:,}")
    print(f"  ⚠️  Not found in Ensembl: {not_found:,}")
    print(f"  ✗ Errors: {errors:,}")
    print(f"{'='*70}")
    
    if not_found > 0:
        print(f"\nNote: {not_found:,} genes not found in Ensembl may be:")
        print(f"  - Retired gene names (aliases)")
        print(f"  - Non-protein-coding features")
        print(f"  - Array control probes")
        print(f"  These can be remapped manually or skipped for now.\n")
    
    return updated, not_found, errors

def get_plasticity_genes():
    """Get only the plasticity genes (most important for tracks)"""
    conn = connect()
    cursor = conn.cursor()
    
    cursor.execute("""
        SELECT DISTINCT gc.gene_symbol
        FROM gene_coordinates gc
        JOIN plasticity_genes pg ON gc.gene_symbol = pg.gene_symbol
        ORDER BY gc.gene_symbol
    """)
    
    genes = [row[0] for row in cursor.fetchall()]
    cursor.close()
    conn.close()
    
    return genes

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-only", action="store_true",
                       help="Preview changes without applying")
    parser.add_argument("--genes", type=str, default=None,
                       help="Comma-separated list of genes to update (default: plasticity genes)")
    parser.add_argument("--all-genes", action="store_true",
                       help="Update ALL genes in database (slow, ~21k genes)")
    args = parser.parse_args()
    
    print("\n" + "="*70)
    print("FIX GENE COORDINATES — Fetch hg38 from Ensembl")
    print("="*70)
    
    if args.genes:
        genes = [g.strip().upper() for g in args.genes.split(",")]
        print(f"\nUpdating {len(genes)} specified genes...")
    elif args.all_genes:
        genes = get_current_genes()
        print(f"\nFound {len(genes):,} genes in database (this will take a while)")
    else:
        genes = get_plasticity_genes()
        print(f"\nFound {len(genes):,} PLASTICITY genes to update (focused approach)")
    
    if not genes:
        print("✗ No genes to update.")
        return
    
    print(f"\nFetching hg38 coordinates from Ensembl REST API...")
    ensembl_data = fetch_ensembl_coordinates(genes)
    
    found = sum(1 for g in genes if "error" not in ensembl_data.get(g, {}))
    print(f"  Successfully fetched: {found:,}/{len(genes):,}")
    
    genes_to_update = [g for g in genes if "error" not in ensembl_data.get(g, {})]
    preview_changes(genes_to_update, ensembl_data)
    
    if args.check_only:
        print("\n✓ Preview complete (check-only mode)")
        return
    
    confirm = input("\nApply updates to gene_coordinates table? (yes/no): ").strip().lower()
    if confirm != "yes":
        print("Cancelled.")
        return
    
    apply_updates(genes_to_update, ensembl_data, check_only=False)
    print("\n✓ Done! Your gene coordinates are now aligned with Ensembl hg38.")
    print("  Methylation tracks should now display at correct genomic locations.\n")

if __name__ == "__main__":
    main()
