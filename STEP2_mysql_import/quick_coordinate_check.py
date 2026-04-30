#!/usr/bin/env python3
"""
Quick Coordinate Check — Identify suspicious coordinates without API calls
Shows genes where coordinates might be wrong (very short spans, unusual sizes, etc.)
"""
import pymysql

DB_CONFIG = {
    "host": "localhost",
    "user": "daria",
    "password": "simba",
    "database": "brain_multiomics"
}

def check_coordinates():
    """Quick sanity check on current coordinates"""
    conn = pymysql.connect(**DB_CONFIG)
    cursor = conn.cursor()
    
    print("\n" + "="*80)
    print("GENE COORDINATE SANITY CHECK")
    print("="*80)
    
    # Get plasticity genes with their current coordinates
    cursor.execute("""
        SELECT gc.gene_symbol, gc.chrom, gc.start_pos, gc.end_pos, 
               (gc.end_pos - gc.start_pos) as span_bp,
               gc.genome_build
        FROM gene_coordinates gc
        JOIN plasticity_genes pg ON gc.gene_symbol = pg.gene_symbol
        ORDER BY (gc.end_pos - gc.start_pos) ASC
        LIMIT 20
    """)
    
    print("\n🔍 Checking 20 plasticity genes with SHORTEST spans:")
    print("   (Very short spans may indicate wrong coordinates)\n")
    
    for gene, chrom, start, end, span, build in cursor.fetchall():
        if span < 1000:
            print(f"  ⚠️  {gene:15} {chrom:5} {start:12,} - {end:12,}  (span: {span:8,} bp) [build={build}]")
        else:
            print(f"  ✓ {gene:15} {chrom:5} {start:12,} - {end:12,}  (span: {span:8,} bp) [build={build}]")
    
    # Check for genes with genome_build set
    print("\n\n📊 Genome build distribution:")
    cursor.execute("""
        SELECT gc.genome_build, COUNT(*) as count
        FROM gene_coordinates gc
        JOIN plasticity_genes pg ON gc.gene_symbol = pg.gene_symbol
        GROUP BY gc.genome_build
    """)
    
    for build, count in cursor.fetchall():
        print(f"  {build or 'NULL':8} : {count:4} genes")
    
    # Sample genes for manual inspection
    print("\n\n📋 Sample of current coordinates (first 10 plasticity genes alphabetically):")
    cursor.execute("""
        SELECT gc.gene_symbol, gc.chrom, gc.start_pos, gc.end_pos, gc.genome_build
        FROM gene_coordinates gc
        JOIN plasticity_genes pg ON gc.gene_symbol = pg.gene_symbol
        ORDER BY gc.gene_symbol
        LIMIT 10
    """)
    
    print(f"\n{'Gene':<15} {'Chromosome':<12} {'Start':<14} {'End':<14} {'Build':<8}")
    print("-" * 70)
    for gene, chrom, start, end, build in cursor.fetchall():
        print(f"{gene:<15} {chrom:<12} {start:<14,} {end:<14,} {build or 'NULL':<8}")
    
    cursor.close()
    conn.close()
    
    print("\n" + "="*80)
    print("\nNOTE: If you see genes with very short spans (<100bp) or NULL genome_build,")
    print("those coordinates may be incorrect and should be manually verified/fixed.")
    print("="*80 + "\n")

if __name__ == "__main__":
    check_coordinates()
