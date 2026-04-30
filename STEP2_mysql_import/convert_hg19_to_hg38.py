#!/usr/bin/env python3
"""
Convert hg19 → hg38 Coordinates using pyliftover

This fixes the coordinate mismatch by converting all plasticity gene coordinates
from hg19 (current) to hg38 (methylation data genome build).
"""
import pymysql
from pyliftover import LiftOver

DB_CONFIG = {
    "host": "localhost",
    "user": "daria",
    "password": "simba",
    "database": "brain_multiomics"
}

def convert_coordinates():
    """Convert hg19 → hg38 coordinates for all plasticity genes"""
    
    # Initialize liftover (downloads chain file on first run)
    print("Loading hg19 → hg38 chain file...")
    lo = LiftOver('hg19', 'hg38')
    print("✓ Chain file loaded\n")
    
    conn = pymysql.connect(**DB_CONFIG)
    cursor = conn.cursor()
    
    # Fetch all plasticity genes in hg19
    print("Fetching 393 plasticity genes from database...")
    cursor.execute("""
        SELECT gc.gene_symbol, gc.chrom, gc.start_pos, gc.end_pos, gc.strand
        FROM gene_coordinates gc
        JOIN plasticity_genes pg ON gc.gene_symbol = pg.gene_symbol
        ORDER BY gc.gene_symbol
    """)
    
    genes = cursor.fetchall()
    print(f"✓ Fetched {len(genes)} genes\n")
    
    # Convert coordinates
    print("Converting coordinates from hg19 to hg38...")
    print("=" * 70)
    
    converted = 0
    failed = 0
    preview_shown = False
    
    for gene_symbol, chrom, start_pos, end_pos, strand in genes:
        # Remove 'chr' prefix if present for liftover
        chrom_clean = chrom.replace('chr', '')
        
        try:
            # Convert start position
            result_start = lo.convert_coordinate(chrom_clean, start_pos)
            if not result_start:
                print(f"  ⚠️  {gene_symbol}: Could not convert start position")
                failed += 1
                continue
            
            new_chrom_start, new_start = result_start[0]  # Take first result
            
            # Convert end position
            result_end = lo.convert_coordinate(chrom_clean, end_pos - 1)  # -1 because end is exclusive
            if not result_end:
                print(f"  ⚠️  {gene_symbol}: Could not convert end position")
                failed += 1
                continue
            
            new_chrom_end, new_end = result_end[0]
            new_end += 1  # +1 to make it exclusive again
            
            # Verify chromosomes match
            if new_chrom_start != new_chrom_end:
                print(f"  ⚠️  {gene_symbol}: Start and end on different chromosomes after liftover")
                failed += 1
                continue
            
            # Show first 5 conversions as preview
            if not preview_shown or converted < 5:
                print(f"  ✓ {gene_symbol}")
                print(f"    hg19: {chrom}:{start_pos:,}-{end_pos:,}")
                print(f"    hg38: chr{new_chrom_start}:{new_start:,}-{new_end:,}")
                if converted == 4:
                    print(f"    ... (converting {len(genes) - 5} more genes)")
                    preview_shown = True
            
            # Update database
            cursor.execute("""
                UPDATE gene_coordinates
                SET chrom = %s, start_pos = %s, end_pos = %s, genome_build = 'hg38'
                WHERE gene_symbol = %s
            """, (f"chr{new_chrom_start}", new_start, new_end, gene_symbol))
            
            converted += 1
        
        except Exception as e:
            print(f"  ✗ {gene_symbol}: {str(e)[:60]}")
            failed += 1
    
    conn.commit()
    cursor.close()
    conn.close()
    
    # Summary
    print("\n" + "=" * 70)
    print("CONVERSION SUMMARY:")
    print(f"  ✓ Successfully converted: {converted:,}/{len(genes)} genes")
    print(f"  ✗ Failed conversions:     {failed:,}")
    print("=" * 70)
    
    if converted > 0:
        print("\n✅ Coordinates successfully converted from hg19 to hg38!")
        print("   Your methylation tracks should now display at correct genomic locations.")
    
    return converted, failed

if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("CONVERT GENE COORDINATES: hg19 → hg38")
    print("=" * 70 + "\n")
    
    converted, failed = convert_coordinates()
    
    if failed == 0 and converted > 0:
        print("\n✓ All coordinates successfully converted!")
        print("  Next: Rebuild your UCSC track files to pick up new coordinates.")
    elif converted == 0:
        print("\n✗ Conversion failed - no genes were updated.")
    else:
        print(f"\n⚠️  {failed} genes had conversion failures - may need manual review.")
