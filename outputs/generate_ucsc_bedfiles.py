#!/usr/bin/env python3
"""
Generate UCSC Genome Browser track files (BED format) for visualization.

Creates separate tracks for:
  - Plasticity genes (curated list)
  - ATAC-seq peaks (open chromatin)
  - H3K27ac enhancers (active enhancers)
  - HIC loop anchors (3D chromatin contacts)
  - Methylation changes (promoter hypermethylation)
  - Integration summary (genes with concordant multi-layer changes)

Usage:
  python generate_ucsc_bedfiles.py [--output-dir outputs/ucsc_tracks]

Output files:
  outputs/ucsc_tracks/
    ├── plasticity_genes.bed
    ├── atac_peaks.bed
    ├── enhancers_h3k27ac.bed
    ├── hic_loop_anchors.bed
    ├── methylation_promoters.bed
    ├── integration_summary.bed
    ├── hub.txt
    ├── genomes.txt
    └── trackDb.txt (for UCSC Hub upload)

Notes:
  - All chromosomes normalized (chr1, chrX, etc.)
  - Mitochondrial/unknown loci excluded
  - BED format: chr start end name score strand [extra fields]
  - Score normalized to 0-999 for UCSC display
"""

import os
import sys
import argparse
import mysql.connector
import pandas as pd
from datetime import datetime
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
        print(f"ERROR: Cannot connect to MySQL database")
        if err.errno == 2003:
            print(f"  - Ensure MySQL is running on {DB_CONFIG['host']}")
        elif err.errno == 1045:
            print(f"  - Check database credentials")
        sys.exit(1)

def normalize_chrom(chrom):
    """Ensure BED-format chromosome names (chr1, chrX, etc.)"""
    if not str(chrom).startswith('chr'):
        return 'chr' + str(chrom)
    return str(chrom)

def write_bed_header(f, name, description, color="100,100,100", visibility="full"):
    """Write UCSC track header"""
    f.write(f"track name=\"{name}\" description=\"{description}\" ")
    f.write(f"color={color} visibility={visibility}\n")

def export_plasticity_genes(conn, output_dir):
    """Export curated brain plasticity genes as BED"""
    cursor = conn.cursor()
    
    print("\n📍 Exporting plasticity genes...")
    cursor.execute("""
        SELECT DISTINCT 
            gc.chrom, gc.start_pos, gc.end_pos, 
            pg.gene_symbol, pg.category
        FROM gene_coordinates gc
        JOIN plasticity_genes pg ON gc.gene_symbol = pg.gene_symbol
        WHERE gc.chrom NOT IN ('MT', 'chrM', 'Un_%')
        ORDER BY gc.chrom, gc.start_pos
    """)
    
    with open(f"{output_dir}/plasticity_genes.bed", "w") as f:
        write_bed_header(f, "PlasticityGenes", 
                        "Brain plasticity genes (365 curated)", 
                        color="100,149,237", visibility="dense")
        
        count = 0
        for chrom, start, end, gene_sym, category in cursor.fetchall():
            chrom = normalize_chrom(chrom)
            # Canonical genes get higher score
            score = 800 if category == "canonical" else 500
            f.write(f"{chrom}\t{int(start)}\t{int(end)}\t{gene_sym}\t{score}\t+\n")
            count += 1
    
    cursor.close()
    print(f"  ✅ Exported {count:,} plasticity genes")

def export_atac_peaks(conn, output_dir):
    """Export ATAC-seq open chromatin peaks"""
    cursor = conn.cursor()
    
    print("\n🔓 Exporting ATAC peaks...")
    cursor.execute("""
        SELECT ap.chrom, ap.chrom_start, ap.chrom_end, ap.peak_id,
               GROUP_CONCAT(DISTINCT agl.gene_symbol SEPARATOR ';') as genes
        FROM atac_peaks ap
        LEFT JOIN atac_gene_links agl ON ap.peak_id = agl.peak_id
        WHERE ap.chrom NOT IN ('MT', 'chrM', 'Un_%')
        GROUP BY ap.chrom, ap.chrom_start, ap.chrom_end, ap.peak_id
        ORDER BY ap.chrom, ap.chrom_start
    """)
    
    with open(f"{output_dir}/atac_peaks.bed", "w") as f:
        write_bed_header(f, "ATAC_Peaks", 
                        "Open chromatin (ATAC-seq, 4603 peaks)", 
                        color="0,176,80", visibility="dense")
        
        count = 0
        for chrom, start, end, peak_id, genes in cursor.fetchall():
            chrom = normalize_chrom(chrom)
            genes_str = genes if genes else "."
            f.write(f"{chrom}\t{int(start)}\t{int(end)}\tATAC_{peak_id}\t500\t+\t\t{genes_str}\n")
            count += 1
    
    cursor.close()
    print(f"  ✅ Exported {count:,} ATAC peaks")

def export_enhancers(conn, output_dir):
    """Export H3K27ac enhancers with gene links"""
    cursor = conn.cursor()
    
    print("\n⭐ Exporting H3K27ac enhancers...")
    cursor.execute("""
        SELECT e.chrom, e.start_pos, e.end_pos, e.enhancer_id,
               GROUP_CONCAT(DISTINCT CONCAT(egl.gene_symbol,'(',egl.region_type,')') SEPARATOR ';') as genes,
               COUNT(DISTINCT egl.gene_symbol) as n_genes
        FROM enhancers e
        LEFT JOIN enhancer_gene_links egl ON e.enhancer_id = egl.enhancer_id
        WHERE e.chrom NOT IN ('MT', 'chrM', 'Un_%')
        GROUP BY e.chrom, e.start_pos, e.end_pos, e.enhancer_id
        ORDER BY e.chrom, e.start_pos
    """)
    
    with open(f"{output_dir}/enhancers_h3k27ac.bed", "w") as f:
        write_bed_header(f, "H3K27ac_Enhancers", 
                        "Active enhancers (H3K27ac ChIP-seq)", 
                        color="255,102,0", visibility="dense")
        
        count = 0
        for chrom, start, end, enh_id, genes, n_genes in cursor.fetchall():
            chrom = normalize_chrom(chrom)
            genes_str = genes if genes else "."
            # Score based on number of linked genes
            score = min(900, 400 + (n_genes or 0) * 100)
            f.write(f"{chrom}\t{int(start)}\t{int(end)}\tENH_{enh_id}\t{score}\t+\t\t{genes_str}\n")
            count += 1
    
    cursor.close()
    print(f"  ✅ Exported {count:,} enhancers")

def export_hic_loop_anchors(conn, output_dir):
    """Export HIC loop endpoints as anchors"""
    cursor = conn.cursor()
    
    print("\n🔗 Exporting HIC loop anchors...")
    cursor.execute("""
        SELECT DISTINCT chrom1, start1, end1, hic_id, contact_strength
        FROM hic_loops
        WHERE chrom1 NOT IN ('MT', 'chrM', 'Un_%')
        UNION
        SELECT DISTINCT chrom2, start2, end2, hic_id, contact_strength
        FROM hic_loops
        WHERE chrom2 NOT IN ('MT', 'chrM', 'Un_%')
        ORDER BY 1, 2
    """)
    
    rows = cursor.fetchall()
    if not rows:
        print(f"  ⚠ No HIC loops found. Run import_hic_loops.py first.")
        cursor.close()
        return
    
    with open(f"{output_dir}/hic_loop_anchors.bed", "w") as f:
        write_bed_header(f, "HIC_Anchors", 
                        "3D chromatin contact loop endpoints", 
                        color="153,51,153", visibility="dense")
        
        count = 0
        for chrom, start, end, hic_id, strength in rows:
            chrom = normalize_chrom(chrom)
            # Scale contact strength to 0-999
            score = min(999, max(100, int((strength or 1.0) * 100)))
            f.write(f"{chrom}\t{int(start)}\t{int(end)}\tHIC_{hic_id}\t{score}\t+\n")
            count += 1
    
    cursor.close()
    print(f"  ✅ Exported {count:,} HIC loop anchors")

def export_methylation_changes(conn, output_dir):
    """Export promoter methylation changes"""
    cursor = conn.cursor()
    
    print("\n🔴 Exporting methylation changes...")
    cursor.execute("""
        SELECT
            gc.chrom,
            GREATEST(gc.start_pos - 200, 0) as start,
            gc.start_pos + 200 as end,
            gc.gene_symbol,
            COUNT(*) as n_probes
        FROM gene_coordinates gc
        JOIN methylation_gene_links mgl ON gc.gene_symbol = mgl.gene_symbol
        WHERE mgl.relation IN ('TSS200', 'TSS1500')
        AND gc.chrom NOT IN ('MT', 'chrM', 'Un_%')
        GROUP BY gc.chrom, gc.start_pos, gc.end_pos, gc.gene_symbol
        ORDER BY gc.chrom, gc.start_pos
    """)
    
    with open(f"{output_dir}/methylation_promoters.bed", "w") as f:
        write_bed_header(f, "Methylation_Promoters", 
                        "Promoter methylation regions", 
                        color="200,0,0", visibility="dense")
        
        count = 0
        for chrom, start, end, gene_sym, n_probes in cursor.fetchall():
            chrom = normalize_chrom(chrom)
            score = min(900, 300 + (n_probes or 1) * 50)
            f.write(f"{chrom}\t{int(start)}\t{int(end)}\tMETH_{gene_sym}\t{score}\t+\n")
            count += 1
    
    cursor.close()
    print(f"  ✅ Exported {count:,} methylation regions")

def export_integration_summary(conn, output_dir):
    """Export genes with multi-layer evidence"""
    cursor = conn.cursor()
    
    print("\n🎯 Exporting integration summary...")
    cursor.execute("""
        SELECT
            gc.chrom,
            gc.start_pos,
            gc.end_pos,
            gc.gene_symbol,
            MAX(CASE WHEN re.expression_id IS NOT NULL THEN 1 ELSE 0 END) as has_rna,
            MAX(CASE WHEN agl.link_id IS NOT NULL THEN 1 ELSE 0 END) as has_atac,
            MAX(CASE WHEN mgl.link_id IS NOT NULL THEN 1 ELSE 0 END) as has_meth,
            MAX(CASE WHEN egl.enhancer_id IS NOT NULL THEN 1 ELSE 0 END) as has_enhancer,
            MAX(CASE WHEN hgl.hic_id IS NOT NULL THEN 1 ELSE 0 END) as has_hic,
            pg.gene_symbol IS NOT NULL as is_plasticity
        FROM gene_coordinates gc
        LEFT JOIN rna_expression re ON gc.gene_symbol = re.gene_symbol
        LEFT JOIN atac_gene_links agl ON gc.gene_symbol = agl.gene_symbol
        LEFT JOIN methylation_gene_links mgl ON gc.gene_symbol = mgl.gene_symbol
        LEFT JOIN enhancer_gene_links egl ON gc.gene_symbol = egl.gene_symbol
        LEFT JOIN hic_gene_links hgl ON gc.gene_symbol = hgl.gene_symbol
        LEFT JOIN plasticity_genes pg ON gc.gene_symbol = pg.gene_symbol
        WHERE gc.chrom NOT IN ('MT', 'chrM', 'Un_%')
        GROUP BY gc.chrom, gc.start_pos, gc.end_pos, gc.gene_symbol
        HAVING (has_rna + has_atac + has_meth) >= 2  -- Multi-layer
        ORDER BY gc.chrom, gc.start_pos
    """)
    
    with open(f"{output_dir}/integration_summary.bed", "w") as f:
        write_bed_header(f, "Multi_Layer_Integration", 
                        "Genes with evidence in 2+ data types", 
                        color="0,100,200", visibility="dense")
        
        count = 0
        for (chrom, start, end, gene_sym, has_rna, has_atac, has_meth, 
             has_enhancer, has_hic, is_plasticity) in cursor.fetchall():
            chrom = normalize_chrom(chrom)
            # Score based on number of layers
            n_layers = sum([has_rna, has_atac, has_meth, has_enhancer, has_hic])
            score = 300 + (n_layers * 100)
            if is_plasticity:
                score += 200
            score = min(999, score)
            
            f.write(f"{chrom}\t{int(start)}\t{int(end)}\t{gene_sym}_ML{n_layers}\t{score}\t+\n")
            count += 1
    
    cursor.close()
    print(f"  ✅ Exported {count:,} multi-layer integration genes")

def generate_hub_files(output_dir):
    """Generate UCSC Hub configuration files"""
    print("\n🌐 Generating UCSC Hub files...")
    
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # hub.txt
    hub_content = f"""hub Alzheimers_Epigenetics
shortLabel Alzheimers Multi-Omics
longLabel Alzheimers Disease Epigenetic Integration (AD vs Control Prefrontal Cortex)
genomesFile genomes.txt
email your_email@example.com
descriptionUrl http://your-institution.org/projects/alzheimers-omics
"""
    
    with open(f"{output_dir}/hub.txt", "w") as f:
        f.write(hub_content)
    
    # genomes.txt
    genomes_content = """genome hg19
trackDb trackDb.txt
databases hg19
"""
    
    with open(f"{output_dir}/genomes.txt", "w") as f:
        f.write(genomes_content)
    
    # trackDb.txt
    trackdb_content = """
# Plasticity Genes (core gene list)
track plasticity_genes
type bed 5
shortLabel Plasticity Genes
longLabel Brain plasticity genes curated from literature (365 genes)
visibility pack
color 100,149,237
priority 1

# ATAC-seq (open chromatin)
track atac_peaks
type bed 9
shortLabel ATAC Peaks
longLabel Open chromatin regions (ATAC-seq, 4603 peaks, GSE129040)
visibility pack
color 0,176,80
priority 2

# H3K27ac Enhancers
track enhancers_h3k27ac
type bed 9
shortLabel H3K27ac Enhancers
longLabel Active enhancers from H3K27ac ChIP-seq (GSE102538)
visibility pack
color 255,102,0
priority 3

# HIC 3D contacts
track hic_loop_anchors
type bed 6
shortLabel HIC Loop Anchors
longLabel 3D chromatin contact loop endpoints
visibility pack
color 153,51,153
priority 4

# Methylation changes
track methylation_promoters
type bed 6
shortLabel Methylation Changes
longLabel Promoter methylation regions (TSS200, TSS1500)
visibility pack
color 200,0,0
priority 5

# Integration summary
track integration_summary
type bed 6
shortLabel Multi-Layer Integration
longLabel Genes with evidence in 2+ data types (RNA + ATAC + Methylation ± Enhancers ± HIC)
visibility pack
color 0,100,200
priority 6
"""
    
    with open(f"{output_dir}/trackDb.txt", "w") as f:
        f.write(trackdb_content)
    
    print(f"  ✅ Generated UCSC Hub files")

def generate_readme(output_dir):
    """Generate README for UCSC tracks"""
    readme_content = """# UCSC Genome Browser Tracks

## Overview
This directory contains BED files for visualizing the Alzheimers Disease multi-omics integration project in the UCSC Genome Browser.

## Files

### Data Tracks
- **plasticity_genes.bed** — 365 curated brain plasticity genes
- **atac_peaks.bed** — Open chromatin regions (ATAC-seq, GSE129040)
- **enhancers_h3k27ac.bed** — Active enhancers (H3K27ac ChIP-seq, GSE102538)
- **hic_loop_anchors.bed** — 3D chromatin contact loop endpoints
- **methylation_promoters.bed** — Promoter methylation regions
- **integration_summary.bed** — Genes with multi-layer evidence

### Hub Configuration
- **hub.txt** — UCSC Hub metadata
- **genomes.txt** — Genome reference
- **trackDb.txt** — Track database definitions

## How to Use

### Option 1: Direct Upload to UCSC
1. Go to [UCSC Genome Browser](https://genome.ucsc.edu)
2. Navigate to: **Manage Custom Tracks** → **My Data** → **Upload**
3. Upload each .bed file individually
4. Set genome to **hg19** (Human Feb 2009)

### Option 2: Create UCSC Hub (requires web server)
1. Upload hub.txt, genomes.txt, trackDb.txt, and .bed files to your web server
2. Update URLs in trackDb.txt to point to your server
3. Go to: **Manage Tracks** → **Track Hubs** → **My Hubs**
4. Paste hub URL (e.g., http://your-server.com/ucsc_tracks/hub.txt)
5. Click **Add Hub**

## BED Format Reference

All files use standard BED 5-6 format:
```
chrom    start    end    name              score    strand
chr1     100000   101000 plasticity_gene   800      +
```

- **chrom**: Chromosome (chr1, chr2, chrX, etc.)
- **start**: 0-based start position
- **end**: 0-based end position (exclusive)
- **name**: Feature ID
- **score**: 0-999 (color intensity in UCSC)
- **strand**: + or -

## Visualization Tips

- **Dense view**: Shows all features without overlap issues
- **Color coding**: Each track has distinct colors for quick identification
- **Score filtering**: Use UCSC's score threshold to highlight high-confidence features
- **Comparison**: Overlay multiple tracks to find coordinated changes

## Data Integration Schema

```
RNA Expression (GSE33000)
         ↓
    ATAC Peaks (GSE129040)
         ↓
    H3K27ac Enhancers (GSE102538)
         ↓
    Methylation (GSE59685)
         ↓
    HIC 3D Contacts (TBD)
         ↓
    [ 365+ Plasticity Genes ]
```

## Questions?
Refer to STEP2_mysql_import/README.md for detailed integration methodology.

Generated: {}
""".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    
    with open(f"{output_dir}/README.md", "w") as f:
        f.write(readme_content)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate UCSC Genome Browser tracks from multi-omics database"
    )
    parser.add_argument("--output-dir", default="outputs/ucsc_tracks",
                       help="Output directory for BED files (default: outputs/ucsc_tracks)")
    
    args = parser.parse_args()
    output_dir = args.output_dir
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    print(f"📁 Output directory: {output_dir}")
    
    # Connect to database
    conn = connect()
    print(f"✅ Connected to MySQL database: {DB_CONFIG['database']}")
    
    try:
        # Export all tracks
        export_plasticity_genes(conn, output_dir)
        export_atac_peaks(conn, output_dir)
        export_enhancers(conn, output_dir)
        export_hic_loop_anchors(conn, output_dir)
        export_methylation_changes(conn, output_dir)
        export_integration_summary(conn, output_dir)
        
        # Generate hub files
        generate_hub_files(output_dir)
        generate_readme(output_dir)
        
        print(f"\n✨ All UCSC tracks generated successfully!")
        print(f"\n📊 Track files:")
        for f in sorted(os.listdir(output_dir)):
            if f.endswith('.bed'):
                size = os.path.getsize(f"{output_dir}/{f}") / 1024  # KB
                print(f"   • {f:40s} ({size:.1f} KB)")
        
        print(f"\n📋 Next steps:")
        print(f"   1. Visit https://genome.ucsc.edu/cgi-bin/hgCustom")
        print(f"   2. Load each .bed file via Custom Tracks")
        print(f"   OR upload hub.txt via Track Hubs (requires web server)")
        
    finally:
        conn.close()
