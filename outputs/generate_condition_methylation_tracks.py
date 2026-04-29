#!/usr/bin/env python3
"""
Generate condition-specific methylation signal tracks (AD vs Control).
Creates bedGraph files for each condition and converts to bigWig.
"""
import mysql.connector
import os
from pathlib import Path
import subprocess

DB_CONFIG = {
    "host": "localhost",
    "user": "daria",
    "password": os.getenv("DB_PASSWORD", "simba"),
    "database": "brain_multiomics"
}

OUTPUT_DIR = Path("outputs/ucsc_tracks")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def export_methylation_to_bedgraph(condition):
    """Export condition-specific methylation to bedGraph format."""
    conn = mysql.connector.connect(**DB_CONFIG)
    cur = conn.cursor()
    
    fname = OUTPUT_DIR / f"cpg_methylation_signal_{condition.lower()}.bedGraph"
    
    with open(fname, "w") as f:
        f.write(f'track type=bedGraph name="CpG_Methylation_{condition}" description="Mean CpG methylation ({condition})" color=200,0,0\n')
        
        # Query per-CpG aggregated methylation for this condition
        # Join via cpg_annotation table to get coordinates
        query = """
        SELECT DISTINCT 
            CONCAT('chr', gc.chromosome) as chrom,
            gc.position as start,
            gc.position + 1 as end,
            AVG(m.beta_value) as mean_beta
        FROM methylation m
        JOIN samples s ON m.sample_id = s.sample_id
        JOIN (
            SELECT DISTINCT cpg_id, chromosome, position 
            FROM cpg_annotation 
            WHERE chromosome IS NOT NULL
        ) gc ON m.cpg_id = gc.cpg_id
        WHERE s.`condition` = %s
        GROUP BY m.cpg_id, gc.chromosome, gc.position
        ORDER BY gc.chromosome, gc.position
        """
        
        try:
            cur.execute(query, (condition,))
            count = 0
            for chrom, start, end, beta in cur.fetchall():
                if chrom and beta is not None:
                    # Ensure values are valid
                    beta = float(beta) if beta else 0.0
                    beta = max(0,  min(1, beta))  # clamp to [0, 1]
                    f.write(f"{chrom}\t{start}\t{end}\t{beta}\n")
                    count += 1
            print(f"✓ {fname}: {count} CpG sites")
        except Exception as e:
            print(f"⚠ Error exporting {condition} methylation: {e}")
    
    conn.close()
    return fname

def convert_bedgraph_to_bigwig(bedgraph_file):
    """Convert bedGraph to bigWig format using bedGraphToBigWig."""
    import pyBigWig
    
    bigwig_file = bedgraph_file.with_suffix('.bw')
    chrom_sizes_file = Path('tools/juicebox/tools/chrom/sizes/hg38.chrom.sizes')
    
    if not chrom_sizes_file.exists():
        print(f"⚠ chrom.sizes not found")
        return None
    
    # Read chrom sizes
    chrom_sizes = {}
    with open(chrom_sizes_file) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])
    
    try:
        bw = pyBigWig.open(str(bigwig_file), "w")
        bw.addHeader([(f"chr{i}", chrom_sizes.get(f"chr{i}", 0)) for i in list(range(1, 23)) + ['X', 'Y', 'M']])
        
        # Parse and add entries from bedGraph
        entries = {}
        with open(bedgraph_file) as f:
            for line in f:
                if line.startswith('track'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 4:
                    try:
                        chrom = parts[0]
                        start = int(parts[1])
                        end = int(parts[2])
                        value = float(parts[3])
                        if chrom not in entries:
                            entries[chrom] = []
                        entries[chrom].append((start, end, value))
                    except:
                        pass
        
        # Sort by chromosome order and write
        chrom_order = [f"chr{i}" for i in list(range(1, 23)) + ['X', 'Y', 'M']]
        for chrom in chrom_order:
            if chrom in entries:
                entries[chrom].sort()
                starts = [e[0] for e in entries[chrom]]
                ends = [e[1] for e in entries[chrom]]
                values = [e[2] for e in entries[chrom]]
                bw.addEntries(starts, ends, values=values, chromosome=chrom)
        
        bw.close()
        print(f"✓ Created bigWig: {bigwig_file}")
        return bigwig_file
    except Exception as e:
        print(f"✗ Error converting to bigWig: {e}")
        return None

if __name__ == "__main__":
    print("=" * 60)
    print("Generating condition-specific methylation tracks")
    print("=" * 60 + "\n")
    
    try:
        for condition in ["Alzheimer", "Control"]:
            print(f"\nProcessing {condition}...")
            bedgraph = export_methylation_to_bedgraph(condition)
            convert_bedgraph_to_bigwig(bedgraph)
        
        print("\n" + "=" * 60)
        print("✓ Condition-specific methylation tracks created")
        print("=" * 60)
    except Exception as e:
        print(f"\n✗ Error: {e}")
