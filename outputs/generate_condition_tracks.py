#!/usr/bin/env python3
"""
Generate condition-specific UCSC tracks (AD vs Control).
Exports separate BED files for each data layer filtered by condition.
"""
import mysql.connector
import os
from pathlib import Path

DB_CONFIG = {
    "host": "localhost",
    "user": "daria",
    "password": os.getenv("DB_PASSWORD", "simba"),
    "database": "brain_multiomics"
}

OUTPUT_DIR = Path("outputs/ucsc_tracks")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def write_bed_header(f, name, description):
    f.write(f'track name="{name}" description="{description}" color=0,0,0 visibility=dense\n')

def export_atac_by_condition():
    """Export condition-specific ATAC peaks."""
    conn = mysql.connector.connect(**DB_CONFIG)
    cur = conn.cursor()
    
    for condition in ["Alzheimer", "Control"]:
        fname = OUTPUT_DIR / f"atac_peaks_{condition.lower()}.bed"
        with open(fname, "w") as f:
            write_bed_header(f, f"ATAC_{condition}", f"ATAC-seq peaks ({condition})")
            
            # Query ATAC peaks for this condition
            query = """
            SELECT DISTINCT 
                a.chrom, a.chrom_start, a.chrom_end, 
                CONCAT('ATAC_', a.peak_id), 
                LEAST(999, ROUND(a.signal_value * 10)),
                '+'
            FROM atac_peaks a
            JOIN samples s ON a.sample_id = s.sample_id
            WHERE s.`condition` = %s
            ORDER BY a.chrom, a.chrom_start
            """
            cur.execute(query, (condition,))
            for row in cur.fetchall():
                f.write("\t".join(str(x) for x in row) + "\n")
        
        print(f"Exported {fname}")
    
    conn.close()

def export_h3k27ac_by_condition():
    """Export condition-specific H3K27ac enhancers."""
    conn = mysql.connector.connect(**DB_CONFIG)
    cur = conn.cursor()
    
    for condition in ["Alzheimer", "Control"]:
        fname = OUTPUT_DIR / f"enhancers_h3k27ac_{condition.lower()}.bed"
        with open(fname, "w") as f:
            write_bed_header(f, f"Enhancers_{condition}", f"H3K27ac enhancers ({condition})")
            
            # Correct column names for enhancers table: start_pos, end_pos
            query = """
            SELECT DISTINCT 
                e.chrom, e.start_pos, e.end_pos, 
                CONCAT('ENH_', e.enhancer_id), 
                LEAST(999, ROUND(COALESCE(e.signal_value, 0) * 10)),
                '+'
            FROM enhancers e
            JOIN samples s ON e.sample_id = s.sample_id
            WHERE s.`condition` = %s
            ORDER BY e.chrom, e.start_pos
            """
            cur.execute(query, (condition,))
            for row in cur.fetchall():
                f.write("\t".join(str(x) for x in row) + "\n")
        
        print(f"✓ Exported {fname}")
    
    conn.close()

def export_methylation_by_condition():
    """Export condition-specific per-CpG methylation signal."""
    conn = mysql.connector.connect(**DB_CONFIG)
    cur = conn.cursor()
    
    for condition in ["Alzheimer", "Control"]:
        fname = OUTPUT_DIR / f"cpg_methylation_signal_{condition.lower()}.bedGraph"
        with open(fname, "w") as f:
            f.write(f'track type=bedGraph name="CpG_Methylation_{condition}" description="Mean CpG methylation ({condition})" color=200,0,0\n')
            
            query = """
            SELECT 
                p.chrom, p.start, p.end, 
                AVG(m.beta_value) as mean_beta
            FROM methylation m
            JOIN samples s ON m.sample_id = s.sample_id
            JOIN (
                SELECT DISTINCT cpg_id, chrom, position as start, position+1 as end 
                FROM cpg_annotation
            ) p ON m.cpg_id = p.cpg_id
            WHERE s.`condition` = %s
            GROUP BY m.cpg_id
            ORDER BY p.chrom, p.start
            """
            try:
                cur.execute(query, (condition,))
                for chrom, start, end, beta in cur.fetchall():
                    f.write(f"{chrom}\t{start}\t{end}\t{beta}\n")
                print(f"Exported {fname}")
            except Exception as e:
                print(f"Error exporting methylation for {condition}: {e}")
    
    conn.close()

def export_promoter_methylation_by_condition():
    """Export condition-specific promoter methylation."""
    conn = mysql.connector.connect(**DB_CONFIG)
    cur = conn.cursor()
    
    for condition in ["Alzheimer", "Control"]:
        fname = OUTPUT_DIR / f"methylation_promoters_{condition.lower()}.bed"
        with open(fname, "w") as f:
            write_bed_header(f, f"Promoter_Methylation_{condition}", f"Promoter methylation ({condition})")
            
            query = """
            SELECT DISTINCT
                mgl.chrom, mgl.start, mgl.end,
                CONCAT('METH_', mgl.gene_symbol),
                LEAST(999, ROUND(AVG(m.beta_value) * 1000)),
                '+'
            FROM methylation m
            JOIN samples s ON m.sample_id = s.sample_id
            JOIN (
                SELECT DISTINCT 
                    cpg_id, gene_symbol,
                    chrom, 
                    tss - 1500 as start,
                    tss + 1500 as end
                FROM gene_coordinates g
                WHERE tss IS NOT NULL
            ) mgl ON m.cpg_id = mgl.cpg_id
            WHERE s.`condition` = %s
            GROUP BY mgl.gene_symbol
            ORDER BY mgl.chrom, mgl.start
            """
            try:
                cur.execute(query, (condition,))
                for row in cur.fetchall():
                    f.write("\t".join(str(x) for x in row) + "\n")
                print(f"Exported {fname}")
            except Exception as e:
                print(f"Error exporting promoter methylation for {condition}: {e}")
    
    conn.close()

if __name__ == "__main__":
    print("Generating condition-specific tracks (AD vs Control)...\n")
    export_atac_by_condition()
    export_h3k27ac_by_condition()
    export_methylation_by_condition()
    export_promoter_methylation_by_condition()
    print("\n✓ All condition-specific BED files exported")
