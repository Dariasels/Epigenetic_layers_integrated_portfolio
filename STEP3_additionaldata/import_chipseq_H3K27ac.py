"""
import_chipseq_H3K27ac.py
==========================
Imports H3K27ac ChIP-seq peaks (GSE102538) into the enhancers table.

About H3K27ac ChIP-seq:
  H3K27 acetylation marks ACTIVE enhancers and promoters.
  A peak means: this region of DNA is actively being used as a
  regulatory element in this cell type under this condition.
  Unlike ATAC-seq (which just says "chromatin is open"),
  H3K27ac specifically says "this region is actively enhancing
  or promoting transcription right now."

File format: narrowPeak or broadPeak (same as ATAC)
  chr  start  end  name  score  strand  signalValue  pValue  qValue  [peak]

How to download the files:
  1. Go to: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102538
  2. Find supplementary files — look for .narrowPeak.gz or .broadPeak.gz
  3. Download all prefrontal cortex samples (NOT hippocampus)
     - Look for sample titles containing "PFC" or "dlPFC"
     - Note the GSM accession for each file

Usage:
  python import_chipseq_H3K27ac.py --file GSM3692183_peaks.narrowPeak.gz \
    --sample_id GSM3692183 \
    --cell_type NeuN+

  # Loop over all files:
  for f in data/chip/*.narrowPeak.gz; do
      gsm=$(basename $f | cut -d_ -f1)
      python import_chipseq_H3K27ac.py --file $f --sample_id $gsm --cell_type NeuN+
  done
"""

import argparse
import mysql.connector
import gzip
import sys
import os

DB_CONFIG = {
    "host":     "localhost",
    "user":     "daria",
    "password": os.getenv("DB_PASSWORD", ""),
    "database": "brain_multiomics"
}

def connect():
    return mysql.connector.connect(**DB_CONFIG)

def open_file(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def ensure_sample(sample_id, condition, cell_type, conn):
    """Add sample to samples table if not already there."""
    cursor = conn.cursor()
    cursor.execute("""
        INSERT IGNORE INTO samples (sample_id, `condition`, tissue, dataset)
        VALUES (%s, %s, %s, %s)
    """, (sample_id, condition, "prefrontal cortex", "GSE102538"))
    conn.commit()
    cursor.close()

def import_peaks(filepath, sample_id, cell_type, conn):
    cursor = conn.cursor()

    sql = """
        INSERT INTO enhancers
            (chrom, start_pos, end_pos, peak_name, signal_value,
             p_value, q_value, sample_id, cell_type, dataset)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
    """

    batch      = []
    BATCH_SIZE = 3000
    inserted   = 0

    with open_file(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue

            cols = line.split("\t")
            if len(cols) < 6:
                continue

            chrom       = cols[0].strip()
            start_pos   = int(cols[1])
            end_pos     = int(cols[2])
            peak_name   = cols[3] if len(cols) > 3 else ""
            signal_val  = float(cols[6]) if len(cols) > 6 else None
            p_val       = float(cols[7]) if len(cols) > 7 else None
            q_val       = float(cols[8]) if len(cols) > 8 else None

            # Only standard chromosomes
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            chrom_num = chrom.replace("chr", "").split("_")[0]
            if chrom_num not in [str(i) for i in range(1, 23)] + ["X", "Y", "M"]:
                continue

            batch.append((chrom, start_pos, end_pos, peak_name,
                          signal_val, p_val, q_val,
                          sample_id, cell_type, "GSE102538"))

            if len(batch) >= BATCH_SIZE:
                cursor.executemany(sql, batch)
                conn.commit()
                inserted += len(batch)
                batch = []
                print(f"  {inserted:,} peaks inserted...", end="\r")

    if batch:
        cursor.executemany(sql, batch)
        conn.commit()
        inserted += len(batch)

    cursor.close()
    print(f"\n  ✓ {inserted:,} H3K27ac peaks inserted for {sample_id} ({cell_type})")
    return inserted

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file",       required=True)
    parser.add_argument("--sample_id",  required=True)
    parser.add_argument("--cell_type",  default="bulk",
                        help="e.g. NeuN+ (neuron), Pu.1+ (microglia), OEG (oligodendrocyte)")
    parser.add_argument("--condition",  default="unknown",
                        help="Alzheimer or Control")
    args = parser.parse_args()

    print(f"\n=== Importing H3K27ac ChIP-seq: {args.sample_id} ({args.cell_type}) ===\n")
    conn = connect()
    ensure_sample(args.sample_id, args.condition, args.cell_type, conn)
    import_peaks(args.file, args.sample_id, args.cell_type, conn)
    conn.close()
    print("\n=== Done ===")
    print("Next: run map_enhancers_to_genes.py to link peaks to genes")

if __name__ == "__main__":
    main()
