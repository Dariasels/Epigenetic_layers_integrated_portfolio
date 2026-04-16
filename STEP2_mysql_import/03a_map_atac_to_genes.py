"""
03a_map_atac_to_genes.py
=========================
Maps ATAC-seq peaks to nearby genes and populates the `atac_gene_links` table.

Strategy:
  For each ATAC peak (chr, start, end), we look for genes whose
  Transcription Start Site (TSS) falls within a defined window.

  Window definitions:
    - Promoter:    peak overlaps within 2000 bp upstream of TSS
    - Distal:      peak is 2000–50000 bp from TSS
    - Intragenic:  peak falls inside the gene body (TSS → TES)

Gene annotation source:
  We use a pre-downloaded TSV file from GENCODE or UCSC.
  Download instructions are printed when you run this script with --download.

Usage:
    # First, download gene annotation:
    python 03a_map_atac_to_genes.py --download

    # Then run the mapping:
    python 03a_map_atac_to_genes.py --annotation gene_annotation.tsv

Requirements:
    pip install mysql-connector-python pandas
"""

import argparse
import pandas as pd
import mysql.connector
import sys
import os
import urllib.request

DB_CONFIG = {
    "host":     "localhost",
    "user":     "daria",
    "password": "simba",   # ← CHANGE THIS
    "database": "brain_multiomics"
}

# Window for "promoter" region (bp upstream of TSS)
PROMOTER_UPSTREAM   = 2000
PROMOTER_DOWNSTREAM = 500
# Maximum distance for "distal" peaks
DISTAL_MAX = 50000

# UCSC GENCODE hg19 gene annotation (basic set, TSV format)
# Using hg19 because GSE59685 methylation data is often hg19
ANNOTATION_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz"
)
ANNOTATION_COLS = [
    "gene_symbol", "transcript_id", "chrom", "strand",
    "txStart", "txEnd", "cdsStart", "cdsEnd",
    "exonCount", "exonStarts", "exonEnds"
]

def connect():
    conn = mysql.connector.connect(**DB_CONFIG)
    return conn

def download_annotation(outfile="gene_annotation.tsv"):
    """Download refFlat gene annotation from UCSC."""
    print(f"  Downloading gene annotation from UCSC ...")
    import gzip, shutil
    tmp = outfile + ".gz"
    urllib.request.urlretrieve(ANNOTATION_URL, tmp)
    with gzip.open(tmp, "rb") as f_in, open(outfile, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(tmp)
    print(f"  ✓ Saved to {outfile}")
    print(f"  Run: python 03a_map_atac_to_genes.py --annotation {outfile}")

def load_annotation(filepath):
    """Load gene annotation file into DataFrame, keep one entry per gene+chrom."""
    print(f"  Loading gene annotation from {filepath} ...")
    df = pd.read_csv(filepath, sep="\t", header=None, names=ANNOTATION_COLS)

    # Keep relevant columns only
    df = df[["gene_symbol", "chrom", "strand", "txStart", "txEnd"]].copy()
    df = df[df["chrom"].str.match(r"^chr[\dXYM]+$")]  # only standard chroms
    df = df.drop_duplicates(subset=["gene_symbol", "chrom"])

    # TSS = txStart for + strand genes, txEnd for - strand genes
    df["tss"] = df.apply(
        lambda r: r["txStart"] if r["strand"] == "+" else r["txEnd"], axis=1
    )

    print(f"  ✓ Loaded {len(df):,} gene records")
    return df

def load_peaks_from_db(conn):
    """Load all ATAC peaks from the database."""
    print("  Loading peaks from database ...")
    cursor = conn.cursor(dictionary=True)
    cursor.execute("SELECT peak_id, chrom, chrom_start, chrom_end FROM atac_peaks")
    peaks = cursor.fetchall()
    cursor.close()
    print(f"  ✓ Loaded {len(peaks):,} peaks")
    return peaks

def classify_region(peak_mid, tss, tx_start, tx_end, strand):
    """Return region type based on peak position relative to gene."""
    dist = abs(peak_mid - tss)

    # Check if peak is within promoter window
    if strand == "+":
        upstream   = tss - PROMOTER_UPSTREAM
        downstream = tss + PROMOTER_DOWNSTREAM
    else:
        upstream   = tss - PROMOTER_DOWNSTREAM
        downstream = tss + PROMOTER_UPSTREAM

    if upstream <= peak_mid <= downstream:
        return "promoter", dist

    # Check if peak is inside gene body
    if tx_start <= peak_mid <= tx_end:
        return "intragenic", dist

    # Check distal
    if dist <= DISTAL_MAX:
        return "distal", dist

    return None, dist

def map_peaks_to_genes(peaks, annotation_df, conn):
    """
    For each peak, find overlapping or nearby genes and insert into atac_gene_links.
    """
    cursor = conn.cursor()

    # Index annotation by chrom for fast lookup
    by_chrom = {chrom: grp for chrom, grp in annotation_df.groupby("chrom")}

    sql = """
        INSERT IGNORE INTO atac_gene_links (peak_id, gene_symbol, distance_bp, region_type)
        VALUES (%s, %s, %s, %s)
    """

    batch        = []
    batch_size   = 3000
    total        = 0
    unmapped     = 0

    for peak in peaks:
        peak_id    = peak["peak_id"]
        chrom      = peak["chrom"]
        peak_mid   = (peak["chrom_start"] + peak["chrom_end"]) // 2

        genes_on_chrom = by_chrom.get(chrom)
        if genes_on_chrom is None:
            unmapped += 1
            continue

        # Filter to genes within DISTAL_MAX of peak midpoint (for speed)
        nearby = genes_on_chrom[
            (genes_on_chrom["tss"] >= peak_mid - DISTAL_MAX) &
            (genes_on_chrom["tss"] <= peak_mid + DISTAL_MAX)
        ]

        mapped_any = False
        for _, gene in nearby.iterrows():
            region_type, dist = classify_region(
                peak_mid, gene["tss"], gene["txStart"], gene["txEnd"], gene["strand"]
            )
            if region_type is not None:
                batch.append((peak_id, gene["gene_symbol"], int(dist), region_type))
                mapped_any = True

        if not mapped_any:
            unmapped += 1

        if len(batch) >= batch_size:
            cursor.executemany(sql, batch)
            conn.commit()
            total += len(batch)
            batch = []
            print(f"    ... {total:,} links inserted", end="\r")

    if batch:
        cursor.executemany(sql, batch)
        conn.commit()
        total += len(batch)

    cursor.close()
    print(f"\n  ✓ atac_gene_links: {total:,} links inserted")
    print(f"  ⚠ {unmapped:,} peaks could not be mapped to any nearby gene")

def main():
    parser = argparse.ArgumentParser(description="Map ATAC peaks to genes")
    parser.add_argument("--annotation", help="Path to gene annotation TSV (refFlat)")
    parser.add_argument("--download",   action="store_true",
                        help="Download UCSC refFlat annotation and exit")
    args = parser.parse_args()

    if args.download:
        download_annotation()
        return

    if not args.annotation:
        print("ERROR: provide --annotation file or use --download to get one")
        sys.exit(1)

    print("\n=== Mapping ATAC peaks to genes ===")
    annotation = load_annotation(args.annotation)
    conn       = connect()
    peaks      = load_peaks_from_db(conn)
    map_peaks_to_genes(peaks, annotation, conn)
    conn.close()
    print("=== Done ===\n")

if __name__ == "__main__":
    main()
