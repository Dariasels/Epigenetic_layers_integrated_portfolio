"""
03b_map_methylation_to_genes.py
================================
Maps methylation CpG probe IDs (e.g. cg00004067) to gene symbols
and populates the `methylation_gene_links` table.

Method:
  We use Illumina's EPIC array annotation file, which tells us
  for each probe: which gene it's in, and what region (TSS200, Body, etc.)

Two ways to get the annotation:
  Option A (R):  Download via Bioconductor (most accurate, recommended)
  Option B (CSV): Download annotation CSV directly from Illumina/GEO

This script handles Option B (CSV) which doesn't require R.
For Option A, see the R code block in comments below.

Usage:
    # Download annotation first (run once):
    python 03b_map_methylation_to_genes.py --download

    # Then map:
    python 03b_map_methylation_to_genes.py --annotation EPIC_annotation.csv

--------------------------------------------------------------------
Option A — R code (run in RStudio if you prefer):

  library(BiocManager)
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  write.csv(anno[, c("Name","chr","pos","UCSC_RefGene_Name","UCSC_RefGene_Group")],
            "EPIC_annotation.csv", row.names=FALSE)

  # Note: GSE59685 used 450k array, not EPIC. If that's the case, use:
  BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  write.csv(anno[, c("Name","chr","pos","UCSC_RefGene_Name","UCSC_RefGene_Group")],
            "450k_annotation.csv", row.names=FALSE)
--------------------------------------------------------------------

Requirements:
    pip install mysql-connector-python pandas
"""

import argparse
import pandas as pd
import mysql.connector
import sys
import urllib.request
import os

DB_CONFIG = {
    "host":     "localhost",
    "user":     "root",
    "password": "your_password",   # ← CHANGE THIS
    "database": "brain_multiomics"
}

# Direct download link for Illumina 450k manifest (GSE59685 uses 450k array)
# This is a large file (~30MB). It comes from GEO platform GPL13534.
MANIFEST_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13534/soft/"
    "GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz"
)

def connect():
    conn = mysql.connector.connect(**DB_CONFIG)
    return conn

def download_annotation(outfile="450k_annotation.csv"):
    """Download and extract Illumina 450k manifest from NCBI GEO."""
    import gzip, shutil
    print(f"  Downloading 450k array annotation (this may take a minute) ...")
    tmp = outfile + ".gz"
    urllib.request.urlretrieve(MANIFEST_URL, tmp)
    with gzip.open(tmp, "rb") as f_in, open(outfile, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(tmp)
    print(f"  ✓ Saved to {outfile}")
    print(f"  Run: python 03b_map_methylation_to_genes.py --annotation {outfile}")

def load_annotation(filepath):
    """
    Load Illumina 450k annotation CSV.
    Key columns:
      - IlmnID (or Name): the cg probe ID
      - UCSC_RefGene_Name:  gene symbol(s), semicolon-separated if multiple
      - UCSC_RefGene_Group: region type (TSS200, TSS1500, Body, 1stExon, etc.)
    """
    print(f"  Loading annotation from {filepath} ...")

    # The Illumina manifest has some comment lines at the top — skip them
    df = pd.read_csv(filepath, skiprows=7, low_memory=False)

    # Rename columns to standard names if needed
    rename = {}
    cols_lower = {c.lower(): c for c in df.columns}
    if "ilmnid" in cols_lower:
        rename[cols_lower["ilmnid"]] = "probe_id"
    elif "name" in cols_lower:
        rename[cols_lower["name"]] = "probe_id"

    if "ucsc_refgene_name" in cols_lower:
        rename[cols_lower["ucsc_refgene_name"]] = "gene_name"
    if "ucsc_refgene_group" in cols_lower:
        rename[cols_lower["ucsc_refgene_group"]] = "gene_group"

    df.rename(columns=rename, inplace=True)
    df = df[["probe_id", "gene_name", "gene_group"]].copy()
    df = df[df["probe_id"].str.startswith("cg", na=False)]  # keep only cg probes

    print(f"  ✓ Loaded {len(df):,} probe annotations")
    return df

def get_probe_ids_in_db(conn):
    """Return the set of probe IDs that are actually in our methylation table."""
    print("  Fetching probe IDs from database ...")
    cursor = conn.cursor()
    cursor.execute("SELECT DISTINCT cpg_id FROM methylation")
    ids = {row[0] for row in cursor.fetchall()}
    cursor.close()
    print(f"  ✓ Found {len(ids):,} unique probes in methylation table")
    return ids

def insert_links(annotation_df, probe_ids_in_db, conn):
    """
    For each probe in our data, extract gene links and insert into
    methylation_gene_links.

    A probe can map to multiple genes (semicolon-separated in annotation).
    """
    cursor = conn.cursor()

    sql = """
        INSERT IGNORE INTO methylation_gene_links (cpg_id, gene_symbol, relation)
        VALUES (%s, %s, %s)
    """

    batch      = []
    batch_size = 5000
    total      = 0
    unmapped   = 0

    for _, row in annotation_df.iterrows():
        probe_id = str(row["probe_id"]).strip()

        # Only process probes we actually have data for
        if probe_id not in probe_ids_in_db:
            continue

        gene_str  = str(row.get("gene_name",  "")).strip()
        group_str = str(row.get("gene_group", "")).strip()

        if not gene_str or gene_str in ("nan", ""):
            unmapped += 1
            continue

        # Genes can be semicolon-separated (probe overlaps multiple genes)
        genes  = gene_str.split(";")
        groups = group_str.split(";") if group_str and group_str != "nan" else []

        for i, gene in enumerate(genes):
            gene = gene.strip()
            if not gene:
                continue
            relation = groups[i].strip() if i < len(groups) else "unknown"
            batch.append((probe_id, gene, relation))

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
    print(f"\n  ✓ methylation_gene_links: {total:,} links inserted")
    print(f"  ⚠ {unmapped} probes had no gene annotation")

def main():
    parser = argparse.ArgumentParser(description="Map CpG probes to genes")
    parser.add_argument("--annotation", help="Path to Illumina 450k annotation CSV")
    parser.add_argument("--download",   action="store_true",
                        help="Download annotation from NCBI and exit")
    args = parser.parse_args()

    if args.download:
        download_annotation()
        return

    if not args.annotation:
        print("ERROR: provide --annotation file or use --download to get one")
        sys.exit(1)

    print("\n=== Mapping methylation probes to genes ===")
    conn           = connect()
    probe_ids      = get_probe_ids_in_db(conn)
    annotation     = load_annotation(args.annotation)
    insert_links(annotation, probe_ids, conn)
    conn.close()
    print("=== Done ===\n")

if __name__ == "__main__":
    main()
