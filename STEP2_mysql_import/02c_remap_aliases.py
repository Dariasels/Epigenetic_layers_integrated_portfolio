"""
remap_aliases.py
=================
Retries the 14,032 unmapped genes from unmapped_genes.txt using
MyGene.info's 'alias' scope, which covers retired/old gene names.

Why this works:
  Names like 'A2BP1', '182-FIP', 'AA019935' are old aliases that have
  since been renamed (e.g. A2BP1 → RBFOX1). MyGene.info knows about
  these under scopes='alias,retired,old_names' but not under 'symbol'.

Run AFTER map_ATACcoordinates_genes3.py has finished.

Usage:
  python remap_aliases.py --genome hg19
  python remap_aliases.py --genome hg19 --input unmapped_genes.txt
"""

import os
import mysql.connector
import requests
import time
import argparse

DB_CONFIG = {
    "host":     "localhost",
    "user":     "daria",
    "password": os.getenv("DB_PASSWORD", ""),
    "database": "brain_multiomics"
}

MYGENE_URL  = "https://mygene.info/v3/query"
BATCH_SIZE  = 200   # smaller batches — alias search is fuzzier
SLEEP_SEC   = 0.3

GENOME_FIELD = {
    "hg19": "genomic_pos_hg19",
    "hg38": "genomic_pos",
}

def connect():
    return mysql.connector.connect(**DB_CONFIG)

def parse_item(item, genome_field):
    if not isinstance(item, dict) or item.get("notfound"):
        return None
    symbol = item.get("symbol") or item.get("query")
    if not symbol:
        return None
    pos = item.get(genome_field)
    if not pos:
        return None
    if isinstance(pos, list):
        standard = [p for p in pos
                    if str(p.get("chr","")).replace("chr","").split(".")[0]
                    in [str(i) for i in range(1,23)] + ["X","Y","M"]]
        pos = standard[0] if standard else pos[0]
    if not isinstance(pos, dict) or "chr" not in pos or "start" not in pos:
        return None
    chrom = str(pos["chr"])
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    try:
        start  = int(pos["start"])
        end    = int(pos.get("end", pos["start"]))
        strand = "+" if pos.get("strand") == 1 else "-"
    except (ValueError, TypeError):
        return None
    return (symbol, chrom, start, end, strand)

def query_alias_batch(genes, genome_field):
    """Query using alias+symbol scopes together."""
    params = {
        "q":       ",".join(genes),
        "scopes":  "symbol,alias",   # alias catches old/retired names
        "fields":  f"{genome_field},symbol,alias",
        "species": "human",
        "size":    len(genes),
    }
    try:
        r = requests.post(MYGENE_URL, data=params, timeout=30)
        if r.status_code == 200:
            data = r.json()
            return data if isinstance(data, list) else []
        elif r.status_code == 400:
            # Split and retry
            if len(genes) <= 1:
                return []
            mid = len(genes) // 2
            return (query_alias_batch(genes[:mid], genome_field) +
                    query_alias_batch(genes[mid:], genome_field))
    except requests.exceptions.RequestException:
        pass
    return []

def run(genome, input_file):
    genome_field = GENOME_FIELD[genome]

    # Load unmapped genes
    try:
        with open(input_file) as f:
            genes = [l.strip() for l in f if l.strip()]
    except FileNotFoundError:
        print(f"✗ Could not find {input_file}")
        print("  Run map_ATACcoordinates_genes3.py first.")
        return

    print(f"Loaded {len(genes):,} unmapped genes from {input_file}")
    print(f"Retrying with alias scope (hg19)...\n")

    conn   = connect()
    cursor = conn.cursor()

    # Only retry genes still not in gene_coordinates
    cursor.execute("SELECT gene_symbol FROM gene_coordinates")
    already_mapped = {row[0] for row in cursor.fetchall()}
    genes = [g for g in genes if g not in already_mapped]
    print(f"After excluding already-mapped: {len(genes):,} to retry\n")

    inserted    = 0
    still_lost  = 0
    alias_wins  = []   # track which old names got resolved

    for i in range(0, len(genes), BATCH_SIZE):
        batch   = genes[i : i + BATCH_SIZE]
        results = query_alias_batch(batch, genome_field)

        for item in results:
            if item.get("notfound"):
                still_lost += 1
                continue
            parsed = parse_item(item, genome_field)
            if parsed is None:
                still_lost += 1
                continue

            symbol, chrom, start, end, strand = parsed
            original_query = item.get("query", symbol)

            try:
                # Insert using the CURRENT canonical symbol (not the old alias)
                cursor.execute("""
                    INSERT IGNORE INTO gene_coordinates
                        (gene_symbol, chrom, start_pos, end_pos, strand, genome_build)
                    VALUES (%s, %s, %s, %s, %s, %s)
                """, (symbol, chrom, start, end, strand, genome))

                # Also insert a row for the OLD alias name so rna_expression
                # rows that use the old name can still be joined
                if symbol != original_query:
                    cursor.execute("""
                        INSERT IGNORE INTO gene_coordinates
                            (gene_symbol, chrom, start_pos, end_pos, strand, genome_build)
                        VALUES (%s, %s, %s, %s, %s, %s)
                    """, (original_query, chrom, start, end, strand, genome))
                    alias_wins.append(f"{original_query} → {symbol}")

                inserted += 1
            except mysql.connector.Error as e:
                print(f"  ✗ DB error for {symbol}: {e}")

        conn.commit()
        pct = min(i + len(batch), len(genes))
        print(f"  [{pct:>6,} / {len(genes):,}]  "
              f"rescued={inserted:,}  still_lost={still_lost:,}",
              end="\r")
        time.sleep(SLEEP_SEC)

    cursor.close()
    conn.close()

    print(f"\n\n{'='*55}")
    print(f"✓ DONE")
    print(f"  Genes rescued via alias     : {inserted:,}")
    print(f"  Still not found             : {still_lost:,}")
    print(f"{'='*55}")

    if alias_wins:
        print(f"\nSample alias resolutions (old name → current symbol):")
        for line in alias_wins[:15]:
            print(f"  {line}")
        if len(alias_wins) > 15:
            print(f"  ... and {len(alias_wins)-15} more")

    if still_lost > 0:
        print(f"\nThe {still_lost:,} remaining unfound genes are likely:")
        print("  - Probesets with no corresponding human gene")
        print("  - Array control probes")
        print("  - Genes with no hg19 coordinates in any database")
        print("  These won't link to ATAC peaks — that's expected and acceptable.")
        print(f"  For context: your array has ~47,000 probes but only ~21,000")
        print(f"  are protein-coding genes. The rest are non-coding, controls, etc.")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", choices=["hg19","hg38"], default="hg19")
    parser.add_argument("--input",  default="unmapped_genes.txt",
                        help="File of unmapped gene names (default: unmapped_genes.txt)")
    args = parser.parse_args()
    run(args.genome, args.input)

if __name__ == "__main__":
    main()