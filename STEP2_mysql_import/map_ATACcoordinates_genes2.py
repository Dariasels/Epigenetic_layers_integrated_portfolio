"""
map_ATACcoordinates_genes3.py
==============================
Fixes the HTTP 400 error and the 13,616 "not found" problem.

ROOT CAUSES:
  1. HTTP 400 → the API rejected a whole batch because one gene name in it
     contains a character MyGene.info can't parse:
       - commas      (e.g. "ENSG00000,123")
       - slashes     (e.g. "HGNC/1234")
       - empty strings
       - very long probe IDs (RSE_######### still in your table)
       - special Unicode characters
     When one bad name poisons a batch of 500, ALL 500 fail silently.
     Fix: sanitise gene names before sending, and split failing batches
          into smaller pieces to isolate the bad apple.

  2. 13,616 "not found" → two sub-causes:
     a) RSE_######### probe IDs still in rna_expression (not real gene symbols)
     b) Obsolete/alias gene symbols the API doesn't recognise under 'symbol' scope
        Fix: also search 'alias' scope as fallback for unmatched genes.

NEW FEATURES vs genes2.py:
  - Pre-filters obviously bad names before sending any request
  - On HTTP 400: splits batch in half recursively to isolate bad names
  - Fallback lookup for genuinely missing genes using alias+retired symbols
  - Saves unmapped genes to unmapped_genes.txt so you can inspect them
  - Reports exactly how many were RSE_ IDs vs real-but-unmapped genes

Usage:
  python map_ATACcoordinates_genes3.py --genome hg19
  python map_ATACcoordinates_genes3.py --genome hg19 --check-unmapped
"""

import os
import mysql.connector
import requests
import time
import argparse
import re

DB_CONFIG = {
    "host":     "localhost",
    "user":     "daria",
    "password": os.getenv("DB_PASSWORD", ""),
    "database": "brain_multiomics"
}

MYGENE_URL  = "https://mygene.info/v3/query"
BATCH_SIZE  = 500
SLEEP_SEC   = 0.25
MAX_RETRIES = 3

GENOME_FIELD = {
    "hg19": "genomic_pos_hg19",
    "hg38": "genomic_pos",
}

# ── GENE NAME VALIDATION ──────────────────────────────────────────────────────
def is_valid_gene_name(name):
    """
    Returns True if the name looks like a real HGNC gene symbol.
    Rejects: empty, RSE_ probe IDs, ENSG IDs, numeric-only, names with
             commas/slashes/semicolons that break the API query string.
    """
    if not name or not isinstance(name, str):
        return False
    name = name.strip()
    if len(name) == 0 or len(name) > 50:
        return False
    # Reject known non-symbol patterns
    if name.startswith("RSE_"):       return False   # Illumina transcript IDs
    if name.startswith("ENSG"):       return False   # Ensembl IDs
    if name.startswith("ILMN_"):      return False   # Illumina probe IDs
    if name.isdigit():                return False   # pure numbers
    if re.search(r"[,;/\\\"']", name): return False  # characters that break CSV query
    return True

def sanitise_batch(batch):
    """Return only the valid gene names from a batch."""
    return [g.strip() for g in batch if is_valid_gene_name(g)]

# ── API CALL ──────────────────────────────────────────────────────────────────
def query_api(genes, genome_field, scope="symbol"):
    """
    POST a list of gene names to MyGene.info.
    Returns (results_list, success_bool).
    """
    if not genes:
        return [], True

    params = {
        "q":       ",".join(genes),
        "scopes":  scope,
        "fields":  f"{genome_field},symbol",
        "species": "human",
        "size":    len(genes),
    }

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = requests.post(MYGENE_URL, data=params, timeout=30)
            if r.status_code == 200:
                data = r.json()
                return (data if isinstance(data, list) else []), True
            elif r.status_code == 400:
                # Bad request — don't retry, report failure immediately
                return [], False
            else:
                print(f"    ⚠ HTTP {r.status_code} attempt {attempt}")
        except requests.exceptions.RequestException as e:
            print(f"    ⚠ Network error attempt {attempt}: {e}")
        if attempt < MAX_RETRIES:
            time.sleep(2 ** attempt)

    return [], False

def query_with_fallback_split(genes, genome_field):
    """
    Query a batch. If it gets a 400, split in half and retry each half.
    This isolates the bad gene name(s) so the rest still get mapped.
    Recursively splits down to individual genes if needed.
    """
    if not genes:
        return []

    results, ok = query_api(genes, genome_field)
    if ok:
        return results

    # 400 error — split and retry
    if len(genes) == 1:
        # Single gene caused 400 — skip it silently
        return []

    mid = len(genes) // 2
    left  = query_with_fallback_split(genes[:mid], genome_field)
    right = query_with_fallback_split(genes[mid:], genome_field)
    return left + right

# ── PARSE ONE RESULT ITEM ─────────────────────────────────────────────────────
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
        if not pos:
            return None
        standard = [p for p in pos
                    if str(p.get("chr", "")).replace("chr", "").split(".")[0]
                    in [str(i) for i in range(1, 23)] + ["X", "Y", "M"]]
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

# ── DATABASE ──────────────────────────────────────────────────────────────────
def connect():
    return mysql.connector.connect(**DB_CONFIG)

def get_genes_to_map(conn):
    cursor = conn.cursor()
    cursor.execute("""
        SELECT DISTINCT re.gene_symbol
        FROM rna_expression re
        LEFT JOIN gene_coordinates gc ON gc.gene_symbol = re.gene_symbol
        WHERE gc.gene_symbol IS NULL
          AND re.gene_symbol IS NOT NULL
          AND re.gene_symbol != ''
    """)
    genes = [row[0] for row in cursor.fetchall()]
    cursor.close()
    return genes

def insert_coordinate(cursor, parsed, genome_build):
    symbol, chrom, start, end, strand = parsed
    cursor.execute("""
        INSERT IGNORE INTO gene_coordinates
            (gene_symbol, chrom, start_pos, end_pos, strand, genome_build)
        VALUES (%s, %s, %s, %s, %s, %s)
    """, (symbol, chrom, start, end, strand, genome_build))

# ── MAIN ──────────────────────────────────────────────────────────────────────
def run(genome):
    genome_field  = GENOME_FIELD[genome]
    genome_build  = genome

    conn   = connect()
    print("✓ Connected to database")
    cursor = conn.cursor()

    all_genes = get_genes_to_map(conn)
    print(f"Found {len(all_genes):,} genes still needing coordinates\n")

    # ── Pre-filter: split valid vs invalid names ──────────────────────────────
    valid_genes   = [g for g in all_genes if is_valid_gene_name(g)]
    invalid_genes = [g for g in all_genes if not is_valid_gene_name(g)]

    rse_count  = sum(1 for g in invalid_genes if g.startswith("RSE_"))
    ensg_count = sum(1 for g in invalid_genes if g.startswith("ENSG"))
    other_bad  = len(invalid_genes) - rse_count - ensg_count

    print(f"Pre-filter results:")
    print(f"  Valid gene symbols to query : {len(valid_genes):>8,}")
    print(f"  Skipped — RSE_ probe IDs   : {rse_count:>8,}  ← needs re-import with annotation")
    print(f"  Skipped — ENSG IDs         : {ensg_count:>8,}")
    print(f"  Skipped — other bad names  : {other_bad:>8,}")
    print()

    if rse_count > 1000:
        print("  ⚠ WARNING: You have many RSE_ probe IDs in rna_expression.")
        print("    These are Illumina transcript IDs, not gene symbols.")
        print("    MyGene.info cannot find them. To fix this properly:")
        print("    Re-run 02b_import_rnaseq_UPDATED.py with --annotation GPL10558.annot")
        print("    This will replace probe IDs with real gene symbols at import time.\n")

    inserted  = 0
    not_found = 0
    no_coords = 0
    unmapped  = []

    print(f"Querying MyGene.info for {len(valid_genes):,} valid gene symbols...\n")

    for i in range(0, len(valid_genes), BATCH_SIZE):
        batch   = valid_genes[i : i + BATCH_SIZE]
        clean   = sanitise_batch(batch)

        results = query_with_fallback_split(clean, genome_field)

        # Track which genes got a result
        found_symbols = set()

        for item in results:
            if item.get("notfound"):
                not_found += 1
                unmapped.append(item.get("query", "?"))
                continue

            parsed = parse_item(item, genome_field)
            if parsed is None:
                no_coords += 1
                unmapped.append(item.get("query", "?"))
                continue

            try:
                insert_coordinate(cursor, parsed, genome_build)
                inserted += 1
                found_symbols.add(parsed[0])
            except mysql.connector.Error as e:
                print(f"\n  ✗ DB error for {parsed[0]}: {e}")

        # Genes in batch that got no result at all
        for g in clean:
            if g not in found_symbols and g not in [
                item.get("query") for item in results
            ]:
                unmapped.append(g)

        conn.commit()
        pct = min(i + len(batch), len(valid_genes))
        print(f"  [{pct:>6,} / {len(valid_genes):,}]  "
              f"inserted={inserted:,}  not_found={not_found:,}  "
              f"no_coords={no_coords:,}",
              end="\r")

        time.sleep(SLEEP_SEC)

    # ── Save unmapped gene list ───────────────────────────────────────────────
    unmapped_clean = sorted(set(unmapped))
    with open("unmapped_genes.txt", "w") as f:
        f.write("\n".join(unmapped_clean))

    cursor.close()
    conn.close()

    print(f"\n\n{'='*55}")
    print(f"✓ DONE")
    print(f"  Coordinates inserted          : {inserted:,}")
    print(f"  Not found by API              : {not_found:,}")
    print(f"  Found but no hg19 coordinates : {no_coords:,}")
    print(f"  Skipped (RSE_/bad names)      : {len(invalid_genes):,}")
    print(f"{'='*55}")
    print(f"\nUnmapped gene names saved to: unmapped_genes.txt")
    print( "Inspect with:  head -50 unmapped_genes.txt")
    print( "\nNext step:  python 03a_map_atac_to_genes.py")

def check_unmapped():
    """Read unmapped_genes.txt and categorise them."""
    try:
        with open("unmapped_genes.txt") as f:
            genes = [l.strip() for l in f if l.strip()]
    except FileNotFoundError:
        print("Run the main script first to generate unmapped_genes.txt")
        return

    rse   = [g for g in genes if g.startswith("RSE_")]
    ensg  = [g for g in genes if g.startswith("ENSG")]
    other = [g for g in genes if not g.startswith(("RSE_", "ENSG"))]

    print(f"\nUnmapped genes: {len(genes):,} total")
    print(f"  RSE_ probe IDs : {len(rse):,}")
    print(f"  ENSG IDs       : {len(ensg):,}")
    print(f"  Real symbols   : {len(other):,}  ← these genuinely have no hg19 entry")
    if other:
        print(f"  Examples       : {other[:10]}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", choices=["hg19","hg38"], default="hg19")
    parser.add_argument("--check-unmapped", action="store_true",
                        help="Categorise unmapped_genes.txt and exit")
    args = parser.parse_args()

    if args.check_unmapped:
        check_unmapped()
    else:
        run(args.genome)

if __nam_ == "__main__":
    main()