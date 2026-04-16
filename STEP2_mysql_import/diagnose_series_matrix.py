"""
diagnose_series_matrix.py
==========================
Run this FIRST on your GSE33000 file to understand its structure
without importing anything. Prints:
  - How many metadata lines are at the top
  - The exact table start/end markers
  - First 5 rows and first 5 columns of the actual data
  - Whether probe IDs look numeric or RSE_-style

Usage:
  python diagnose_series_matrix.py --file GSE33000_series_matrix.txt
  python diagnose_series_matrix.py --file GSE33000_series_matrix.txt.gz
"""

import argparse
import gzip
import io
import sys

def diagnose(filepath):
    print(f"\n{'='*60}")
    print(f"Diagnosing: {filepath}")
    print(f"{'='*60}\n")

    if filepath.endswith(".gz"):
        opener = lambda: gzip.open(filepath, "rt", encoding="utf-8", errors="replace")
    else:
        opener = lambda: open(filepath, "r", encoding="utf-8", errors="replace")

    meta_lines   = 0
    in_table     = False
    data_rows    = []
    found_start  = False
    found_end    = False

    with opener() as fh:
        for i, line in enumerate(fh):
            line_stripped = line.strip()

            if not in_table:
                # Show first 10 metadata lines
                if i < 10:
                    print(f"  Line {i+1:>4}: {line_stripped[:120]}")
                if line_stripped.startswith("!series_matrix_table_begin"):
                    found_start = True
                    in_table    = True
                    print(f"\n  >>> Table starts at line {i+1} <<<\n")
                    continue
                meta_lines += 1
            else:
                if line_stripped.startswith("!series_matrix_table_end"):
                    found_end = True
                    break
                data_rows.append(line_stripped)

    print(f"\n--- Summary ---")
    print(f"  Metadata lines before table:  {meta_lines}")
    print(f"  Table start marker found:     {found_start}")
    print(f"  Table end marker found:       {found_end}")
    print(f"  Data rows in table:           {len(data_rows)}")

    if not found_start:
        print("\n  *** WARNING: No '!series_matrix_table_begin' found! ***")
        print("  This might not be a series matrix file, OR it's a")
        print("  supplementary file (e.g. GSE33000_HBTRC_SetA_logNorm.tsv.gz).")
        print("  Check the GEO page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33000")
        print("  Look for a file ending in '_series_matrix.txt.gz'")
        return

    if data_rows:
        print(f"\n--- First 5 rows × 5 columns of the data table ---")
        header = data_rows[0].split("\t")
        print(f"  Total columns: {len(header)}")
        print(f"  First 5 column names: {header[:5]}")
        print(f"  Last  3 column names: {header[-3:]}")
        print()

        # Check probe ID format
        if len(data_rows) > 1:
            first_id = data_rows[1].split("\t")[0].strip('"')
            print(f"  First probe ID: '{first_id}'")
            if first_id.startswith("RSE_"):
                print("  ⚠  Probe IDs start with RSE_ — these are Illumina transcript IDs")
                print("     NOT standard gene symbols. Use GPL10558 annotation to map them.")
            elif first_id.isdigit() or first_id.replace(".", "").isdigit():
                print("  ✓  Probe IDs look numeric — these match GPL10558 'ID' column directly")
            else:
                print(f"  ? Unknown probe ID format: '{first_id}'")

        # Show first few data rows
        for i, row in enumerate(data_rows[1:6]):
            cols = row.split("\t")
            display = [c.strip('"') for c in cols[:5]]
            print(f"  Row {i+1}: {display}")

    print(f"\n--- What to do next ---")
    if found_start:
        print("  1. Download GPL10558 annotation:")
        print("     wget 'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz'")
        print("     gunzip GPL10558.annot.gz")
        print()
        print("  2. Run the import script:")
        print("     python 02b_import_rnaseq_UPDATED.py \\")
        print("       --file GSE33000_series_matrix.txt \\")
        print("       --annotation GPL10558.annot")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", required=True)
    args = parser.parse_args()
    diagnose(args.file)

if __name__ == "__main__":
    main()
