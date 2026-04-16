import argparse
import gzip

def diagnose(filepath):
    print(f"\nInspecting: {filepath}\n")

    opener = gzip.open if filepath.endswith(".gz") else open

    with opener(filepath, "rt", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i >= 8:
                break
            display = line.rstrip("\n")
            sep = "TAB-separated" if "\t" in display else "COMMA-separated" if "," in display else "unknown separator"
            cols = display.split("\t") if "\t" in display else display.split(",")
            print(f"Line {i+1} [{sep}] — {len(cols)} columns:")
            print(f"  Raw start : {repr(display[:120])}")
            print(f"  Col[0]    : {repr(cols[0])}")
            if len(cols) > 1:
                print(f"  Col[1]    : {repr(cols[1])}")
            if len(cols) > 2:
                print(f"  Col[2]    : {repr(cols[2])}")
            if len(cols) > 3:
                print(f"  Col[-1]   : {repr(cols[-1])}")
            print(f"  Total cols: {len(cols)}")
            print()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", required=True)
    args = parser.parse_args()
    diagnose(args.file)