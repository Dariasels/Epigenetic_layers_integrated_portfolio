#!/usr/bin/env python3
"""
Lift the GSE102538 H3K27ac enhancer BED tracks from hg19 to hg38.

Inputs:
  - outputs/ucsc_tracks/enhancers_h3k27ac.bed
  - outputs/ucsc_tracks/enhancers_h3k27ac_ad.bed
  - outputs/ucsc_tracks/enhancers_h3k27ac_control.bed

Outputs:
  - outputs/ucsc_tracks/h3k27ac_hg38/*.bed
  - outputs/ucsc_tracks/h3k27ac_hg38/*.bb
  - outputs/ucsc_tracks/hub_h3k27ac_hg38.txt
  - outputs/ucsc_tracks/genomes_h3k27ac_hg38.txt
  - outputs/ucsc_tracks/trackDb_h3k27ac_hg38.txt

Requires:
  - pyliftover
  - bedToBigBed
  - hg38.chrom.sizes
"""
from __future__ import annotations

from pathlib import Path
import urllib.request
import subprocess

from pyliftover import LiftOver

ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT / "outputs" / "ucsc_tracks"
OUT_DIR = SRC_DIR / "h3k27ac_hg38"
CHROMSIZES = ROOT / "tools" / "juicebox" / "tools" / "chrom" / "sizes" / "hg38.chrom.sizes"
BEDTOBIGBED = Path("/home/daria/anaconda3/bin/bedToBigBed")
CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
CHAIN_FILE = Path("/tmp/hg19ToHg38.over.chain.gz")
STANDARD_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}

TRACKS = [
    ("enhancers_h3k27ac.bed", "enhancers_h3k27ac", "H3K27ac Enhancers - All", "H3K27ac enhancer track (all samples) lifted to hg38", "255,102,0"),
    ("enhancers_h3k27ac_ad.bed", "enhancers_h3k27ac_ad", "H3K27ac Enhancers - AD", "H3K27ac enhancer signal in Alzheimer's samples lifted to hg38", "204,0,0"),
    ("enhancers_h3k27ac_control.bed", "enhancers_h3k27ac_control", "H3K27ac Enhancers - Control", "H3K27ac enhancer signal in control samples lifted to hg38", "0,102,204"),
]


def parse_bed(path: Path):
    with path.open() as fh:
        header = fh.readline().rstrip("\n")
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 6:
                cols = line.split()
            if len(cols) < 6:
                continue
            yield cols


def liftover_interval(lo: LiftOver, chrom: str, start: int, end: int):
    chrom = chrom.replace("chr", "")

    s = []
    for pos in range(start, start + 200):
        s = lo.convert_coordinate(f"chr{chrom}", pos)
        if s:
            break

    e = []
    for pos in range(end - 1, max(start, end - 200), -1):
        e = lo.convert_coordinate(f"chr{chrom}", pos)
        if e:
            break

    if not s or not e:
        return None

    s_chrom, s_pos = s[0][0], s[0][1]
    e_chrom, e_pos = e[0][0], e[0][1]
    if s_chrom != e_chrom:
        return None

    if s_chrom.startswith("chr"):
        out_chrom = s_chrom
    else:
        out_chrom = f"chr{s_chrom}"
    new_start = min(s_pos, e_pos)
    new_end = max(s_pos, e_pos) + 1
    return out_chrom, new_start, new_end


def make_track_header(name: str, description: str, color: str) -> str:
    return f'track name="{name}" description="{description}" color={color} visibility=dense'


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not CHAIN_FILE.exists():
        print(f"Downloading chain file: {CHAIN_URL}")
        urllib.request.urlretrieve(CHAIN_URL, CHAIN_FILE)
    lo = LiftOver(str(CHAIN_FILE))
    print(f"Loaded hg19 -> hg38 liftover from {CHAIN_FILE}")

    for src_name, short_name, label, desc, color in TRACKS:
        src = SRC_DIR / src_name
        out_bed = OUT_DIR / src_name
        out_sorted = OUT_DIR / src_name.replace('.bed', '.sorted.bed')
        out_bb = OUT_DIR / src_name.replace('.bed', '.bb')

        lifted_rows = []
        total = 0
        skipped = 0
        for cols in parse_bed(src):
            total += 1
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            lifted = liftover_interval(lo, chrom, start, end)
            if not lifted:
                skipped += 1
                continue
            new_chrom, new_start, new_end = lifted
            if new_chrom not in STANDARD_CHROMS:
                skipped += 1
                continue
            cols[0], cols[1], cols[2] = new_chrom, str(new_start), str(new_end)
            lifted_rows.append(cols)

        lifted_rows.sort(key=lambda r: (r[0], int(r[1]), int(r[2])))
        with out_bed.open("w") as fh:
            fh.write(make_track_header(label, desc, color) + "\n")
            for row in lifted_rows:
                fh.write("\t".join(row) + "\n")

        # For bigBed, remove the track line first
        with out_sorted.open("w") as fh:
            for row in lifted_rows:
                fh.write("\t".join(row[:6]) + "\n")

        subprocess.run([str(BEDTOBIGBED), str(out_sorted), str(CHROMSIZES), str(out_bb)], check=True)
        print(f"{src_name}: kept {len(lifted_rows):,}/{total:,}, skipped {skipped:,}")

    (SRC_DIR / "hub_h3k27ac_hg38.txt").write_text(
        "hub H3K27acEnhancers_hg38\n"
        "shortLabel H3K27ac Enhancers hg38\n"
        "longLabel GSE102538 H3K27ac enhancer tracks lifted to hg38 (AD vs Control)\n"
        "genomesFile genomes_h3k27ac_hg38.txt\n"
        "email your_email@example.com\n"
        "descriptionUrl https://github.com/Dariasels/Epigenetic_layers_integrated_portfolio\n"
    )
    (SRC_DIR / "genomes_h3k27ac_hg38.txt").write_text(
        "genome hg38\n"
        "trackDb trackDb_h3k27ac_hg38.txt\n"
    )
    (SRC_DIR / "trackDb_h3k27ac_hg38.txt").write_text(
        "# H3K27ac enhancer tracks from GSE102538 lifted to hg38\n"
        "# Genome: hg38\n\n"
        "track enhancers_h3k27ac\n"
        "shortLabel H3K27ac Enhancers - All\n"
        "longLabel H3K27ac ChIP-seq enhancers - all samples averaged\n"
        "type bigBed 6\n"
        "bigDataUrl h3k27ac_hg38/enhancers_h3k27ac.bb\n"
        "color 255,102,0\n"
        "visibility dense\n"
        "priority 1\n\n"
        "track enhancers_h3k27ac_ad\n"
        "shortLabel H3K27ac Enhancers - AD\n"
        "longLabel H3K27ac enhancer signal in Alzheimer's samples\n"
        "type bigBed 6\n"
        "bigDataUrl h3k27ac_hg38/enhancers_h3k27ac_ad.bb\n"
        "color 204,0,0\n"
        "visibility dense\n"
        "priority 1.1\n\n"
        "track enhancers_h3k27ac_control\n"
        "shortLabel H3K27ac Enhancers - Control\n"
        "longLabel H3K27ac enhancer signal in control samples\n"
        "type bigBed 6\n"
        "bigDataUrl h3k27ac_hg38/enhancers_h3k27ac_control.bb\n"
        "color 0,102,204\n"
        "visibility dense\n"
        "priority 1.2\n"
    )


if __name__ == "__main__":
    main()
