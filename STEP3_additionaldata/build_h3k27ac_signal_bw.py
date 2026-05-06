#!/usr/bin/env python3
"""
Build signal-style bigWig tracks from the lifted hg38 H3K27ac enhancer BEDs.
Uses the BED score column as the signal value.
"""
from __future__ import annotations

from pathlib import Path
import pyBigWig

ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT / "outputs" / "ucsc_tracks" / "h3k27ac_hg38"
CHROMSIZES = ROOT / "tools" / "juicebox" / "tools" / "chrom" / "sizes" / "hg38.chrom.sizes"

TRACKS = [
    ("enhancers_h3k27ac.bed", "enhancers_h3k27ac_signal.bw", "H3K27ac Enhancers - All", "H3K27ac enhancer signal (all samples, hg38)", "255,102,0"),
    ("enhancers_h3k27ac_ad.bed", "enhancers_h3k27ac_ad_signal.bw", "H3K27ac Enhancers - AD", "H3K27ac enhancer signal in Alzheimer's samples (hg38)", "204,0,0"),
    ("enhancers_h3k27ac_control.bed", "enhancers_h3k27ac_control_signal.bw", "H3K27ac Enhancers - Control", "H3K27ac enhancer signal in control samples (hg38)", "0,102,204"),
]

TRACKDB_FILE = SRC_DIR.parent / "trackDb_h3k27ac_hg38.txt"


def load_chrom_sizes() -> dict[str, int]:
    sizes: dict[str, int] = {}
    with CHROMSIZES.open() as fh:
        for line in fh:
            chrom, size = line.strip().split()[:2]
            sizes[chrom] = int(size)
    return sizes


def parse_bed(path: Path):
    with path.open() as fh:
        first = fh.readline()
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                cols = line.split()
            if len(cols) < 6:
                continue
            yield cols


def write_bigwig(in_bed: Path, out_bw: Path, chrom_sizes: dict[str, int]):
    events: dict[str, list[tuple[int, int, float]]] = {}
    for cols in parse_bed(in_bed):
        chrom = cols[0]
        if chrom not in chrom_sizes:
            continue
        start = int(cols[1])
        end = int(cols[2])
        if end <= start:
            continue
        try:
            value = float(cols[4])
        except Exception:
            value = 0.0
        # +1 at start, -1 at end; store the value separately
        events.setdefault(chrom, []).append((start, 1, value))
        events[chrom].append((end, -1, value))

    bw = pyBigWig.open(str(out_bw), "w")
    bw.addHeader([(chrom, chrom_sizes[chrom]) for chrom in sorted(chrom_sizes.keys())])

    for chrom in sorted(events.keys()):
        if chrom not in chrom_sizes:
            continue
        chrom_events = sorted(events[chrom], key=lambda x: (x[0], -x[1]))
        active: list[float] = []
        cur_pos = None
        starts: list[int] = []
        ends: list[int] = []
        values: list[float] = []

        for pos, kind, value in chrom_events:
            if cur_pos is not None and pos > cur_pos and active:
                starts.append(cur_pos)
                ends.append(pos)
                values.append(sum(active) / len(active))
            cur_pos = pos
            if kind == 1:
                active.append(value)
            else:
                # remove one matching value
                try:
                    active.remove(value)
                except ValueError:
                    pass

        if starts:
            bw.addEntries([chrom] * len(starts), starts, ends=ends, values=values)
    bw.close()


def main():
    chrom_sizes = load_chrom_sizes()
    for in_name, out_name, _, _, _ in TRACKS:
        in_bed = SRC_DIR / in_name
        out_bw = SRC_DIR / out_name
        write_bigwig(in_bed, out_bw, chrom_sizes)
        print(f"wrote {out_bw}")

    trackdb = TRACKDB_FILE.read_text()
    if "enhancers_h3k27ac_ad_signal" not in trackdb:
        trackdb += (
            "\ntrack enhancers_h3k27ac_signal\n"
            "shortLabel H3K27ac Enhancers - Signal\n"
            "longLabel H3K27ac enhancer signal (all samples, hg38)\n"
            "type bigWig\n"
            "bigDataUrl h3k27ac_hg38/enhancers_h3k27ac_signal.bw\n"
            "color 255,102,0\n"
            "visibility full\n"
            "graphTypeDefault bar\n"
            "viewLimits 0:400\n"
            "autoScale on\n"
            "windowingFunction mean\n"
            "priority 0.9\n\n"
            "track enhancers_h3k27ac_ad_signal\n"
            "shortLabel H3K27ac AD Signal\n"
            "longLabel H3K27ac enhancer signal in Alzheimer's samples (hg38)\n"
            "type bigWig\n"
            "bigDataUrl h3k27ac_hg38/enhancers_h3k27ac_ad_signal.bw\n"
            "color 204,0,0\n"
            "visibility full\n"
            "graphTypeDefault bar\n"
            "viewLimits 0:400\n"
            "autoScale on\n"
            "windowingFunction mean\n"
            "priority 1.05\n\n"
            "track enhancers_h3k27ac_control_signal\n"
            "shortLabel H3K27ac Control Signal\n"
            "longLabel H3K27ac enhancer signal in control samples (hg38)\n"
            "type bigWig\n"
            "bigDataUrl h3k27ac_hg38/enhancers_h3k27ac_control_signal.bw\n"
            "color 0,102,204\n"
            "visibility full\n"
            "graphTypeDefault bar\n"
            "viewLimits 0:400\n"
            "autoScale on\n"
            "windowingFunction mean\n"
            "priority 1.1\n"
        )
        TRACKDB_FILE.write_text(trackdb)
        print("updated trackDb_h3k27ac_hg38.txt")

    hub = SRC_DIR / "hub_h3k27ac_hg38.txt"
    if hub.exists():
        print(f"hub file: {hub}")


if __name__ == "__main__":
    main()
