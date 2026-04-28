#!/usr/bin/env python3
from pathlib import Path
import subprocess

ROOT = Path('/home/daria/Epigenetic_layers_integrated_portfolio')
OUTDIR = ROOT / 'outputs' / 'ucsc_tracks'
BEDTOBIGBED = Path('/tmp/bedToBigBed')
CHROMSIZES = Path('/tmp/hg38.chrom.sizes')

order = {}
sizes = {}
with CHROMSIZES.open() as f:
    for i, line in enumerate(f):
        chrom, size = line.split()[:2]
        order[chrom] = i
        sizes[chrom] = int(size)

tracks = [
    'plasticity_genes',
    'atac_peaks',
    'enhancers_h3k27ac',
    'methylation_promoters',
]

for track in tracks:
    src = OUTDIR / f'{track}.bed'
    clean = OUTDIR / f'{track}.clean.bed6'
    bb = OUTDIR / f'{track}.bb'

    rows = []
    kept = 0
    dropped = 0

    with src.open() as f:
        for line in f:
            if not line.strip() or line.startswith('track '):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6:
                dropped += 1
                continue
            chrom, start_s, end_s = parts[0], parts[1], parts[2]
            if chrom not in sizes:
                dropped += 1
                continue
            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                dropped += 1
                continue
            if start < 0 or end > sizes[chrom] or start >= end:
                dropped += 1
                continue
            rows.append((order[chrom], chrom, start, end, parts[3], parts[4], parts[5]))
            kept += 1

    rows.sort(key=lambda x: (x[0], x[2], x[3]))
    with clean.open('w') as out:
        for _, chrom, start, end, name, score, strand in rows:
            out.write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n')

    subprocess.run([
        str(BEDTOBIGBED),
        str(clean),
        str(CHROMSIZES),
        str(bb),
    ], check=True)

    print(f'{track}: kept {kept:,}, dropped {dropped:,}, wrote {bb.name}')

print('Done')
