#!/usr/bin/env python3
from pathlib import Path
import csv
import subprocess

ROOT = Path('/home/daria/Epigenetic_layers_integrated_portfolio')
ANNOT = ROOT / 'STEP2_mysql_import' / 'GPL13534_annotation.csv'
OUTDIR = ROOT / 'outputs' / 'ucsc_tracks'
BED = OUTDIR / 'cpg_probes.bed'
BB = OUTDIR / 'cpg_probes.bb'
CHROMSIZES = Path('/tmp/hg38.chrom.sizes')
BEDTOBIGBED = Path('/tmp/bedToBigBed')

sizes = {}
with CHROMSIZES.open() as f:
    for line in f:
        chrom, size = line.split()[:2]
        sizes[chrom] = int(size)

kept = 0
with ANNOT.open(newline='') as f, BED.open('w') as out:
    reader = csv.DictReader(f)
    out.write('track name="CpG_Probes" description="Illumina 450K CpG probes" color=120,120,120 visibility=dense\n')
    for row in reader:
        probe = row['Name'].strip().strip('"')
        chrom = row['chr'].strip().strip('"')
        pos = row['pos'].strip().strip('"')
        if not probe.startswith('cg'):
            continue
        if chrom not in sizes:
            continue
        try:
            p = int(pos)
        except ValueError:
            continue
        if p <= 0 or p > sizes[chrom]:
            continue
        start = p - 1
        end = p
        out.write(f'{chrom}\t{start}\t{end}\t{probe}\t500\t+\n')
        kept += 1

# sort and convert
sorted_bed = OUTDIR / 'cpg_probes.sorted.bed6'
with BED.open() as src, sorted_bed.open('w') as dst:
    lines = [l for l in src if l.strip() and not l.startswith('track ')]
    lines.sort(key=lambda l: (l.split('\t')[0], int(l.split('\t')[1])))
    dst.writelines(lines)

subprocess.run([str(BEDTOBIGBED), str(sorted_bed), str(CHROMSIZES), str(BB)], check=True)
print(f'kept {kept:,} CpG probes')
