#!/usr/bin/env python3
"""Build UCSC RNA condition tracks from gene-level summary data.

This exports AD and control mean expression directly from the MySQL summary
tables by joining rna_expression -> samples -> gene_coordinates.

The result is a pair of bar-style bigWig tracks plus bedGraph companions for
inspection and debugging.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import mysql.connector
import pyBigWig

ROOT = Path('/home/daria/Epigenetic_layers_integrated_portfolio')
OUTDIR = ROOT / 'outputs' / 'ucsc_tracks' / 'rna_hg38'
CHROMSIZES = ROOT / 'tools' / 'juicebox' / 'tools' / 'chrom' / 'sizes' / 'hg38.chrom.sizes'
BEDTOBIGBED = Path('/tmp/bedToBigBed')
BEDTOBIGBED_FALLBACK = Path('/home/daria/anaconda3/bin/bedToBigBed')

DB_CONFIG = {
    'host': 'localhost',
    'user': 'daria',
    'password': 'simba',
    'database': 'brain_multiomics',
}

TRACKS = {
    'ad': {
        'bw_name': 'rna_expression_ad_signal.bw',
        'bedgraph_name': 'rna_expression_ad_signal.bedGraph',
        'short_label': 'RNA Expression AD',
        'long_label': "Mean RNA microarray expression in Alzheimer's samples (gene-level)",
        'color': '204,0,0',
    },
    'control': {
        'bw_name': 'rna_expression_control_signal.bw',
        'bedgraph_name': 'rna_expression_control_signal.bedGraph',
        'short_label': 'RNA Expression Control',
        'long_label': 'Mean RNA microarray expression in control samples (gene-level)',
        'color': '0,102,204',
    },
}

PLASTICITY_LOG2FC_BED = OUTDIR / 'plasticity_genes_rna_log2fc.bed'
PLASTICITY_LOG2FC_BED6 = OUTDIR / 'plasticity_genes_rna_log2fc.clean.bed'
PLASTICITY_LOG2FC_BB = OUTDIR / 'plasticity_genes_rna_log2fc.bb'


def connect():
    return mysql.connector.connect(**DB_CONFIG)


def load_chrom_sizes() -> tuple[dict[str, int], dict[str, int]]:
    sizes: dict[str, int] = {}
    order: dict[str, int] = {}
    with CHROMSIZES.open() as fh:
        for idx, line in enumerate(fh):
            chrom, size = line.split()[:2]
            sizes[chrom] = int(size)
            order[chrom] = idx
    return sizes, order


def fetch_summary_rows() -> list[tuple[str, int, int, str, str, float, int]]:
    conn = connect()
    try:
        cursor = conn.cursor()
        cursor.execute(
            """
            SELECT
                gc.chrom,
                gc.start_pos,
                gc.end_pos,
                rs.gene_symbol,
                rs.condition_key,
                rs.mean_expression,
                rs.n_samples
            FROM (
                SELECT
                    re.gene_symbol,
                    CASE
                        WHEN LOWER(s.condition) IN ('ad', 'alzheimer') THEN 'ad'
                        WHEN LOWER(s.condition) = 'control' THEN 'control'
                        ELSE NULL
                    END AS condition_key,
                    AVG(re.expression_value) AS mean_expression,
                    COUNT(DISTINCT re.sample_id) AS n_samples
                FROM rna_expression re
                INNER JOIN samples s ON s.sample_id = re.sample_id
                WHERE LOWER(s.condition) IN ('ad', 'alzheimer', 'control')
                GROUP BY re.gene_symbol, condition_key
                HAVING condition_key IS NOT NULL
            ) rs
            INNER JOIN gene_coordinates gc ON gc.gene_symbol = rs.gene_symbol
            WHERE gc.chrom NOT IN ('MT', 'chrM', 'Un_%')
            AND NOT EXISTS (
                SELECT 1
                FROM gene_coordinates gc2
                WHERE gc2.gene_symbol = gc.gene_symbol
                  AND gc2.chrom NOT IN ('MT', 'chrM', 'Un_%')
                  AND (
                    gc2.start_pos < gc.start_pos
                    OR (gc2.start_pos = gc.start_pos AND gc2.end_pos < gc.end_pos)
                  )
            )
            ORDER BY gc.chrom, gc.start_pos, rs.gene_symbol, rs.condition_key
            """
        )
        rows = cursor.fetchall()
        cursor.close()
        return rows
    finally:
        conn.close()


def fetch_plasticity_log2fc_rows() -> list[tuple[str, int, int, str, float | None, float | None]]:
    conn = connect()
    try:
        cursor = conn.cursor()
        cursor.execute(
            """
            SELECT
                gc.chrom,
                gc.start_pos,
                gc.end_pos,
                ps.gene_symbol,
                ps.rna_AD,
                ps.rna_Ctrl,
                ps.rna_delta
            FROM plasticity_summary ps
            INNER JOIN gene_coordinates gc ON gc.gene_symbol = ps.gene_symbol
              AND gc.chrom NOT IN ('MT', 'chrM', 'Un_%')
            AND NOT EXISTS (
                SELECT 1
                FROM gene_coordinates gc2
                WHERE gc2.gene_symbol = gc.gene_symbol
                  AND gc2.chrom NOT IN ('MT', 'chrM', 'Un_%')
                  AND (
                    gc2.start_pos < gc.start_pos
                    OR (gc2.start_pos = gc.start_pos AND gc2.end_pos < gc.end_pos)
                  )
            )
            ORDER BY gc.chrom, gc.start_pos, ps.gene_symbol
            """
        )
        rows = cursor.fetchall()
        cursor.close()
        return rows
    finally:
        conn.close()


def normalize_chrom(chrom: str) -> str:
    return chrom if chrom.startswith('chr') else f'chr{chrom}'


def write_plasticity_log2fc_labels(sizes: dict[str, int], order: dict[str, int]) -> None:
    rows = fetch_plasticity_log2fc_rows()

    entries: list[tuple[str, int, int, str, int]] = []
    for chrom, start, end, gene_symbol, ad_mean, control_mean, rna_delta in rows:
        chrom = normalize_chrom(str(chrom))
        if chrom not in sizes:
            continue
        start = int(start)
        end = int(end)
        if start < 0 or end <= start:
            continue

        if rna_delta is None:
            label = f'{gene_symbol}:NA'
            score = 100
        else:
            label = f'{gene_symbol}:{float(rna_delta):.3f}'
            score = min(999, max(100, int(abs(float(rna_delta)) * 250)))

        entries.append((chrom, start, end, label, score))

    entries.sort(key=lambda item: (order[item[0]], item[1], item[3]))
    with PLASTICITY_LOG2FC_BED.open('w') as out:
        out.write(
            'track name="Plasticity_RNA_log2FC" '
            'description="Plasticity genes with RNA log2FC values (AD-control) labels" '
            'color=120,80,200 visibility=pack\n'
        )
        for chrom, start, end, label, score in entries:
            out.write(f'{chrom}\t{start}\t{end}\t{label}\t{score}\t+\n')

    with PLASTICITY_LOG2FC_BED6.open('w') as out:
        for chrom, start, end, label, score in entries:
            out.write(f'{chrom}\t{start}\t{end}\t{label}\t{score}\t+\n')

    bed_to_big_bed = BEDTOBIGBED if BEDTOBIGBED.exists() else BEDTOBIGBED_FALLBACK
    if bed_to_big_bed.exists():
        subprocess.run([
            str(bed_to_big_bed),
            str(PLASTICITY_LOG2FC_BED6),
            str(CHROMSIZES),
            str(PLASTICITY_LOG2FC_BB),
        ], check=True)
    else:
        raise FileNotFoundError(
            f'missing bedToBigBed converter at {BEDTOBIGBED} or {BEDTOBIGBED_FALLBACK}'
        )

    print(
        f'wrote {PLASTICITY_LOG2FC_BED.name} and {PLASTICITY_LOG2FC_BB.name} '
        f'({len(entries):,} plasticity genes with labels)'
    )


def write_tracks() -> None:
    sizes, order = load_chrom_sizes()
    OUTDIR.mkdir(parents=True, exist_ok=True)

    rows = fetch_summary_rows()
    by_condition: dict[str, list[tuple[str, int, int, str, float, int]]] = {'ad': [], 'control': []}
    for chrom, start, end, gene_symbol, condition, mean_expression, n_samples in rows:
        chrom = normalize_chrom(str(chrom))
        if chrom not in sizes:
            continue
        start = int(start)
        end = int(end)
        if start < 0 or end <= start:
            continue
        if condition not in by_condition:
            continue
        by_condition[condition].append((chrom, start, end, gene_symbol, float(mean_expression), int(n_samples)))

    for condition, config in TRACKS.items():
        entries = by_condition[condition]
        entries.sort(key=lambda item: (order[item[0]], item[1], item[3]))

        bw_path = OUTDIR / config['bw_name']
        bedgraph_path = OUTDIR / config['bedgraph_name']

        bw = pyBigWig.open(str(bw_path), 'w')
        try:
            bw.addHeader(list(sizes.items()))
            if entries:
                chroms = [chrom for chrom, _, _, _, _, _ in entries]
                starts = [start for _, start, _, _, _, _ in entries]
                ends = [end for _, _, end, _, _, _ in entries]
                values = [mean_expression for _, _, _, _, mean_expression, _ in entries]
                bw.addEntries(chroms, starts, ends=ends, values=values)
        finally:
            bw.close()

        with bedgraph_path.open('w') as out:
            out.write(
                f'track type=bedGraph name="{config["short_label"].replace(" ", "_")}" '
                f'description="{config["long_label"]}" color={config["color"]}\n'
            )
            for chrom, start, end, gene_symbol, mean_expression, n_samples in entries:
                out.write(f'{chrom}\t{start}\t{end}\t{mean_expression:.4f}\n')

        print(f'wrote {bw_path.name} and {bedgraph_path.name} ({len(entries):,} gene intervals)')

    write_plasticity_log2fc_labels(sizes, order)


if __name__ == '__main__':
    write_tracks()