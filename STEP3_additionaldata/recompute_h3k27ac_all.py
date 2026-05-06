#!/usr/bin/env python3
"""
Recompute 'All' enhancer track as weighted mean of AD and Control.
Uses n_AD=24, n_CTRL=23 for proper weighting.
Then regenerates bigWig tracks with normalized viewLimits.
"""
from pathlib import Path
import pymysql

ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT / "outputs" / "ucsc_tracks" / "h3k27ac_hg38"
CHROMSIZES = ROOT / "tools" / "juicebox" / "tools" / "chrom" / "sizes" / "hg38.chrom.sizes"

N_AD = 24
N_CTRL = 23
TOTAL_N = N_AD + N_CTRL

def load_chrom_sizes():
    sizes = {}
    with CHROMSIZES.open() as fh:
        for line in fh:
            chrom, size = line.strip().split()[:2]
            sizes[chrom] = int(size)
    return sizes

def parse_bed(path):
    """Parse BED file, skip track header, yield data columns."""
    with path.open() as fh:
        for line in fh:
            if line.startswith('track'):
                continue
            if not line.strip():
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 6:
                cols = line.split()
            if len(cols) < 6:
                continue
            yield cols

def recompute_all_bed():
    """
    Read AD and Control BEDs, compute weighted mean for All, write new .bed file.
    """
    ad_bed = SRC_DIR / "enhancers_h3k27ac_ad.bed"
    ctrl_bed = SRC_DIR / "enhancers_h3k27ac_control.bed"
    all_bed = SRC_DIR / "enhancers_h3k27ac.bed"
    
    # Load AD and Control data indexed by (chrom, start, end, name)
    ad_data = {}
    for cols in parse_bed(ad_bed):
        chrom, start, end, name = cols[0], int(cols[1]), int(cols[2]), cols[3]
        try:
            score = float(cols[4])
        except:
            score = 0.0
        key = (chrom, start, end)
        ad_data[key] = (name, score)
    
    ctrl_data = {}
    for cols in parse_bed(ctrl_bed):
        chrom, start, end, name = cols[0], int(cols[1]), int(cols[2]), cols[3]
        try:
            score = float(cols[4])
        except:
            score = 0.0
        key = (chrom, start, end)
        ctrl_data[key] = (name, score)
    
    # Merge and compute weighted averages
    all_data = []
    for key in sorted(set(ad_data.keys()) | set(ctrl_data.keys())):
        chrom, start, end = key
        ad_name, ad_score = ad_data.get(key, ("", 0.0))
        ctrl_name, ctrl_score = ctrl_data.get(key, ("", 0.0))
        
        # Compute weighted mean
        weighted_score = (ad_score * N_AD + ctrl_score * N_CTRL) / TOTAL_N
        
        # Extract base name (remove _AD or _CTRL suffix)
        base_name = ad_name or ctrl_name
        base_name = base_name.split('_')[0] if '_' in base_name else base_name
        if '|' in base_name:
            base_name = base_name.split('|')[0]
        
        all_data.append((chrom, start, end, f"{base_name}|AVG={weighted_score:.0f}", weighted_score, '+'))
    
    # Write to file
    with all_bed.open('w') as fh:
        fh.write('track name="H3K27ac Enhancers - All" description="H3K27ac enhancer signal (all samples, weighted) lifted to hg38" color=255,102,0 visibility=dense\n')
        for chrom, start, end, name, score, strand in all_data:
            fh.write(f"{chrom}\t{start}\t{end}\t{name}\t{int(score)}\t{strand}\n")
    
    print(f"Wrote {len(all_data)} enhancers to {all_bed}")
    return all_data

def compute_global_percentile_99(enhancers_data):
    """Compute 99th percentile of all enhancer scores."""
    scores = [score for _, _, _, _, score, _ in enhancers_data]
    
    # Also read AD and Control scores
    ad_bed = SRC_DIR / "enhancers_h3k27ac_ad.bed"
    ctrl_bed = SRC_DIR / "enhancers_h3k27ac_control.bed"
    
    for cols in parse_bed(ad_bed):
        try:
            scores.append(float(cols[4]))
        except:
            pass
    for cols in parse_bed(ctrl_bed):
        try:
            scores.append(float(cols[4]))
        except:
            pass
    
    scores.sort()
    p99_idx = int(len(scores) * 0.99)
    p99_val = scores[min(p99_idx, len(scores) - 1)]
    
    print(f"99th percentile of all scores: {p99_val:.1f} (from {len(scores)} total values)")
    return p99_val

if __name__ == "__main__":
    print("Recomputing 'All' enhancer track as weighted mean (n_AD=24, n_CTRL=23)...")
    all_enhancers = recompute_all_bed()
    
    p99 = compute_global_percentile_99(all_enhancers)
    
    # Next step: regenerate bigWigs using build_h3k27ac_signal_bw.py
    print(f"\nNext: run build_h3k27ac_signal_bw.py to regenerate bigWigs")
    print(f"Then: update trackDb_h3k27ac_hg38.txt with viewLimits 0:{int(p99)}")
