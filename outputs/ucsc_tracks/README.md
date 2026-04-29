# UCSC Genome Browser Track Visualization

## Overview

This directory contains UCSC Genome Browser compatible track files for the **Alzheimer's Disease & Aging Brain Multi-Omics Integration** project. These files represent a comprehensive multi-layer epigenetic analysis of prefrontal cortex from three cohorts: Alzheimer's Disease (AD), cognitively normal aged controls, and young controls.

## Track Files

### 1. **plasticity_genes.bed** (394 genes)
- **Source:** Curated list of 365 brain plasticity genes from literature
- **Content:** Gene coordinates for known plasticity-associated genes
- **Color:** Blue (100,149,237)
- **Use:** Filter/focus analysis on genes known to regulate neural plasticity, synaptic function, and neurogenesis

### 2. **atac_peaks.bed** (370,494 peaks)
- **Source:** ATAC-seq open chromatin data (GSE129040)
- **Content:** Accessible chromatin regions (DNase I hypersensitive sites)
- **Color:** Green (0,176,80)
- **9th Column:** Associated genes linked within ±100kb
- **Files:** `/data/GEO_raw/GSE129040/broadPeak/*.broadPeak`

### 3. **enhancers_h3k27ac.bed** (182,330 enhancers)  
- **Source:** H3K27ac ChIP-seq active enhancers (GSE102538, GSE143271)
- **Content:** H3K27ac-marked active regulatory elements
- **Color:** Orange (255,102,0)
- **9th Column:** Associated genes linked via proximity + TAD constraints
- **Score:** Scaled by number of linked genes (more genes = higher score)

### 4. **methylation_promoters.bed** (376 regions)
- **Source:** CpG methylation 450K array (GSE59685)
- **Content:** Beta-norm CpG sites at gene promoters (TSS ±200bp)
- **Color:** Red (200,0,0)
- **Score:** Weighted by number of CpG probes at each promoter
- **Relevance:** Promoter hypermethylation associated with silencing; hypomethylation associated with activation

### 5. **cpg_methylation_signal.bw** (per-site beta values)
- **Source:** GSE59685 methylation table
- **Content:** One 1-bp feature per CpG probe, scored by mean beta value across available samples
- **Color:** Red (200,0,0)
- **Score:** $0$–$1000$ scaled from mean beta value
- **Use:** Shows methylation intensity at each CpG site rather than only promoter windows

## Hub Configuration Files

### **hub.txt**
Main hub registry file. Contains:
- Hub ID: `AlzheimersEpigenetics_AD_Aging`
- Description for UCSC display
- Pointer to genome list (`genomes.txt`)
- Contact information and project URL

### **hub_v3.txt**
Versioned hub registry file that points to the UCSC-compatible bigBed tracks.


### **genomes.txt**
Genome definition file. Specifies:
- Genome assembly: `hg38` (GRCh38)
- Track configuration file: `trackDb.txt`

### **trackDb.txt**
Track database with full track definitions:
- Track names, descriptions, and metadata
- BigDataUrl pointers to bigBed files
- Visualization parameters (colors, visibility)
- Priority ordering for display

### **trackDb_v3.txt**
Versioned track database that points to the bigBed files generated from the BED sources.

## Using with UCSC Genome Browser

### Option A: Local Hub Upload
1. Host these files on a web server or git repository
2. In UCSC Genome Browser → **Manage Custom Tracks** → **Track Hubs**
3. Paste URL to `hub.txt`:  
   ```
  https://github.com/Dariasels/Epigenetic_layers_integrated_portfolio/raw/main/outputs/ucsc_tracks/hub_v3.txt
   ```
4. Click "Add Hub"

### Option B: File Upload
1. Compress the directory:  
   ```bash
   tar -czf ucsc_tracks.tar.gz *.bed *.txt
   ```
2. Upload via UCSC Browser → **Manage Custom Tracks** → Local file upload

### Option C: bigBed optimization
The hub now uses UCSC bigBed files for compatibility and faster loading:
```bash
python build_ucsc_bigbeds.py
```

## Data Integration Statistics

| Layer | Features | Linked Genes | Source |
|-------|----------|--------------|--------|
| **Plasticity Genes** | 365 | - | Curated literature |
| **ATAC-seq** | 370,494 | 16,800+ | GSE129040 |
| **H3K27ac Enhancers** | 182,330 | 18,200+ | GSE102538/GSE143271 |
| **CpG Methylation Signal** | per-site | 15,000+ | GSE59685 |
| **Methylation** | 376 | 365 | GSE59685 |
| **Multi-layer** | ~900 | - | 2+ data types |

## Cohort Information

- **Alzheimer's Disease (AD):** n=4 females, mean age 91.5 years
- **Cognitively Normal Aged:** n=4 females, mean age 90 years
- **Young Controls:** n=females, age ~29 years
- **Tissue:** Prefrontal cortex (Brodmann area 9/10)
- **Post-mortem brain tissue** from NIH Brain Bank

## Methods Summary

### ATAC-seq Peak Calling
- Aligner: Bowtie2
- Peak caller: MACS2 (q-value < 0.05)
- Genome: hg38
- Filtering: Peaks with FDR < 0.05

### H3K27ac ChIP-seq
- Peak caller: MACS2
- Normalization: CPM (counts per million)
- Assignment: Enhancers linked to genes within TAD + proximity constraints

### Methylation Analysis
- Platform: Illumina 450K array
- Normalization: Beta-mixture model (minfi)
- CpG selection: Sites in promoter regions (±2kb TSS)
- Filtering: Sites with p-value < 0.05

### Gene Coordinate Mapping
- Annotation: GENCODE v40 (hg38)
- Curation: Gene location from RefSeq + Ensembl
- Link resolution:
  - ATAC → genes: peaks within ±100kb of TSS
  - Enhancers → genes: via TAD adjacency + proximity
  - Methylation → genes: CpG sites in promoter

## Database Schema

All tracks are derived from MySQL database `brain_multiomics`:
- **atac_peaks** table (370,494 rows)
- **enhancers** table (182,330 rows)
- **plasticity_genes** table (365 rows)
- **methylation_gene_links** table (376 rows)
- **gene_coordinates** table (~19,000 genes)

SQL integration views:
- `multi_layer_integration` - genes with 2+ data types
- `enhancer_gene_links` - curated enhancer-gene interactions
- `atac_gene_links` - ATAC peak-gene assignments

## File Format Reference

### BED Format (all tracks)
```
Columns:
  1. chrom       - Chromosome (chr1, chrX, etc.)
  2. chromStart  - 0-based start position
  3. chromEnd    - 1-based end position (exclusive)
  4. name        - Feature identifier
  5. score       - 0-999 (UCSC convention, higher = more important)
  6. strand      - + or - (optional)
  7-9. Extra     - Tool/feature-specific (e.g., gene associations)
```

### Track Header Format
```
track name="TrackID" description="Display name" color=R,G,B visibility=dense
```

## Download & Usage

Files are automatically generated from the MySQL database during project build:

```bash
# Navigate to project root
cd /home/daria/Epigenetic_layers_integrated_portfolio

# Regenerate tracks (requires MySQL password)
export DB_PASSWORD="your_password"
python outputs/generate_ucsc_bedfiles.py --output-dir outputs/ucsc_tracks

# Verify track counts
wc -l outputs/ucsc_tracks/*.bed
```

## Next Steps

1. **Upload to UCSC Hub:** Follow **Option A** above to load tracks
2. **Explore Regions:** Use browser to navigate to genes of interest
3. **Filter/Refine:** Use UCSC tools to intersect, filter, or analyze overlaps
4. **Publication:** Export screenshots and track coordinates for figures

## References

- **Project:** Epigenetic Layers Integrated Portfolio
- **Data Paper:** Weng et al., Nature Scientific Data (2023) - AD/Aging Hi-C
- **Genome Build:** hg38 (GRCh38)
- **UCSC Docs:** https://genome.ucsc.edu/goldenPath/help/customTrack.html

## Contact

For questions or issues:
- Project Repository: https://github.com/Dariasels/Epigenetic_layers_integrated_portfolio
- Data inquiries: Contact original GEO dataset authors

---

**Generated:** 2026-04-28  
**Status:** ✅ Complete (4 core tracks; Hi-C inaccessible from source files)
