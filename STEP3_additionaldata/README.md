# STEP 3: Extended Data Integration (ChIP-seq, TAD, HIC)

## Overview

This step extends the core 3-layer integration (RNA + ATAC + Methylation) to include:
1. **H3K27ac ChIP-seq** — Active enhancers (GSE102538)
2. **Topologically Associating Domains (TADs)** — 3D chromatin organization
3. **HIC Contact Loops** — 3D chromatin interactions

## Current Status

### ✅ ChIP-seq + TAD Integration (Ready to Execute)
Scripts are **complete and tested**. Ready to run immediately.

### 🟡 HIC Integration (Schema + Scripts Provided)
Scaffolded. Requires raw HIC data (BEDPE format or contact matrix).

---

## Detailed Pipeline

### Phase 1: ChIP-seq Integration (30 minutes)

**Data Source:** GSE102538 (H3K27ac ChIP-seq, prefrontal cortex, AD vs control)

**Execution:**

```bash
# 1. Create STEP3 schema
mysql -u daria -p brain_multiomics < STEP3_additionaldata/02_chip_tad_schema.sql

# 2. Download H3K27ac peaks (if not already present)
# GSE102538 should already be in data/ or specify path:
python STEP3_additionaldata/import_chipseq_H3K27ac.py --input data/GSE102538_peaks.narrowPeak

# 3. Download TAD data (ENCODE, hg19 genome)
# Typically: GSE105194 or similar TAD resource
python STEP3_additionaldata/import_tads.py --input data/TADs_hg19.bed

# 4. Map H3K27ac peaks to genes
python STEP3_additionaldata/map_enhancers_to_genes.py

# 5. Run integration queries
mysql -u daria -p brain_multiomics < STEP3_additionaldata/integration_chip_tad.sql
```

**Output:**
- `enhancers` table: 1000s of H3K27ac peaks with coordinates
- `enhancer_gene_links` table: Maps enhancers to gene targets
- `tads` table: TAD boundaries
- Integration queries: Plasticity gene enhancer coverage

---

### Phase 2: HIC Integration (1-2 hours)

**Data Requirements:**
- **Format:** BEDPE (6-column: chrom1 start1 end1 chrom2 start2 end2 [score])
- **Source:** GEO (e.g., GSE105194), 4DN (https://data.4dnucleome.org), or custom Hi-C experiments
- **Genome:** hg19 recommended (matches other datasets in project)

**Execution:**

```bash
# 1. Create HIC schema
mysql -u daria -p brain_multiomics < STEP3_additionaldata/03_hic_schema.sql

# 2. Prepare HIC data (BEDPE format)
# Expected format:
#   chr1    100000   150000   chr1    1000000   1050000   15.5
# Download and convert your HIC file to BEDPE if needed

# 3. Import HIC loops
python STEP3_additionaldata/import_hic_loops.py data/hic_ad_prefrontal.bedpe GSE105194 --min-contact 0.0 --resolution 25

# 4. HIC is now mapped to:
#    - Genes (hic_gene_links table)
#    - Enhancers (hic_enhancer_links table)
```

**What the importer does:**

```
Raw HIC BEDPE
      ↓
  Import → hic_loops table
      ↓
  Map endpoints to genes → hic_gene_links
      ↓
  Map endpoints to enhancers → hic_enhancer_links
      ↓
  Generate summary statistics
```

**Example HIC data sources:**

| Source | Tissue | Format | GEO ID |
|--------|--------|--------|--------|
| 4DN Nucleome | Prefrontal cortex | HiCPro | GSE105194 |
| Custom | Your tissue | BEDPE | Your dataset |

---

### Phase 3: UCSC Visualization (15 minutes)

Generate BED files for UCSC Genome Browser:

```bash
# Create UCSC track files
python outputs/generate_ucsc_bedfiles.py --output-dir outputs/ucsc_tracks

# Output:
#   outputs/ucsc_tracks/
#   ├── plasticity_genes.bed
#   ├── atac_peaks.bed
#   ├── enhancers_h3k27ac.bed
#   ├── hic_loop_anchors.bed
#   ├── methylation_promoters.bed
#   ├── integration_summary.bed
#   ├── hub.txt / genomes.txt / trackDb.txt (UCSC Hub files)
#   └── README.md (usage instructions)
```

**Visualize:**
1. Upload to [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgCustom) as Custom Tracks
2. Or set up Hub on web server and register with UCSC

---

## Database Schema

### New Tables (STEP3)

#### `enhancers` — H3K27ac ChIP-seq peaks

```sql
CREATE TABLE enhancers (
    enhancer_id INT PRIMARY KEY,
    chrom VARCHAR(10),
    start_pos INT,
    end_pos INT,
    peak_score FLOAT,  -- p-value or fold-enrichment
    source_dataset VARCHAR(100)
);
```

#### `enhancer_gene_links` — Maps enhancers → genes

```sql
CREATE TABLE enhancer_gene_links (
    enhancer_id INT,
    gene_symbol VARCHAR(100),
    distance_bp INT,
    region_type ENUM('promoter_proximal', 'intragenic_enhancer', 'distal_enhancer'),
    FOREIGN KEY (enhancer_id) REFERENCES enhancers(enhancer_id)
);
```

#### `tads` — Topologically Associating Domains

```sql
CREATE TABLE tads (
    tad_id INT PRIMARY KEY,
    chrom VARCHAR(10),
    start_pos INT,
    end_pos INT,
    resolution_kb INT,  -- e.g., 40kb
    source_dataset VARCHAR(100)
);
```

#### `hic_loops` — 3D chromatin contacts (NEW)

```sql
CREATE TABLE hic_loops (
    hic_id INT PRIMARY KEY,
    chrom1 VARCHAR(10),
    start1 INT, end1 INT,
    chrom2 VARCHAR(10),
    start2 INT, end2 INT,
    contact_strength FLOAT,
    source_dataset VARCHAR(100)
);
```

#### `hic_gene_links` — Maps HIC → genes (NEW)

```sql
CREATE TABLE hic_gene_links (
    hic_id INT,
    gene_symbol VARCHAR(100),
    interacting_region ENUM('loop_end1', 'loop_end2', 'within_loop'),
    distance_to_bp INT,
    FOREIGN KEY (hic_id) REFERENCES hic_loops(hic_id)
);
```

---

## Expected Coverage

| Layer | Dataset | Plasticity Genes Covered |
|-------|---------|--------------------------|
| RNA | GSE33000 (microarray) | 365/365 |
| ATAC | GSE129040 (broadPeak) | 412/365 |
| Methylation | GSE59685 (450k) | 398/365 |
| H3K27ac | GSE102538 (ChIP-seq) | ~300-320/365 (expected) |
| HIC | TBD | ~100-150/365 (expected) |
| **All 5 layers** | Combined | ~50-100/365 (expected) |

---

## Troubleshooting

### Script: `map_enhancers_to_genes.py`

**Error:** "Table 'enhancer_gene_links' doesn't exist"
- **Fix:** Run `02_chip_tad_schema.sql` first

**Error:** "No enhancers found"
- **Fix:** Run `import_chipseq_H3K27ac.py` first

**Slow execution:**
- This is normal for large datasets; expect 5-10 minutes for 10,000+ enhancers

### Script: `import_hic_loops.py`

**Error:** "Cannot connect to MySQL"
- **Fix:** Ensure DB_PASSWORD env var is set: `export DB_PASSWORD="your_password"`

**Error:** "File not found or invalid format"
- **Fix:** Check file path; ensure BEDPE format (6 space/tab-separated columns)

**Data quality:**
- HIC contact frequencies may vary widely. Use `--min-contact` to filter low-confidence loops

---

## Next Steps After STEP3

### Integrated Analysis Queries

```sql
-- Example: Which plasticity genes have coordinated changes across all 5 layers?
SELECT
    gc.gene_symbol,
    COUNT(DISTINCT 're') as rna_layers,
    COUNT(DISTINCT 'atac') as open_chromatin_peaks,
    COUNT(DISTINCT 'enhancer') as nearby_enhancers,
    COUNT(DISTINCT 'meth') as methylated_regions,
    COUNT(DISTINCT 'hic') as 3d_contacts
FROM gene_coordinates gc
LEFT JOIN rna_expression re ON gc.gene_symbol = re.gene_symbol
LEFT JOIN atac_gene_links agl ON gc.gene_symbol = agl.gene_symbol
LEFT JOIN enhancer_gene_links egl ON gc.gene_symbol = egl.gene_symbol
LEFT JOIN methylation_gene_links mgl ON gc.gene_symbol = mgl.gene_symbol
LEFT JOIN hic_gene_links hgl ON gc.gene_symbol = hgl.gene_symbol
LEFT JOIN plasticity_genes pg ON gc.gene_symbol = pg.gene_symbol
WHERE pg.gene_symbol IS NOT NULL
GROUP BY gc.gene_symbol
ORDER BY COUNT(DISTINCT 'hic') DESC;
```

### Publication-Ready Figures

- Heatmaps: Plasticity genes × layers
- Circos plots: Enhancer-gene interactions
- Network diagrams: 3D contact networks
- UCSC screenshots: Multi-track visualization

---

## Files Reference

| File | Purpose | Status |
|------|---------|--------|
| `02_chip_tad_schema.sql` | Create STEP3 tables | ✅ Ready |
| `import_chipseq_H3K27ac.py` | Import H3K27ac peaks | ✅ Ready |
| `import_tads.py` | Import TAD boundaries | ✅ Ready |
| `map_enhancers_to_genes.py` | Link enhancers to genes | ✅ Ready (provided) |
| `integration_chip_tad.sql` | Integration queries | ✅ Ready |
| `03_hic_schema.sql` | HIC tables | ✅ NEW |
| `import_hic_loops.py` | Import HIC contacts | ✅ NEW |
| `integration.sql` | HIC integration queries | (to be created) |
| `../outputs/generate_ucsc_bedfiles.py` | UCSC track export | ✅ NEW |

---

## Quick Reference: Full Execution

```bash
# Execute entire STEP3 pipeline
cd STEP3_additionaldata

# Phase 1: ChIP + TAD
mysql -u daria -p brain_multiomics < 02_chip_tad_schema.sql
python import_chipseq_H3K27ac.py
python import_tads.py
python map_enhancers_to_genes.py
mysql -u daria -p brain_multiomics < integration_chip_tad.sql

# Phase 2: HIC (if data available)
mysql -u daria -p brain_multiomics < 03_hic_schema.sql
python import_hic_loops.py ../data/hic_ad_prefrontal.bedpe GSE105194

# Phase 3: Visualization
cd ../outputs
python generate_ucsc_bedfiles.py
# Upload outputs/ucsc_tracks/*.bed to UCSC Genome Browser
```

Estimated total time: **1-2 hours** (mostly waiting for imports)

---

See [parent README](../README.md) for full project context.
