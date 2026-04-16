# Epigenetic Layers Integrated

Alzheimer vs control x brain plasticity genes x expression x ATAC-seq x methylation

## Project overview
This is an end-to-end data integration project that builds a queryable MySQL resource from three AD prefrontal cortex molecular layers:
- Expression (`GSE33000`, microarray)
- Chromatin accessibility (`GSE129040`, ATAC broadPeak)
- DNA methylation (`GSE59685`, Illumina CpG beta values)

The pipeline imports each layer into MySQL, maps regulatory features to genes, and focuses analysis on curated brain plasticity genes to support biologically interpretable AD comparisons.

## Why this project is useful
- Demonstrates practical multi-omics integration across different assay types and identifier systems.
- Shows robust data engineering from GEO metadata parsing to relational modeling.
- Produces reusable SQL analyses for cross-layer questions (expression, promoter accessibility, promoter methylation).

## Portfolio highlights
- Built a reproducible import pipeline in Python + MySQL.
- Added mapping layers (`atac_gene_links`, `methylation_gene_links`) to connect assays at gene level.
- Implemented focused methylation loading to avoid huge unnecessary tables.
- Created integration queries for AD vs control comparisons on plasticity genes.

## Results snapshot
Key outputs from the SQL analysis workflow (`STEP2_mysql_import/insert.sql`):

| Metric | Result |
|---|---:|
| Plasticity genes with RNA layer coverage | 365 |
| Plasticity genes with ATAC layer coverage | 412 |
| Plasticity genes with methylation layer coverage | 398 |
| RNA rows overlapping plasticity genes | 556,818 |
| ATAC peaks linked to plasticity genes | 4,603 |
| Plasticity genes with both AD and control expression data | 365 |
| Plasticity genes hypermethylated at promoters in AD | 141 |
| Genes covered by all three layers | 344 to 356* |

`*` The triple-layer overlap count differs across exploratory queries due to slightly different filters used in those SQL runs.

## Current repository structure

```text
multiomics_data/
├── README.md
├── .gitignore
├── STEP1_download-meta-data/
│   ├── clean_metadata_full.csv
│   ├── clean_metadata_with_internal_ids.csv
│   ├── GSE*_series_matrix.txt
│   ├── GPL4372.annot
│   └── tranpose_rna_parse_metadata_resolveissuename.py
├── STEP2_mysql_import/
│   ├── 01_mysql_schema_update.sql
│   ├── 02a_import_metadata.py
│   ├── 02b_import_rnaseq_UPDATED.py
│   ├── 02c_import_atac.py
│   ├── 02d_works_import_methylation_smart.py
│   ├── 03a_map_atac_to_genes.py
│   ├── 03b_map_methylation_to_genes_skip.py
│   ├── 04a_import_plasticitygenes_geo_into_mysql.py
│   ├── 04_combine_plasticity_gene_lists.py
│   ├── 05_integration_queries.sql
│   └── helper/legacy scripts (see section below)
├── STEP3_additionaldata/   # planned (ChIP-seq / TAD / Hi-C integration)
├── scripts_used/           # metadata preparation helpers
├── Data_files/             # raw large files (ignored in git)
└── metadata/               # downloaded metadata cache (ignored in git)
```

## What each main script does

### Step 1: Metadata preparation
- `scripts_used/import_metadata.py`: downloads GEO `series_matrix.txt.gz` metadata files for selected GSE datasets.
- `scripts_used/tranpose_parse_metadata.py`: parses `!Sample_` fields from series-matrix files and creates a combined sample metadata table.
- `scripts_used/import_metadata33000.py`: enriches metadata by scraping GSM pages (disease, age, sex) for GSE33000-style fields.
- `scripts_used/import_metadata59685.py`: similar enrichment for GSE59685 field naming (including `ad.disease.status`).
- `STEP1_download-meta-data/tranpose_rna_parse_metadata_resolveissuename.py`: finalized parser/importer for GSE33000 series matrix + GPL4372 probe annotation; writes expression into `rna_expression`.

### Step 2: MySQL integration pipeline
- `STEP2_mysql_import/01_mysql_schema_update.sql`: adds/updates schema objects used by integration (new link tables and integrated view).
- `STEP2_mysql_import/02a_import_metadata.py`: imports cleaned metadata into `metadata_import` and `samples`.
- `STEP2_mysql_import/02b_import_rnaseq_UPDATED.py`: imports GSE33000 expression values and maps probes to gene symbols via annotation file.
- `STEP2_mysql_import/02c_import_atac.py`: imports all ATAC broadPeak files into `atac_peaks` and registers sample IDs.
- `STEP2_mysql_import/02d_works_import_methylation_smart.py`: methylation importer with plasticity-focused strategy to avoid huge full-table loads.
- `STEP2_mysql_import/03a_map_atac_to_genes.py`: maps ATAC peaks to nearby genes (promoter/distal/intragenic) into `atac_gene_links`.
- `STEP2_mysql_import/03b_map_methylation_to_genes_skip.py`: maps methylation probes to genes and regions (`TSS200`, `Body`, etc.) into `methylation_gene_links`.
- `STEP2_mysql_import/04a_import_plasticitygenes_geo_into_mysql.py`: imports curated plasticity gene list into `plasticity_genes`.
- `STEP2_mysql_import/04_combine_plasticity_gene_lists.py`: adds category/source annotations for plasticity genes.
- `STEP2_mysql_import/05_integration_queries.sql`: analysis queries for AD vs control across all integrated layers.

### Step 3: Planned extension
- `STEP3_additionaldata/`: reserved for ChIP/TAD/Hi-C integration; not yet part of current main run.

## Run order (current pipeline)

```bash
# 1) Apply schema updates
mysql -u <user> -p brain_multiomics < STEP2_mysql_import/01_mysql_schema_update.sql

# 2) Import sample metadata
python STEP2_mysql_import/02a_import_metadata.py \
  --file STEP1_download-meta-data/clean_metadata_full.csv

# 3) Import expression (GSE33000)
python STEP2_mysql_import/02b_import_rnaseq_UPDATED.py \
  --file STEP1_download-meta-data/GSE33000_series_matrix.txt \
  --annotation STEP2_mysql_import/GPL10558.annot

# 4) Import ATAC peaks (all files in folder)
python STEP2_mysql_import/02c_import_atac.py

# 5) Import methylation (plasticity-aware strategy)
python STEP2_mysql_import/02d_works_import_methylation_smart.py \
  --file Data_files/GSE59685_betas.csv \
  --annotation STEP2_mysql_import/GPL13534_annotation.csv \
  --strategy plasticity

# 6) Map ATAC peaks to genes
python STEP2_mysql_import/03a_map_atac_to_genes.py \
  --annotation STEP2_mysql_import/gene_annotation.tsv

# 7) Map CpGs to genes
python STEP2_mysql_import/03b_map_methylation_to_genes_skip.py \
  --annotation STEP2_mysql_import/GPL13534_annotation.csv

# 8) Load and annotate plasticity genes
python STEP2_mysql_import/04a_import_plasticitygenes_geo_into_mysql.py
python STEP2_mysql_import/04_combine_plasticity_gene_lists.py

# 9) Run integration queries
mysql -u <user> -p brain_multiomics < STEP2_mysql_import/05_integration_queries.sql

# 10) Run additional analysis SQL used for reporting/presentation
# Note: insert.sql currently starts with: USE epigenetic_manifest;
# If your database is brain_multiomics, edit that line first.
mysql -u <user> -p < STEP2_mysql_import/insert.sql
```

## Additional SQL analysis file
- `STEP2_mysql_import/insert.sql`: collection of exploratory and reporting SQL queries used during results interpretation and presentation.
- Includes checks for table coverage, plasticity-gene overlap, AD vs control comparisons, and multi-layer gene-focused queries.
- This file was added from your MySQL working directory (`/home/daria/mysql/databases_epigenetic_manifest`).

## Legacy or duplicate scripts (cleanup candidates)
These look like earlier iterations or alternatives and are likely not in the final run order:
- `STEP2_mysql_import/02b_import_rnaseq_legacy_autoparse.py`
- `STEP2_mysql_import/importrnaseq.py`
- `STEP2_mysql_import/map_ATACcoordinates_genes.py`
- `STEP2_mysql_import/map_ATACcoordinates_genes2.py`
- `STEP2_mysql_import/02c_remap_aliases.py`
- `STEP2_mysql_import/diagnose_methylation.py`
- `STEP2_mysql_import/diagnose_series_matrix.py`

Keep them if you still use them for debugging or one-off recovery. Otherwise, move them into `archive/` or remove.

### What `02b_import_rnaseq_legacy_autoparse.py` does
- This is an older fallback importer for expression files.
- It auto-detects whether the first row is a header.
- If no header exists, it creates `Sample_0`, `Sample_1`, ... column names.
- It performs aggressive filtering of likely invalid gene identifiers (`RSE_*`, `LOC*`, numeric-only values).
- It inserts sample IDs into `samples` and expression values into `rna_expression`.
- It is kept as a legacy backup, not the recommended primary importer.

## Notes for public GitHub
- Raw files in `Data_files/` and `metadata/` are intentionally excluded by `.gitignore`.
- Database passwords are now read from environment variables (`DB_PASSWORD`) instead of hard-coded values.

### Environment setup
Set your DB password in the shell before running import scripts:

```bash
export DB_PASSWORD='your_mysql_password'
```

### Public push checklist
Use this flow for a clean public push:

```bash
git status
git add .
git commit -m "Update portfolio README and SQL analysis scripts"
git push
```

If a password was previously committed, rotate that password in MySQL even after cleanup.
