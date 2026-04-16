# Brain Multi-Omics Integration Pipeline
## Epigenetic layers in Alzheimer's Disease — Prefrontal Cortex

---

## What this project does

You are integrating **three epigenetic/transcriptomic datasets** from Alzheimer's Disease
prefrontal cortex samples, then filtering to genes related to **brain plasticity**.

| Dataset   | Type              | File format      | Notes                                      |
|-----------|-------------------|------------------|--------------------------------------------|
| GSE33000  | Gene Expression   | Microarray (TSV) | Already normalized — values are log2-like  |
| GSE129040 | Chromatin Access  | broadPeak        | ATAC-seq peaks, need gene mapping          |
| GSE59685  | DNA Methylation   | Beta values CSV  | CpG probe IDs, need mapping via EPIC array |

> **Important:** GSE33000 is a **microarray** dataset, not RNA-seq.
> The data is already processed/normalized. You do NOT need to run alignment.
> The `Gene` column in your file already has gene symbols — you just need to
> clean and import them.

---

## Folder structure

```
brain_multiomics/
├── 01_setup/
│   └── 01_mysql_schema_update.sql      # Updated database schema
├── 02_import/
│   ├── 02a_import_metadata.py          # Import clean_metadata_full.csv
│   ├── 02b_import_rnaseq.py            # Import GSE33000 expression values
│   ├── 02c_import_atac.py              # Import broadPeak ATAC files
│   └── 02d_import_methylation.py       # Import GSE59685 beta values
├── 03_mapping/
│   ├── 03a_map_atac_to_genes.py        # Map ATAC peaks → genes (window-based)
│   └── 03b_map_methylation_to_genes.py # Map CpG probes → genes (annotation file)
├── 04_plasticity/
│   └── 04_filter_plasticity_genes.py   # Filter all tables to plasticity genes
├── 05_integration/
│   └── 05_integration_queries.sql      # SQL queries to join all layers
└── README.md                           # This file
```

---

## Step-by-step instructions

### Prerequisites
- Python 3.8+
- MySQL running with database `brain_multiomics` already created
- Python packages: `pip install mysql-connector-python pandas`
- R (for methylation probe annotation): `install.packages("BiocManager"); BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")`

### Run order

```bash
# 1. Update the database schema
mysql -u root -p brain_multiomics < 01_setup/01_mysql_schema_update.sql

# 2. Import metadata first (other tables depend on it)
python 02_import/02a_import_metadata.py

# 3. Import expression data
python 02_import/02b_import_rnaseq.py

# 4. Import ATAC peaks (run for each broadPeak file)
python 02_import/02c_import_atac.py --file path/to/yourfile.broadPeak --sample_id SAMPLE_ID

# 5. Import methylation beta values
python 02_import/02d_import_methylation.py

# 6. Map ATAC peaks to genes (requires gene annotation)
python 03_mapping/03a_map_atac_to_genes.py

# 7. Map methylation probes to genes
python 03_mapping/03b_map_methylation_to_genes.py

# 8. Filter everything to plasticity genes
python 04_plasticity/04_filter_plasticity_genes.py

# 9. Explore integrated data
mysql -u root -p brain_multiomics < 05_integration/05_integration_queries.sql
```

---

## Key concepts (beginner-friendly)

**Why do we need probe mapping?**
- ATAC peaks are genomic coordinates (chr:start-end). We need to know which *gene*
  is nearby (within ~2000 bp = "promoter window").
- Methylation CpG probes (cg00004067 etc.) are on the Illumina EPIC array.
  A lookup table tells us which gene each probe is in.
- RNA expression already has gene symbols in the `Gene` column — easy!

**What are plasticity genes?**
Genes involved in synaptic plasticity, LTP/LTD, neurogenesis, and memory.
Examples: BDNF, ARC, NRXN1, SHANK3, CAMK2A, HOMER1, etc.
We use a curated list from published brain plasticity gene databases.

**What does "integration" mean here?**
Joining the three tables by gene symbol so you can ask:
"For gene BDNF, what is its expression level, is the promoter open (ATAC peak nearby),
and is its CpG island methylated — and does this differ between AD and controls?"
