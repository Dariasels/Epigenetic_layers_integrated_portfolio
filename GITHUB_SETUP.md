# GitHub Setup & Contribution Guide

## Initial GitHub Setup

### 1. Create Remote Repository

```bash
# Initialize git in your local directory
cd /user/gent/501/vsc50116/Epigenetic_layers_integrated_portfolio
git init
git config user.name "Your Name"
git config user.email "your.email@example.com"
```

### 2. Create Repository on GitHub

Steps:
1. Go to [GitHub.com](https://github.com)
2. Log in or create an account
3. Click **+** (top right) → **New repository**
4. Name: `epigenetic-layers-integrated-portfolio`
5. Description: "Multi-omics integration of brain plasticity genes in Alzheimer's Disease (RNA-seq + ATAC-seq + Methylation + ChIP-seq + HIC)"
6. Choose **Private** (add public later if desired)
7. DO NOT initialize with README (we have one)
8. Click **Create repository**

### 3. Add Remote & Initial Push

```bash
# Add remote
git remote add origin https://github.com/YOUR_USERNAME/epigenetic-layers-integrated-portfolio.git

# Verify
git remote -v

# Initial commit
git add .
git commit -m "Initial commit: Multi-layer epigenetic integration pipeline

- STEP 1: Metadata preparation from GEO
  - GSE33000 (RNA-seq, Agilent 44K microarray)
  - GSE129040 (ATAC-seq open chromatin)
  - GSE59685 (Methylation 450k array)
  
- STEP 2: MySQL integration (complete)
  - Core 3-layer database schema
  - Import scripts for all data types
  - Gene mapping (ATAC→genes, Meth→genes)
  - Plasticity gene curation (365+ genes)
  - Multi-layer integration views
  
- STEP 3: Extended integration (ready)
  - H3K27ac ChIP-seq enhancers (GSE102538)
  - TAD integration scripts
  - HIC contact mapping (3D chromatin)
  - UCSC Genome Browser visualization
  
Coverage: 365+ brain plasticity genes across all layers"

# Push to GitHub
git branch -M main
git push -u origin main
```

---

## Ongoing Workflow

### Before Making Changes

```bash
# Update local repository
git pull origin main
```

### Making Changes

```bash
# Create feature branch (for major work)
git checkout -b feature/hic-integration
# or: git checkout -b fix/schema-bug

# Make changes...
git add .
git commit -m "Add HIC import functionality

- Implement import_hic_loops.py with BEDPE support
- Map HIC endpoints to genes within 50kb window
- Link HIC loops to H3K27ac enhancers
- Generate summary statistics
- Add comprehensive error handling"
```

### Pushing Changes

```bash
# Push feature branch
git push origin feature/hic-integration

# Create Pull Request on GitHub (optional, for review)
# Then merge to main when ready:
git checkout main
git pull origin main
git merge feature/hic-integration
git push origin main
```

---

## Repository Structure (for GitHub)

```
epigenetic-layers-integrated-portfolio/
├── README.md                           # Project overview
├── .gitignore                          # Ignore databases, credentials, large files
├── LICENSE                             # Recommend: MIT or GPL v3
│
├── STEP1_download-meta-data/           # Metadata preparation
│   ├── README.md
│   ├── clean_metadata_full.csv
│   └── tranpose_rna_parse_metadata_resolveissuename.py
│
├── STEP2_mysql_import/                 # Core 3-layer integration
│   ├── README.md
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
│   ├── annotation_import_methylation.r
│   ├── gene_annotation.tsv
│   ├── GPL13534_annotation.csv
│   ├── plasticity_genes.csv
│   └── legacy/
│
├── STEP3_additionaldata/               # Extended integration (ChIP + TAD + HIC)
│   ├── README.md
│   ├── 02_chip_tad_schema.sql
│   ├── 03_hic_schema.sql              # NEW
│   ├── import_chipseq_H3K27ac.py
│   ├── import_tads.py
│   ├── import_hic_loops.py            # NEW
│   ├── map_enhancers_to_genes.py
│   ├── map_atac_to_genes.py
│   └── integration_chip_tad.sql
│
├── outputs/                            # (Not in git — large files)
│   └── generate_ucsc_bedfiles.py      # NEW
│
├── scripts_used/                       # Historical scripts
├── errorhandling.py
└── .gitignore
```

---

## GitHub Best Practices

### Commit Messages

Use clear, descriptive commit messages:

```
Good:
  git commit -m "Add HIC integration with BEDPE format support

  - Parse BEDPE files with contact strength normalization
  - Map HIC loop endpoints to genes (±50kb window)
  - Link HIC loops to H3K27ac enhancers
  - Add summary statistics query
  - Handles edge cases: chromosome normalization, duplicates"

Bad:
  git commit -m "Update script"
  git commit -m "Fix bugs"
```

### Branching Strategy

Use feature branches for clarity:

```bash
# Feature: new functionality
git checkout -b feature/hic-integration

# Fix: bug fix
git checkout -b fix/chromosome-normalization

# Doc: documentation updates
git checkout -b doc/add-troubleshooting-guide
```

### Keeping Repository Clean

```bash
# Remove sensitive data (if accidentally added)
git rm --cached .env
echo ".env" >> .gitignore
git commit -m "Remove .env file (credentials)"

# Delete old branches
git branch -d feature/old-feature
git push origin --delete feature/old-feature
```

---

## Sharing & Collaboration

### Making Repository Public

1. Go to GitHub repository settings
2. Change **Privacy** to **Public**
3. Add LICENSE file

### Adding Collaborators

1. Settings → **Collaborators**
2. Search for username
3. Grant access (write/admin)

### Using Pull Requests (for reviews)

```bash
# Team member forks repo
# Makes changes in feature branch
# Submits Pull Request on GitHub
# You review, request changes, merge when approved
```

---

## Backup & Safety

```bash
# Create local backup
tar -czf epigenetic-portfolio-backup-$(date +%Y%m%d).tar.gz .git

# Clone to second location
git clone https://github.com/YOUR_USERNAME/epigenetic-layers-integrated-portfolio.git backup_copy

# Regular pushes ensure GitHub is backup copy
git push origin main  # Always push important changes
```

---

## Useful GitHub Commands

```bash
# View commit history
git log --oneline -n 20

# See changes since last commit
git diff

# View remote branches
git branch -r

# Fetch latest without merging
git fetch origin

# Squash commits (before push)
git rebase -i HEAD~3

# Tag versions
git tag -a v1.0 -m "Initial 3-layer integration complete"
git push origin v1.0
```

---

## Publishing to GitHub for Data Sharing

Once integration is complete:

1. **README.md** — Already complete, explains full pipeline
2. **LICENSE** — Add MIT or GPL v3 (science-friendly)
3. **Issues** — Enable issue tracking for questions/bugs
4. **Releases** — Create releases for each milestone (v1.0, v1.1, etc.)
5. **Citation** — Add CITATION.cff file for academic citation

### Example CITATION.cff

```yaml
cff-version: 1.2.0
title: "Multi-omics Integration of Brain Plasticity Genes in Alzheimer's Disease"
authors:
  - given-names: Your
    family-names: Name
    affiliation: "Your Institution"
date-released: 2026-04-27
version: 1.0.0
repository-code: "https://github.com/YOUR_USERNAME/epigenetic-layers-integrated-portfolio"
license: MIT
```

---

## Questions?

GitHub documentation: https://docs.github.com
Git tutorial: https://git-scm.com/doc

Good luck with the project! 🚀
