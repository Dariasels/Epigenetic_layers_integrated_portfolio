-- ══════════════════════════════════════════════════════════════════
-- STEP 1: Update existing empty tables to match what we will import
-- ══════════════════════════════════════════════════════════════════

USE brain_multiomics;

-- Update enhancers table (was too minimal before)
-- Drop and recreate with proper columns for H3K27ac ChIP-seq peaks
DROP TABLE IF EXISTS enhancer_gene_links;
DROP TABLE IF EXISTS enhancers;

CREATE TABLE enhancers (
    enhancer_id     INT AUTO_INCREMENT PRIMARY KEY,
    chrom           VARCHAR(10)   NOT NULL,
    start_pos       INT UNSIGNED  NOT NULL,
    end_pos         INT UNSIGNED  NOT NULL,
    peak_name       VARCHAR(100),
    signal_value    FLOAT,          -- fold enrichment over input
    p_value         FLOAT,          -- -log10 p-value from peak caller
    q_value         FLOAT,          -- -log10 q-value (FDR corrected)
    sample_id       VARCHAR(50),
    cell_type       VARCHAR(100),   -- e.g. NeuN+ (neuronal), Pu.1+ (microglia)
    dataset         VARCHAR(50)     -- GSE102538
);

-- Recreate enhancer_gene_links bridge table
CREATE TABLE enhancer_gene_links (
    link_id         INT AUTO_INCREMENT PRIMARY KEY,
    enhancer_id     INT           NOT NULL,
    gene_symbol     VARCHAR(50)   NOT NULL,
    distance_bp     INT,           -- distance from enhancer midpoint to gene TSS
    region_type     VARCHAR(50),   -- 'promoter_proximal', 'distal_enhancer'
    FOREIGN KEY (enhancer_id)   REFERENCES enhancers(enhancer_id)
                                ON DELETE CASCADE,
    FOREIGN KEY (gene_symbol)   REFERENCES plasticity_genes(gene_symbol)
                                ON DELETE CASCADE
);

-- Update tads table (was too minimal)
DROP TABLE IF EXISTS tads;

CREATE TABLE tads (
    tad_id          INT AUTO_INCREMENT PRIMARY KEY,
    chrom           VARCHAR(10)   NOT NULL,
    start_pos       INT UNSIGNED  NOT NULL,
    end_pos         INT UNSIGNED  NOT NULL,
    tad_size_kb     INT,           -- size in kilobases (for quick filtering)
    source          VARCHAR(100)   DEFAULT 'GSE105194 ENCODE hg19',
    INDEX idx_tad_chrom (chrom),
    INDEX idx_tad_pos   (chrom, start_pos, end_pos)
);

-- New bridge table: which TAD does each plasticity gene fall in?
CREATE TABLE IF NOT EXISTS tad_gene_links (
    link_id         INT AUTO_INCREMENT PRIMARY KEY,
    tad_id          INT           NOT NULL,
    gene_symbol     VARCHAR(50)   NOT NULL,
    FOREIGN KEY (tad_id)        REFERENCES tads(tad_id)
                                ON DELETE CASCADE
    -- No FK to plasticity_genes here — we map ALL genes, not just plasticity
    -- so we can later ask which TAD contains BOTH a plasticity gene and an enhancer
);

-- New bridge table: which TAD does each ATAC peak fall in?
-- This is the key integration query: same TAD = same regulatory neighbourhood
CREATE TABLE IF NOT EXISTS tad_atac_links (
    link_id         INT AUTO_INCREMENT PRIMARY KEY,
    tad_id          INT           NOT NULL,
    peak_id         INT           NOT NULL,
    FOREIGN KEY (tad_id)  REFERENCES tads(tad_id)     ON DELETE CASCADE,
    FOREIGN KEY (peak_id) REFERENCES atac_peaks(peak_id) ON DELETE CASCADE
);

-- ══════════════════════════════════════════════════════════════════
-- STEP 2: Add GSM3692183 metadata to samples table
-- GSE102538: H3K27ac ChIP-seq, prefrontal cortex, AD vs control
-- Cell types: NeuN+ (neurons), Pu.1+ (microglia), NeuN-/Pu.1- (OEG)
-- ══════════════════════════════════════════════════════════════════

INSERT IGNORE INTO samples (sample_id, `condition`, tissue, dataset)
VALUES
  ('GSM3692183', 'Alzheimer', 'prefrontal cortex', 'GSE102538'),
  ('GSM3692184', 'Alzheimer', 'prefrontal cortex', 'GSE102538'),
  ('GSM3692185', 'Alzheimer', 'prefrontal cortex', 'GSE102538'),
  ('GSM3692186', 'Control',   'prefrontal cortex', 'GSE102538'),
  ('GSM3692187', 'Control',   'prefrontal cortex', 'GSE102538'),
  ('GSM3692188', 'Control',   'prefrontal cortex', 'GSE102538');
-- Note: add all GSM IDs from the series — these are examples.
-- Run the import script which will auto-insert them from the broadPeak files.

-- Verify
SELECT dataset, `condition`, COUNT(*) FROM samples
WHERE dataset = 'GSE102538'
GROUP BY dataset, `condition`;
