-- ============================================================
-- BRAIN MULTI-OMICS DATABASE — SCHEMA UPDATE
-- Run this AFTER your original CREATE TABLE statements
-- It adds missing columns and new helper tables
-- ============================================================

USE brain_multiomics;

-- ------------------------------------------------------------
-- 1. ADD sample_id TO atac_peaks
--    (Your original table was missing this — we need to know
--     which sample each peak came from)
-- ------------------------------------------------------------
ALTER TABLE atac_peaks
    ADD COLUMN sample_id VARCHAR(50) AFTER peak_id,
    ADD COLUMN peak_name VARCHAR(100) AFTER chrom_end,
    ADD FOREIGN KEY fk_atac_sample (sample_id) REFERENCES samples(sample_id);


-- ------------------------------------------------------------
-- 2. NEW TABLE: plasticity_genes
--    Stores the curated list of brain plasticity genes
--    that we will use to filter all other tables
-- ------------------------------------------------------------
CREATE TABLE IF NOT EXISTS plasticity_genes (
    gene_symbol  VARCHAR(50) PRIMARY KEY,
    category     VARCHAR(100),   -- e.g. "synaptic plasticity", "LTP", "neurogenesis"
    source       VARCHAR(100)    -- where the gene came from (e.g. "SynGO", "manual curation")
);


-- ------------------------------------------------------------
-- 3. NEW TABLE: atac_gene_links
--    After mapping ATAC peaks to nearby genes (script 03a),
--    we store those links here
-- ------------------------------------------------------------
CREATE TABLE IF NOT EXISTS atac_gene_links (
    link_id      INT AUTO_INCREMENT PRIMARY KEY,
    peak_id      INT,
    gene_symbol  VARCHAR(50),
    distance_bp  INT,            -- how far the peak is from the gene TSS
    region_type  VARCHAR(50),    -- "promoter", "distal", "intragenic"
    FOREIGN KEY (peak_id) REFERENCES atac_peaks(peak_id)
);


-- ------------------------------------------------------------
-- 4. NEW TABLE: methylation_gene_links
--    After mapping CpG probes to genes (script 03b),
--    we store those links here
-- ------------------------------------------------------------
CREATE TABLE IF NOT EXISTS methylation_gene_links (
    link_id      INT AUTO_INCREMENT PRIMARY KEY,
    cpg_id       VARCHAR(50),    -- e.g. cg00004067
    gene_symbol  VARCHAR(50),
    relation     VARCHAR(50)     -- "TSS200", "TSS1500", "Body", "1stExon" etc.
                                 -- this comes from Illumina EPIC annotation
);


-- ------------------------------------------------------------
-- 5. NEW VIEW: integrated_plasticity
--    Once all data is imported and mapped, this view lets you
--    query all three layers together for plasticity genes.
--    Run this AFTER all import and mapping scripts are done.
-- ------------------------------------------------------------
CREATE OR REPLACE VIEW integrated_plasticity AS
SELECT
    pg.gene_symbol,
    pg.category                         AS plasticity_category,

    -- RNA expression
    re.sample_id                        AS rna_sample_id,
    s_rna.condition                     AS rna_condition,
    re.expression_value,

    -- Methylation
    m.sample_id                         AS meth_sample_id,
    s_meth.condition                    AS meth_condition,
    m.cpg_id,
    mgl.relation                        AS cpg_region,
    m.beta_value,

    -- ATAC chromatin accessibility
    ap.sample_id                        AS atac_sample_id,
    s_atac.condition                    AS atac_condition,
    ap.chrom,
    ap.chrom_start,
    ap.chrom_end,
    ap.signal_value                     AS atac_signal,
    agl.region_type                     AS atac_region_type

FROM plasticity_genes pg

-- Join RNA expression via gene symbol
LEFT JOIN rna_expression re
    ON re.gene_symbol = pg.gene_symbol
LEFT JOIN samples s_rna
    ON s_rna.sample_id = re.sample_id

-- Join methylation via cpg→gene link table
LEFT JOIN methylation_gene_links mgl
    ON mgl.gene_symbol = pg.gene_symbol
LEFT JOIN methylation m
    ON m.cpg_id = mgl.cpg_id
LEFT JOIN samples s_meth
    ON s_meth.sample_id = m.sample_id

-- Join ATAC peaks via peak→gene link table
LEFT JOIN atac_gene_links agl
    ON agl.gene_symbol = pg.gene_symbol
LEFT JOIN atac_peaks ap
    ON ap.peak_id = agl.peak_id
LEFT JOIN samples s_atac
    ON s_atac.sample_id = ap.sample_id;


-- ------------------------------------------------------------
-- Confirm everything was created
-- ------------------------------------------------------------
SHOW TABLES;
DESCRIBE atac_peaks;
DESCRIBE plasticity_genes;
DESCRIBE atac_gene_links;
DESCRIBE methylation_gene_links;
