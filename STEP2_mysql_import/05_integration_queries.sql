-- ============================================================
-- 05_integration_queries.sql
-- Useful queries to explore the integrated multi-omics data
-- Run these in MySQL after all import and mapping scripts are done
-- ============================================================

USE brain_multiomics;

-- ─────────────────────────────────────────────────────────────
-- Q1. How many plasticity genes do we have data for?
-- ─────────────────────────────────────────────────────────────
SELECT
    pg.category,
    COUNT(DISTINCT pg.gene_symbol)                      AS total_genes_in_category,
    COUNT(DISTINCT re.gene_symbol)                      AS genes_with_RNA,
    COUNT(DISTINCT agl.gene_symbol)                     AS genes_with_ATAC,
    COUNT(DISTINCT mgl.gene_symbol)                     AS genes_with_methylation
FROM plasticity_genes pg
LEFT JOIN rna_expression re        ON re.gene_symbol  = pg.gene_symbol
LEFT JOIN atac_gene_links agl      ON agl.gene_symbol = pg.gene_symbol
LEFT JOIN methylation_gene_links mgl ON mgl.gene_symbol = pg.gene_symbol
GROUP BY pg.category
ORDER BY total_genes_in_category DESC;


-- ─────────────────────────────────────────────────────────────
-- Q2. Average expression of plasticity genes: AD vs control
-- ─────────────────────────────────────────────────────────────
SELECT
    re.gene_symbol,
    pg.category,
    s.condition,
    COUNT(DISTINCT re.sample_id)       AS n_samples,
    AVG(re.expression_value)           AS mean_expression,
    MIN(re.expression_value)           AS min_expression,
    MAX(re.expression_value)           AS max_expression
FROM rna_expression re
INNER JOIN plasticity_genes pg ON pg.gene_symbol = re.gene_symbol
INNER JOIN samples s            ON s.sample_id   = re.sample_id
WHERE s.condition IN ('AD', 'control')
GROUP BY re.gene_symbol, pg.category, s.condition
ORDER BY re.gene_symbol, s.condition;


-- ─────────────────────────────────────────────────────────────
-- Q3. BDNF specifically — expression AD vs control
-- ─────────────────────────────────────────────────────────────
SELECT
    s.condition,
    s.dataset,
    COUNT(DISTINCT re.sample_id)    AS n_samples,
    AVG(re.expression_value)        AS mean_expression
FROM rna_expression re
INNER JOIN samples s ON s.sample_id = re.sample_id
WHERE re.gene_symbol = 'BDNF'
  AND s.condition IN ('AD', 'control')
GROUP BY s.condition, s.dataset;


-- ─────────────────────────────────────────────────────────────
-- Q4. Average methylation (beta value) of plasticity gene
--     promoters: AD vs control
-- (TSS200 = within 200bp of TSS — the "core promoter")
-- ─────────────────────────────────────────────────────────────
SELECT
    mgl.gene_symbol,
    pg.category,
    s.condition,
    mgl.relation                       AS cpg_region,
    COUNT(DISTINCT m.cpg_id)           AS n_cpg_probes,
    COUNT(DISTINCT m.sample_id)        AS n_samples,
    AVG(m.beta_value)                  AS mean_beta,
    STDDEV(m.beta_value)               AS sd_beta
FROM methylation m
INNER JOIN methylation_gene_links mgl  ON mgl.cpg_id     = m.cpg_id
INNER JOIN plasticity_genes pg         ON pg.gene_symbol  = mgl.gene_symbol
INNER JOIN samples s                   ON s.sample_id     = m.sample_id
WHERE mgl.relation IN ('TSS200', 'TSS1500', '1stExon')   -- promoter regions
  AND s.condition IN ('AD', 'control')
GROUP BY mgl.gene_symbol, pg.category, s.condition, mgl.relation
ORDER BY mgl.gene_symbol, mgl.relation, s.condition;


-- ─────────────────────────────────────────────────────────────
-- Q5. ATAC-seq: chromatin accessibility near plasticity genes
--     How many peaks are in promoters of plasticity genes?
-- ─────────────────────────────────────────────────────────────
SELECT
    agl.gene_symbol,
    pg.category,
    agl.region_type,
    s.condition,
    COUNT(DISTINCT ap.peak_id)         AS n_peaks,
    COUNT(DISTINCT ap.sample_id)       AS n_samples,
    AVG(ap.signal_value)               AS mean_signal
FROM atac_peaks ap
INNER JOIN atac_gene_links agl  ON agl.peak_id    = ap.peak_id
INNER JOIN plasticity_genes pg  ON pg.gene_symbol = agl.gene_symbol
INNER JOIN samples s            ON s.sample_id    = ap.sample_id
WHERE s.condition IN ('AD', 'control')
GROUP BY agl.gene_symbol, pg.category, agl.region_type, s.condition
ORDER BY agl.gene_symbol, agl.region_type, s.condition;


-- ─────────────────────────────────────────────────────────────
-- Q6. FULL INTEGRATION VIEW — one gene, all three layers
--     Example: BDNF across all samples
-- This uses the integrated_plasticity view created in the schema
-- ─────────────────────────────────────────────────────────────
SELECT
    gene_symbol,
    plasticity_category,

    -- RNA
    rna_condition,
    AVG(expression_value)              AS mean_expression,

    -- Methylation (promoter only)
    meth_condition,
    cpg_region,
    AVG(beta_value)                    AS mean_beta,

    -- ATAC
    atac_condition,
    atac_region_type,
    AVG(atac_signal)                   AS mean_atac_signal

FROM integrated_plasticity
WHERE gene_symbol = 'BDNF'
GROUP BY
    gene_symbol, plasticity_category,
    rna_condition, meth_condition, cpg_region,
    atac_condition, atac_region_type;


-- ─────────────────────────────────────────────────────────────
-- Q7. Find genes with BOTH high methylation AND low expression
--     in AD (potential epigenetic silencing candidates)
-- ─────────────────────────────────────────────────────────────
SELECT
    ad_meth.gene_symbol,
    pg.category,
    ad_meth.mean_beta_AD,
    ctrl_meth.mean_beta_ctrl,
    (ad_meth.mean_beta_AD - ctrl_meth.mean_beta_ctrl)     AS delta_beta,
    ad_rna.mean_expr_AD,
    ctrl_rna.mean_expr_ctrl

FROM plasticity_genes pg

-- Methylation in AD
INNER JOIN (
    SELECT mgl.gene_symbol, AVG(m.beta_value) AS mean_beta_AD
    FROM methylation m
    JOIN methylation_gene_links mgl ON mgl.cpg_id = m.cpg_id
    JOIN samples s ON s.sample_id = m.sample_id
    WHERE s.condition = 'AD' AND mgl.relation IN ('TSS200','TSS1500')
    GROUP BY mgl.gene_symbol
) ad_meth ON ad_meth.gene_symbol = pg.gene_symbol

-- Methylation in control
INNER JOIN (
    SELECT mgl.gene_symbol, AVG(m.beta_value) AS mean_beta_ctrl
    FROM methylation m
    JOIN methylation_gene_links mgl ON mgl.cpg_id = m.cpg_id
    JOIN samples s ON s.sample_id = m.sample_id
    WHERE s.condition = 'control' AND mgl.relation IN ('TSS200','TSS1500')
    GROUP BY mgl.gene_symbol
) ctrl_meth ON ctrl_meth.gene_symbol = pg.gene_symbol

-- Expression in AD
INNER JOIN (
    SELECT re.gene_symbol, AVG(re.expression_value) AS mean_expr_AD
    FROM rna_expression re
    JOIN samples s ON s.sample_id = re.sample_id
    WHERE s.condition = 'AD'
    GROUP BY re.gene_symbol
) ad_rna ON ad_rna.gene_symbol = pg.gene_symbol

-- Expression in control
INNER JOIN (
    SELECT re.gene_symbol, AVG(re.expression_value) AS mean_expr_ctrl
    FROM rna_expression re
    JOIN samples s ON s.sample_id = re.sample_id
    WHERE s.condition = 'control'
    GROUP BY re.gene_symbol
) ctrl_rna ON ctrl_rna.gene_symbol = pg.gene_symbol

-- Filter: hypermethylated (delta_beta > 0.1) AND expressed lower in AD
HAVING delta_beta > 0.1
   AND mean_expr_AD < mean_expr_ctrl

ORDER BY delta_beta DESC;


-- ─────────────────────────────────────────────────────────────
-- Q8. Open chromatin (ATAC) + low methylation in AD promoters
--     These may be genes that are ACTIVE / upregulated in AD
-- ─────────────────────────────────────────────────────────────
SELECT
    pg.gene_symbol,
    pg.category,
    AVG(m.beta_value)       AS mean_beta_AD,
    AVG(ap.signal_value)    AS mean_atac_signal_AD
FROM plasticity_genes pg
INNER JOIN methylation_gene_links mgl ON mgl.gene_symbol = pg.gene_symbol
INNER JOIN methylation m              ON m.cpg_id = mgl.cpg_id
INNER JOIN samples sm                 ON sm.sample_id = m.sample_id
INNER JOIN atac_gene_links agl        ON agl.gene_symbol = pg.gene_symbol
INNER JOIN atac_peaks ap              ON ap.peak_id = agl.peak_id
INNER JOIN samples sa                 ON sa.sample_id = ap.sample_id
WHERE sm.condition = 'AD'
  AND sa.condition = 'AD'
  AND mgl.relation IN ('TSS200', 'TSS1500')
  AND agl.region_type = 'promoter'
GROUP BY pg.gene_symbol, pg.category
HAVING mean_beta_AD < 0.3           -- low methylation = open for transcription
   AND mean_atac_signal_AD > 5      -- high ATAC signal = accessible chromatin
ORDER BY mean_atac_signal_AD DESC;
