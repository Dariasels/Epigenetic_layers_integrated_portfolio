-- ══════════════════════════════════════════════════════════════════
-- integration_chip_tad.sql
-- Queries integrating H3K27ac ChIP-seq and TAD data with existing layers
-- ══════════════════════════════════════════════════════════════════

USE brain_multiomics;

-- ─────────────────────────────────────────────────────────────────
-- Q1. Basic coverage check after import
-- ─────────────────────────────────────────────────────────────────

-- How many H3K27ac peaks per condition and cell type?
SELECT s.`condition`, e.cell_type,
       COUNT(*) AS n_peaks,
       AVG(e.signal_value) AS mean_signal
FROM enhancers e
JOIN samples s ON s.sample_id = e.sample_id
GROUP BY s.`condition`, e.cell_type
ORDER BY s.`condition`, e.cell_type;

-- How many TADs were imported?
SELECT COUNT(*) AS n_tads,
       AVG(tad_size_kb) AS avg_size_kb,
       MIN(tad_size_kb) AS min_kb,
       MAX(tad_size_kb) AS max_kb
FROM tads;

-- How many plasticity genes are inside a TAD?
SELECT COUNT(DISTINCT tgl.gene_symbol) AS plasticity_genes_in_tads
FROM tad_gene_links tgl
JOIN plasticity_genes pg ON pg.gene_symbol = tgl.gene_symbol;

-- ─────────────────────────────────────────────────────────────────
-- Q2. H3K27ac signal at plasticity gene promoters: AD vs Control
-- The key question: are AD plasticity gene promoters less acetylated?
-- (lower H3K27ac = less active promoter = less transcription)
-- ─────────────────────────────────────────────────────────────────
SELECT
    pg.gene_symbol,
    pg.category,
    AVG(CASE WHEN s.`condition` = 'Alzheimer' THEN e.signal_value END) AS h3k27ac_AD,
    AVG(CASE WHEN s.`condition` = 'Control'   THEN e.signal_value END) AS h3k27ac_Ctrl,
    AVG(CASE WHEN s.`condition` = 'Alzheimer' THEN e.signal_value END) -
    AVG(CASE WHEN s.`condition` = 'Control'   THEN e.signal_value END) AS h3k27ac_delta,
    COUNT(DISTINCT e.enhancer_id) AS n_peaks
FROM plasticity_genes pg
JOIN enhancer_gene_links egl ON egl.gene_symbol = pg.gene_symbol
JOIN enhancers e             ON e.enhancer_id   = egl.enhancer_id
JOIN samples s               ON s.sample_id     = e.sample_id
WHERE egl.region_type = 'promoter_proximal'
GROUP BY pg.gene_symbol, pg.category
HAVING h3k27ac_AD IS NOT NULL AND h3k27ac_Ctrl IS NOT NULL
ORDER BY h3k27ac_delta ASC   -- most depleted in AD first
LIMIT 30;

-- ─────────────────────────────────────────────────────────────────
-- Q3. The TAD co-occupancy question
-- Which TADs contain BOTH a plasticity gene AND an H3K27ac peak?
-- These are the most likely functional regulatory interactions.
-- ─────────────────────────────────────────────────────────────────
SELECT
    t.chrom,
    t.start_pos,
    t.end_pos,
    t.tad_size_kb,
    COUNT(DISTINCT tgl.gene_symbol) AS plasticity_genes_in_tad,
    COUNT(DISTINCT ap.peak_id)      AS atac_peaks_in_tad,
    COUNT(DISTINCT e.enhancer_id)   AS h3k27ac_peaks_in_tad,
    GROUP_CONCAT(DISTINCT tgl.gene_symbol ORDER BY tgl.gene_symbol) AS gene_list
FROM tads t
JOIN tad_gene_links tgl   ON tgl.tad_id  = t.tad_id
JOIN plasticity_genes pg  ON pg.gene_symbol = tgl.gene_symbol
LEFT JOIN tad_atac_links tal ON tal.tad_id = t.tad_id
LEFT JOIN atac_peaks ap   ON ap.peak_id   = tal.peak_id
LEFT JOIN enhancer_gene_links egl ON egl.gene_symbol = tgl.gene_symbol
LEFT JOIN enhancers e     ON e.enhancer_id = egl.enhancer_id
GROUP BY t.tad_id, t.chrom, t.start_pos, t.end_pos, t.tad_size_kb
HAVING plasticity_genes_in_tad > 0
   AND (atac_peaks_in_tad > 0 OR h3k27ac_peaks_in_tad > 0)
ORDER BY plasticity_genes_in_tad DESC, h3k27ac_peaks_in_tad DESC
LIMIT 20;

-- ─────────────────────────────────────────────────────────────────
-- Q4. FULL 5-LAYER INTEGRATION for a single gene (e.g. BDNF)
-- RNA expression + Methylation + ATAC + H3K27ac + TAD membership
-- ─────────────────────────────────────────────────────────────────
SELECT
    pg.gene_symbol,
    pg.category,

    -- RNA
    ps.rna_AD,
    ps.rna_Ctrl,
    ps.rna_delta,

    -- Methylation (promoter)
    ps.meth_AD,
    ps.meth_Ctrl,
    ps.meth_delta,

    -- ATAC (chromatin accessibility)
    ps.atac_AD,
    ps.atac_Ctrl,
    ps.atac_delta,

    -- H3K27ac (active enhancer/promoter mark)
    (SELECT AVG(e.signal_value)
     FROM enhancers e
     JOIN enhancer_gene_links egl ON e.enhancer_id = egl.enhancer_id
     JOIN samples s ON e.sample_id = s.sample_id
     WHERE egl.gene_symbol = pg.gene_symbol
       AND egl.region_type = 'promoter_proximal'
       AND s.`condition` = 'Alzheimer') AS h3k27ac_AD,

    (SELECT AVG(e.signal_value)
     FROM enhancers e
     JOIN enhancer_gene_links egl ON e.enhancer_id = egl.enhancer_id
     JOIN samples s ON e.sample_id = s.sample_id
     WHERE egl.gene_symbol = pg.gene_symbol
       AND egl.region_type = 'promoter_proximal'
       AND s.`condition` = 'Control') AS h3k27ac_Ctrl,

    -- TAD: which TAD does this gene live in?
    (SELECT CONCAT(t.chrom, ':', t.start_pos, '-', t.end_pos,
                   ' (', t.tad_size_kb, 'kb)')
     FROM tad_gene_links tgl
     JOIN tads t ON t.tad_id = tgl.tad_id
     WHERE tgl.gene_symbol = pg.gene_symbol
     LIMIT 1) AS tad_location,

    -- How many other plasticity genes share this TAD?
    (SELECT COUNT(DISTINCT tgl2.gene_symbol) - 1
     FROM tad_gene_links tgl
     JOIN tad_gene_links tgl2 ON tgl2.tad_id = tgl.tad_id
     JOIN plasticity_genes pg2 ON pg2.gene_symbol = tgl2.gene_symbol
     WHERE tgl.gene_symbol = pg.gene_symbol) AS plasticity_neighbours_in_tad

FROM plasticity_genes pg
JOIN plasticity_summary ps ON ps.gene_symbol = pg.gene_symbol
WHERE pg.gene_symbol = 'BDNF';   -- change this to any gene

-- ─────────────────────────────────────────────────────────────────
-- Q5. Add H3K27ac delta to plasticity_summary table
-- Run this after the H3K27ac import to enrich the summary table
-- ─────────────────────────────────────────────────────────────────

-- First add columns if they don't exist
ALTER TABLE plasticity_summary
  ADD COLUMN IF NOT EXISTS h3k27ac_AD    FLOAT,
  ADD COLUMN IF NOT EXISTS h3k27ac_Ctrl  FLOAT,
  ADD COLUMN IF NOT EXISTS h3k27ac_delta FLOAT,
  ADD COLUMN IF NOT EXISTS tad_id        INT;

-- Then populate H3K27ac
UPDATE plasticity_summary ps
SET
    h3k27ac_AD = (
        SELECT AVG(e.signal_value)
        FROM enhancers e
        JOIN enhancer_gene_links egl ON e.enhancer_id = egl.enhancer_id
        JOIN samples s ON e.sample_id = s.sample_id
        WHERE egl.gene_symbol = ps.gene_symbol
          AND egl.region_type = 'promoter_proximal'
          AND s.`condition` = 'Alzheimer'
    ),
    h3k27ac_Ctrl = (
        SELECT AVG(e.signal_value)
        FROM enhancers e
        JOIN enhancer_gene_links egl ON e.enhancer_id = egl.enhancer_id
        JOIN samples s ON e.sample_id = s.sample_id
        WHERE egl.gene_symbol = ps.gene_symbol
          AND egl.region_type = 'promoter_proximal'
          AND s.`condition` = 'Control'
    );

UPDATE plasticity_summary
SET h3k27ac_delta = h3k27ac_AD - h3k27ac_Ctrl
WHERE h3k27ac_AD IS NOT NULL AND h3k27ac_Ctrl IS NOT NULL;

-- Populate TAD membership
UPDATE plasticity_summary ps
SET tad_id = (
    SELECT tgl.tad_id
    FROM tad_gene_links tgl
    WHERE tgl.gene_symbol = ps.gene_symbol
    LIMIT 1
);

-- Final 5-layer summary check
SELECT
    gene_symbol, category,
    ROUND(rna_delta, 3)      AS rna_delta,
    ROUND(meth_delta, 4)     AS meth_delta,
    ROUND(atac_delta, 3)     AS atac_delta,
    ROUND(h3k27ac_delta, 3)  AS h3k27ac_delta,
    tad_id
FROM plasticity_summary
WHERE rna_delta < 0
  AND meth_delta > 0
  AND h3k27ac_delta < 0
ORDER BY rna_delta ASC
LIMIT 20;
