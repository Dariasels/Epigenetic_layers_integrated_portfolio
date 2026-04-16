USE epigenetic_manifest;
show databases;
show tables;

select * from plasticity_multiomic_deltas;
select * from plasticity_genes; ##yes
select * from atac_gene_links; ##yes
select * from atac_peaks; ##yes
select * from enhancer_gene_links; ##no
select * from enhancers; ##no
select * from gene_coordinates; ##yes
select * from genes; ##no
select * from hic_loops; ##no
select * from metadata_import; ##yes
select * from methylation; ##yes
select * from methylation_gene_links; ##yes
select * from rna_expression; ##yes
select * from samples; ##yess


-- SET FOREIGN_KEY_CHECKS = 0;
-- TRUNCATE TABLE rna_expression;
-- SET FOREIGN_KEY_CHECKS = 1;

  
-- ALTER TABLE metadata_import
--   DROP COLUMN disease,
--   DROP COLUMN health_state,
--   DROP COLUMN brain_region,
--   DROP COLUMN tissue;
--   ALTER TABLE metadata_import
--   ADD COLUMN `condition` VARCHAR(20);
  
  
  SELECT 'RNA' AS layer, COUNT(DISTINCT gene_symbol) AS gene_count 
FROM rna_expression 
WHERE gene_symbol IN (SELECT gene_name FROM TargetGene)

UNION ALL

SELECT 'ATAC' AS layer, COUNT(DISTINCT gene_symbol) AS gene_count 
FROM atac_expression 
WHERE gene_symbol IN (SELECT gene_name FROM TargetGene)

UNION ALL

SELECT 'Methylation' AS layer, COUNT(DISTINCT gene_symbol) AS gene_count 
FROM methylation_expression 
WHERE gene_symbol IN (SELECT gene_name FROM TargetGene);
-- How many plasticity genes appear in each layer?
SELECT
    'RNA expression'  AS layer,
    COUNT(DISTINCT re.gene_symbol) AS genes_with_data
FROM rna_expression re
INNER JOIN plasticity_genes pg ON pg.gene_symbol = re.gene_symbol
UNION ALL
SELECT
    'ATAC-seq',
    COUNT(DISTINCT agl.gene_symbol)
FROM atac_gene_links agl
INNER JOIN plasticity_genes pg ON pg.gene_symbol = agl.gene_symbol
UNION ALL
SELECT
    'Methylation',
    COUNT(DISTINCT mgl.gene_symbol)
FROM methylation_gene_links mgl
INNER JOIN plasticity_genes pg ON pg.gene_symbol = mgl.gene_symbol;

###result: RNA expr: 365; ATAC-seq:412;Mehtylation:398!!


  -- Check specifically for BDNF -- does it exist in rna_expression?
SELECT DISTINCT gene_symbol FROM rna_expression 
WHERE gene_symbol LIKE '%BDNF%' 
   OR gene_symbol LIKE '%bdnf%'
   OR gene_symbol LIKE '%Bdnf%';
 

-- Which plasticity genes are covered by ALL three layers? 344
SELECT pg.gene_symbol, pg.category
FROM plasticity_genes pg
WHERE pg.gene_symbol IN (SELECT DISTINCT gene_symbol FROM rna_expression)
  AND pg.gene_symbol IN (SELECT DISTINCT gene_symbol FROM atac_gene_links)
  AND pg.gene_symbol IN (SELECT DISTINCT gene_symbol FROM methylation_gene_links)
ORDER BY pg.category, pg.gene_symbol;


## 356 genes overlap the three layers!
SELECT COUNT(*) AS triple_layer_count
FROM plasticity_genes pg
WHERE pg.gene_symbol IN (SELECT DISTINCT gene_symbol FROM rna_expression)
  AND pg.gene_symbol IN (SELECT DISTINCT gene_symbol FROM atac_gene_links)
  AND pg.gene_symbol IN (SELECT DISTINCT gene_symbol FROM methylation_gene_links);


## select plasticity genes for remaining layers
-- RNA: how many rows are plasticity genes? '556818'

SELECT COUNT(*) AS overlapping_rows
FROM rna_expression re
INNER JOIN plasticity_genes pg ON re.gene_symbol = pg.gene_symbol;

-- ATAC: how many peaks link to plasticity genes? '4603', '412'

SELECT COUNT(DISTINCT agl.peak_id) AS plasticity_peaks,
       COUNT(DISTINCT agl.gene_symbol) AS plasticity_genes
FROM atac_gene_links agl
INNER JOIN plasticity_genes pg ON pg.gene_symbol = agl.gene_symbol;


###QUESTION 1: which plastciity genes show DE In AD ?  
SELECT re.gene_symbol, pg.category,
    AVG(CASE WHEN UPPER(s.condition) = 'Alzheimer'      THEN re.expression_value END) AS mean_AD,
    AVG(CASE WHEN UPPER(s.condition) = 'CONTROL' THEN re.expression_value END) AS mean_ctrl
FROM rna_expression re
JOIN plasticity_genes pg ON re.gene_symbol = pg.gene_symbol
JOIN samples s           ON re.sample_id   = s.sample_id
GROUP BY re.gene_symbol, pg.category
HAVING mean_AD IS NOT NULL AND mean_ctrl IS NOT NULL
ORDER BY ABS(mean_AD - mean_ctrl) DESC
LIMIT 100;

SELECT COUNT(*) AS total_plasticity_genes_with_data
FROM (
    SELECT re.gene_symbol
    FROM rna_expression re
    JOIN plasticity_genes pg ON re.gene_symbol = pg.gene_symbol
    JOIN samples s           ON re.sample_id   = s.sample_id
    GROUP BY re.gene_symbol, pg.category
    HAVING AVG(CASE WHEN UPPER(s.condition) = 'Alzheimer' THEN re.expression_value END) IS NOT NULL
       AND AVG(CASE WHEN UPPER(s.condition) = 'CONTROL'   THEN re.expression_value END) IS NOT NULL
) AS validated_genes;
## 365

-- This should return a count > 0
-- SELECT COUNT(*) 
-- FROM rna_expression re 
-- JOIN samples s ON re.sample_id = s.sample_id;

###Q2 — Which plasticity gene promoters are hypermethylated in AD? 141
SELECT mgl.gene_symbol, pg.category,
    AVG(CASE WHEN s.condition='Alzheimer'      THEN m.beta_value END) AS beta_AD,
    AVG(CASE WHEN s.condition='control' THEN m.beta_value END) AS beta_ctrl
FROM methylation m
JOIN methylation_gene_links mgl ON mgl.cpg_id = m.cpg_id
JOIN plasticity_genes pg        ON pg.gene_symbol = mgl.gene_symbol
JOIN samples s                  ON s.sample_id = m.sample_id
WHERE mgl.relation IN ('TSS200','TSS1500')
GROUP BY mgl.gene_symbol, pg.category
HAVING beta_AD > beta_ctrl
ORDER BY (beta_AD - beta_ctrl) DESC;

#Q3 — Which plasticity genes have ALL THREE layers pointing the same direction in AD?
#(hypermethylated promoter + closed chromatin + lower expression = epigenetically silenced)
#vey long - never finished

SELECT pg.gene_symbol, pg.category,
    AVG(CASE WHEN s_r.condition='AD' THEN re.expression_value END)      AS expr_AD,
    AVG(CASE WHEN s_r.condition='control' THEN re.expression_value END) AS expr_ctrl,
    AVG(CASE WHEN s_m.condition='AD' THEN m.beta_value END)             AS meth_AD,
    AVG(CASE WHEN s_a.condition='AD' THEN ap.signal_value END)          AS atac_AD
FROM plasticity_genes pg
JOIN rna_expression re          ON re.gene_symbol  = pg.gene_symbol
JOIN samples s_r                ON s_r.sample_id   = re.sample_id
JOIN methylation_gene_links mgl ON mgl.gene_symbol = pg.gene_symbol
JOIN methylation m              ON m.cpg_id        = mgl.cpg_id
JOIN samples s_m                ON s_m.sample_id   = m.sample_id
JOIN atac_gene_links agl        ON agl.gene_symbol = pg.gene_symbol
JOIN atac_peaks ap              ON ap.peak_id      = agl.peak_id
JOIN samples s_a                ON s_a.sample_id   = ap.sample_id
WHERE mgl.relation IN ('TSS200','TSS1500')
  AND agl.region_type = 'promoter'
GROUP BY pg.gene_symbol, pg.category;

#werkt ni
WITH TargetGene AS (
    SELECT 'BDNF' as gene_name -- <--- Change your gene here
),
RNA_Avg AS (
    SELECT 
        AVG(CASE WHEN s.condition = 'AD' THEN re.expression_value END) as rna_ad,
        AVG(CASE WHEN s.condition = 'control' THEN re.expression_value END) as rna_ctrl
    FROM rna_expression re
    JOIN samples s ON re.sample_id = s.sample_id
    WHERE re.gene_symbol = (SELECT gene_name FROM TargetGene)
),
Meth_Avg AS (
    SELECT 
        AVG(CASE WHEN s.condition = 'AD' THEN m.beta_value END) as meth_ad,
        AVG(CASE WHEN s.condition = 'control' THEN m.beta_value END) as meth_ctrl
    FROM methylation m
    JOIN samples s ON m.sample_id = s.sample_id
    JOIN methylation_gene_links mgl ON m.cpg_id = mgl.cpg_id
    WHERE mgl.gene_symbol = (SELECT gene_name FROM TargetGene)
    AND mgl.relation IN ('TSS200', 'TSS1500')
),
ATAC_Avg AS (
    SELECT 
        AVG(CASE WHEN s.condition = 'AD' THEN ap.signal_value END) as atac_ad,
        AVG(CASE WHEN s.condition = 'control' THEN ap.signal_value END) as atac_ctrl
    FROM atac_peaks ap
    JOIN samples s ON ap.sample_id = s.sample_id
    JOIN atac_gene_links agl ON ap.peak_id = agl.peak_id
    WHERE agl.gene_symbol = (SELECT gene_name FROM TargetGene)
    AND agl.region_type = 'promoter'
)
SELECT 
    (SELECT gene_name FROM TargetGene) AS gene,
    rna_ad, rna_ctrl, 
    meth_ad, meth_ctrl, 
    atac_ad, atac_ctrl
FROM RNA_Avg, Meth_Avg, ATAC_Avg;

-- Are any plasticity genes linked to ATAC peaks?
SELECT pg.gene_symbol, pg.category,
       COUNT(agl.peak_id)      AS n_peaks,
       MIN(agl.distance_bp)    AS closest_peak_bp,
       GROUP_CONCAT(DISTINCT agl.region_type) AS region_types
FROM plasticity_genes pg
JOIN atac_gene_links agl ON agl.gene_symbol = pg.gene_symbol
GROUP BY pg.gene_symbol, pg.category
ORDER BY n_peaks DESC
LIMIT 20;

-- What region types do you have overall?
SELECT region_type, COUNT(*) AS n_links
FROM atac_gene_links
GROUP BY region_type;

-- How many total peaks, samples, and gene links do you have?
SELECT
  (SELECT COUNT(*)        FROM atac_peaks)     AS total_peaks,
  (SELECT COUNT(DISTINCT sample_id) FROM atac_peaks) AS n_samples,
  (SELECT COUNT(*)        FROM atac_gene_links) AS total_gene_links,
  (SELECT COUNT(DISTINCT gene_symbol) FROM atac_gene_links) AS unique_genes;
  
  
  
  ##everything!
  #1st we need to idex!
--   SHOW INDEX FROM rna_expression;
  -- This tells MySQL to build the index without locking the whole table
##now index the other tables: 
-- ALTER TABLE rna_expression 
-- ADD INDEX idx_rna_gene (gene_symbol),
-- ADD INDEX idx_rna_sample (sample_id),
-- ALGORITHM=INPLACE, LOCK=NONE;
-- ##this worked
-- -- Index the Link table (Gene Symbol to Peak ID)
-- ALTER TABLE atac_gene_links ADD INDEX idx_agl_gene (gene_symbol);
-- ALTER TABLE atac_gene_links ADD INDEX idx_agl_peak (peak_id);

-- -- Index the Data table (Sample ID and Peak ID)
-- ALTER TABLE atac_peaks ADD INDEX idx_ap_sample (sample_id);
-- ALTER TABLE atac_peaks ADD INDEX idx_ap_peak (peak_id);

## make metadata table ready to work
-- Ensure the sample_id and the filtering columns are fast
ALTER TABLE samples ADD INDEX idx_samples_id (sample_id);
ALTER TABLE samples ADD INDEX idx_samples_tissue (tissue);

SHOW INDEX FROM methylation;
SHOW INDEX FROM methylation_gene_links;
show index from atac_gene_links;

-- Methylation Bridge
ALTER TABLE methylation_gene_links ADD INDEX idx_mgl_gene (gene_symbol);
ALTER TABLE methylation_gene_links ADD INDEX idx_mgl_cpg (cpg_id);

-- ATAC Bridge
ALTER TABLE atac_gene_links ADD INDEX idx_agl_gene (gene_symbol);
ALTER TABLE atac_gene_links ADD INDEX idx_agl_peak (peak_id);
#done
  ALTER TABLE atac_peaks ADD INDEX idx_ap_sample (sample_id);
ALTER TABLE atac_peaks ADD INDEX idx_ap_peak (peak_id);

show index from atac_peaks;

  -- -- Delete the old ones (don't worry if they don't exist, it will just give a warning)
-- DROP INDEX idx_rna_gene ON rna_expression;
-- DROP INDEX idx_rna_sample ON rna_expression;



#####This is a Meta-Analytic Multi-Omic Profile. yeas worksss!!
SELECT 
    pg.gene_symbol, 
    pg.category,
    -- RNA Stats (Average of all AD samples vs all Control samples)
    (SELECT AVG(re.expression_value) 
     FROM rna_expression re 
     JOIN samples s ON re.sample_id = s.sample_id 
     WHERE re.gene_symbol = pg.gene_symbol 
       AND s.condition = 'Alzheimer' 
       AND s.tissue LIKE '%prefrontal cortex%') AS rna_AD,

    (SELECT AVG(re.expression_value) 
     FROM rna_expression re 
     JOIN samples s ON re.sample_id = s.sample_id 
     WHERE re.gene_symbol = pg.gene_symbol 
       AND s.condition = 'Control' 
       AND s.tissue LIKE '%prefrontal cortex%') AS rna_Ctrl,

    -- Methylation Stats (Average of all AD samples vs all Control samples)
    (SELECT AVG(m.beta_value) 
     FROM methylation m 
     JOIN methylation_gene_links mgl ON m.cpg_id = mgl.cpg_id
     JOIN samples s ON m.sample_id = s.sample_id 
     WHERE mgl.gene_symbol = pg.gene_symbol 
       AND s.condition = 'Alzheimer' 
       AND s.tissue LIKE '%prefrontal cortex%') AS meth_AD,

    (SELECT AVG(m.beta_value) 
     FROM methylation m 
     JOIN methylation_gene_links mgl ON m.cpg_id = mgl.cpg_id
     JOIN samples s ON m.sample_id = s.sample_id 
     WHERE mgl.gene_symbol = pg.gene_symbol 
       AND s.condition = 'Control' 
       AND s.tissue LIKE '%prefrontal cortex%') AS meth_Ctrl

FROM plasticity_genes pg
-- Let's look at the first 20 to see the trends
LIMIT 20;

##include atac
SELECT 
    pg.gene_symbol, 
    pg.category,

    -- RNA
    (SELECT AVG(re.expression_value) 
     FROM rna_expression re 
     JOIN samples s ON re.sample_id = s.sample_id 
     WHERE re.gene_symbol = pg.gene_symbol 
       AND s.condition = 'Alzheimer' 
       AND s.tissue LIKE '%prefrontal cortex%') AS rna_AD,

    (SELECT AVG(re.expression_value) 
     FROM rna_expression re 
     JOIN samples s ON re.sample_id = s.sample_id 
     WHERE re.gene_symbol = pg.gene_symbol 
       AND s.condition = 'Control' 
       AND s.tissue LIKE '%prefrontal cortex%') AS rna_Ctrl,

    -- Methylation (all CpGs mapped to this gene, all regions)
    (SELECT AVG(m.beta_value) 
     FROM methylation m 
     JOIN methylation_gene_links mgl ON m.cpg_id = mgl.cpg_id
     JOIN samples s ON m.sample_id = s.sample_id 
     WHERE mgl.gene_symbol = pg.gene_symbol 
	   AND mgl.relation IN ('TSS200', 'TSS1500')
       AND s.condition = 'Alzheimer' 
       AND s.tissue LIKE '%prefrontal cortex%') AS meth_AD,

    (SELECT AVG(m.beta_value) 
     FROM methylation m 
     JOIN methylation_gene_links mgl ON m.cpg_id = mgl.cpg_id
     JOIN samples s ON m.sample_id = s.sample_id 
     WHERE mgl.gene_symbol = pg.gene_symbol 
       AND s.condition = 'Control' 
       AND s.tissue LIKE '%prefrontal cortex%') AS meth_Ctrl,

    -- ATAC: average signal of peaks near this gene's promoter, AD samples
    (SELECT AVG(ap.signal_value)
     FROM atac_peaks ap
     JOIN atac_gene_links agl ON ap.peak_id = agl.peak_id
     JOIN samples s ON ap.sample_id = s.sample_id
     WHERE agl.gene_symbol = pg.gene_symbol
       AND agl.region_type = 'promoter'
       AND s.condition = 'Alzheimer') AS atac_AD,

    -- ATAC: average signal of peaks near this gene's promoter, Control samples
    (SELECT AVG(ap.signal_value)
     FROM atac_peaks ap
     JOIN atac_gene_links agl ON ap.peak_id = agl.peak_id
     JOIN samples s ON ap.sample_id = s.sample_id
     WHERE agl.gene_symbol = pg.gene_symbol
       AND agl.region_type = 'promoter'
       AND s.condition = 'Control') AS atac_Ctrl

FROM plasticity_genes pg
LIMIT 50; ## this works perfectly!


--  #run for all - save to table
--  CREATE TABLE plasticity_multiomic_deltas AS
-- SELECT 
--     pg.gene_symbol, 
--     pg.category,
--     -- RNA Stats & Delta
--     ROUND(rna_results.ad_avg, 4) AS rna_AD,
--     ROUND(rna_results.ctrl_avg, 4) AS rna_Ctrl,
--     ROUND(rna_results.ad_avg - rna_results.ctrl_avg, 4) AS rna_delta,
--     
--     -- Methylation Stats & Delta
--     ROUND(meth_results.ad_avg, 4) AS meth_AD,
--     ROUND(meth_results.ctrl_avg, 4) AS meth_Ctrl,
--     ROUND(meth_results.ad_avg - meth_results.ctrl_avg, 4) AS meth_delta

-- FROM plasticity_genes pg
-- -- Subquery for RNA
-- LEFT JOIN (
--     SELECT re.gene_symbol,
--            AVG(CASE WHEN s.condition = 'Alzheimer' THEN re.expression_value END) AS ad_avg,
--            AVG(CASE WHEN s.condition = 'Control' THEN re.expression_value END) AS ctrl_avg
--     FROM rna_expression re
--     JOIN samples s ON re.sample_id = s.sample_id
--     WHERE s.tissue LIKE '%prefrontal cortex%'
--     GROUP BY re.gene_symbol
-- ) AS rna_results ON pg.gene_symbol = rna_results.gene_symbol

-- -- Subquery for Methylation
-- LEFT JOIN (
--     SELECT mgl.gene_symbol,
--            AVG(CASE WHEN s.condition = 'Alzheimer' THEN m.beta_value END) AS ad_avg,
--            AVG(CASE WHEN s.condition = 'Control' THEN m.beta_value END) AS ctrl_avg
--     FROM methylation m
--     JOIN methylation_gene_links mgl ON m.cpg_id = mgl.cpg_id
--     JOIN samples s ON m.sample_id = s.sample_id
--     WHERE s.tissue LIKE '%prefrontal cortex%'
--     GROUP BY mgl.gene_symbol
-- ) AS meth_results ON pg.gene_symbol = meth_results.gene_symbol

-- WHERE rna_results.ad_avg IS NOT NULL 
--   AND meth_results.ad_avg IS NOT NULL
-- ORDER BY rna_delta ASC;
--  
-- drop table plasticity_summary;
-- drop table plasticity_multiomic_deltas;
#make materialized summary table
-- Step 1: create the summary table
CREATE TABLE plasticity_summary (
    gene_symbol   VARCHAR(50) PRIMARY KEY,
    category      VARCHAR(100),

    rna_AD        FLOAT,
    rna_Ctrl      FLOAT,
    rna_delta     FLOAT,          -- AD minus Control (positive = higher in AD)

    meth_AD       FLOAT,
    meth_Ctrl     FLOAT,
    meth_delta    FLOAT,          -- positive = more methylated in AD

    atac_AD       FLOAT,
    atac_Ctrl     FLOAT,
    atac_delta    FLOAT,          -- positive = more open chromatin in AD

    -- convenience flags for quick filtering
    -- is_silenced   TINYINT(1),     -- 1 if lower RNA + higher meth + lower ATAC in AD
--     is_activated  TINYINT(1),     -- 1 if higher RNA + lower meth + higher ATAC in AD
    layers_with_data TINYINT(1),  -- how many of the 3 layers have non-NULL values

    last_updated  TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Step 2: populate it (this is the slow query, run once)
INSERT INTO plasticity_summary
    (gene_symbol, category,
     rna_AD, rna_Ctrl, rna_delta,
     meth_AD, meth_Ctrl, meth_delta,
     atac_AD, atac_Ctrl, atac_delta,
     layers_with_data)

SELECT
    pg.gene_symbol,
    pg.category,

    -- RNA
    (SELECT AVG(re.expression_value)
     FROM rna_expression re
     JOIN samples s ON re.sample_id = s.sample_id
     WHERE re.gene_symbol = pg.gene_symbol
       AND s.condition = 'Alzheimer'
       AND s.tissue LIKE '%prefrontal cortex%') AS rna_AD,

    (SELECT AVG(re.expression_value)
     FROM rna_expression re
     JOIN samples s ON re.sample_id = s.sample_id
     WHERE re.gene_symbol = pg.gene_symbol
       AND s.condition = 'Control'
       AND s.tissue LIKE '%prefrontal cortex%') AS rna_Ctrl,

    -- rna_delta computed inline
    (SELECT AVG(re.expression_value)
     FROM rna_expression re JOIN samples s ON re.sample_id = s.sample_id
     WHERE re.gene_symbol = pg.gene_symbol AND s.condition = 'Alzheimer'
       AND s.tissue LIKE '%prefrontal cortex%')
    -
    (SELECT AVG(re.expression_value)
     FROM rna_expression re JOIN samples s ON re.sample_id = s.sample_id
     WHERE re.gene_symbol = pg.gene_symbol AND s.condition = 'Control'
       AND s.tissue LIKE '%prefrontal cortex%') AS rna_delta,

    -- Methylation
    (SELECT AVG(m.beta_value)
     FROM methylation m
     JOIN methylation_gene_links mgl ON m.cpg_id = mgl.cpg_id
     JOIN samples s ON m.sample_id = s.sample_id
     WHERE mgl.gene_symbol = pg.gene_symbol
       AND mgl.relation IN ('TSS200', 'TSS1500')
       AND s.condition = 'Alzheimer'
       AND s.tissue LIKE '%prefrontal cortex%') AS meth_AD,

    (SELECT AVG(m.beta_value)
     FROM methylation m
     JOIN methylation_gene_links mgl ON m.cpg_id = mgl.cpg_id
     JOIN samples s ON m.sample_id = s.sample_id
     WHERE mgl.gene_symbol = pg.gene_symbol
       AND mgl.relation IN ('TSS200', 'TSS1500')
       AND s.condition = 'Control'
       AND s.tissue LIKE '%prefrontal cortex%') AS meth_Ctrl,

    (SELECT AVG(m.beta_value)
     FROM methylation m JOIN methylation_gene_links mgl ON m.cpg_id = mgl.cpg_id
     JOIN samples s ON m.sample_id = s.sample_id
     WHERE mgl.gene_symbol = pg.gene_symbol AND mgl.relation IN ('TSS200','TSS1500') AND s.condition = 'Alzheimer'
       AND s.tissue LIKE '%prefrontal cortex%')
    -
    (SELECT AVG(m.beta_value)
     FROM methylation m JOIN methylation_gene_links mgl ON m.cpg_id = mgl.cpg_id
     JOIN samples s ON m.sample_id = s.sample_id
     WHERE mgl.gene_symbol = pg.gene_symbol AND mgl.relation IN ('TSS200','TSS1500') AND s.condition = 'Control'
       AND s.tissue LIKE '%prefrontal cortex%') AS meth_delta,

    -- ATAC
    (SELECT AVG(ap.signal_value)
     FROM atac_peaks ap
     JOIN atac_gene_links agl ON ap.peak_id = agl.peak_id
     JOIN samples s ON ap.sample_id = s.sample_id
     WHERE agl.gene_symbol = pg.gene_symbol
       AND agl.region_type = 'promoter'
       AND s.condition = 'Alzheimer') AS atac_AD,

    (SELECT AVG(ap.signal_value)
     FROM atac_peaks ap
     JOIN atac_gene_links agl ON ap.peak_id = agl.peak_id
     JOIN samples s ON ap.sample_id = s.sample_id
     WHERE agl.gene_symbol = pg.gene_symbol
       AND agl.region_type = 'promoter'
       AND s.condition = 'Control') AS atac_Ctrl,

    (SELECT AVG(ap.signal_value)
     FROM atac_peaks ap JOIN atac_gene_links agl ON ap.peak_id = agl.peak_id
     JOIN samples s ON ap.sample_id = s.sample_id
     WHERE agl.gene_symbol = pg.gene_symbol AND agl.region_type = 'promoter'
       AND s.condition = 'Alzheimer')
    -
    (SELECT AVG(ap.signal_value)
     FROM atac_peaks ap JOIN atac_gene_links agl ON ap.peak_id = agl.peak_id
     JOIN samples s ON ap.sample_id = s.sample_id
     WHERE agl.gene_symbol = pg.gene_symbol AND agl.region_type = 'promoter'
       AND s.condition = 'Control') AS atac_delta,

    -- is_silenced: RNA down + methylation up + ATAC down in AD
    -- uses the delta columns computed above (recomputed inline here)
     -- placeholder, we update these below
    NULL

FROM plasticity_genes pg;

-- Step 3: fill in the convenience flags
UPDATE plasticity_summary
SET
    -- is_silenced = CASE
--         WHEN rna_delta  < 0          -- lower expression in AD
--          AND meth_delta > 0.05       -- meaningfully more methylated in AD
--          AND (atac_delta < 0         -- less open chromatin in AD
--               OR atac_AD IS NULL)    -- or no peak at all in AD
--         THEN 1 ELSE 0
--     END,

--     is_activated = CASE
--         WHEN rna_delta  > 0          -- higher expression in AD
--          AND meth_delta < -0.05      -- meaningfully less methylated in AD
--          AND atac_delta > 0          -- more open chromatin in AD
--         THEN 1 ELSE 0
--     END,

    layers_with_data = (
        CASE WHEN rna_AD  IS NOT NULL THEN 1 ELSE 0 END +
        CASE WHEN meth_AD IS NOT NULL THEN 1 ELSE 0 END +
        CASE WHEN atac_AD IS NOT NULL THEN 1 ELSE 0 END
    );
    
    
    SELECT gene_symbol, 
       ROUND(rna_delta, 4)  AS rna_delta,
       ROUND(meth_delta, 4) AS meth_delta,
       atac_delta
FROM plasticity_summary
WHERE rna_delta IS NOT NULL
ORDER BY rna_delta ASC
LIMIT 20;

-- TRUNCATE TABLE plasticity_summary;

    -- Look up one gene
SELECT * FROM plasticity_summary WHERE gene_symbol = 'BDNF';
SELECT * FROM plasticity_summary;

-- ALTER TABLE plasticity_summary
--   DROP COLUMN is_silenced,
--   DROP COLUMN is_activated;

-- Genes with lower expression AND higher methylation in AD (the meaningful signal)
SELECT gene_symbol, category,
       ROUND(rna_delta, 4)  AS rna_delta,
       ROUND(meth_delta, 4) AS meth_delta,
       ROUND(atac_delta, 4) AS atac_delta
FROM plasticity_summary
WHERE rna_delta  < 0       -- lower expression in AD
  AND meth_delta > 0       -- higher methylation in AD
--   AND layers_with_data >= 
ORDER BY rna_delta ASC
LIMIT 30;


-- Genes with higher expression AND lower methylation in AD
SELECT gene_symbol, category,
       ROUND(rna_delta, 4)  AS rna_delta,
       ROUND(meth_delta, 4) AS meth_delta,
       ROUND(atac_delta, 4) AS atac_delta
FROM plasticity_summary
WHERE rna_delta  > 0
  AND meth_delta < -0.01 #changed from zero
--   AND layers_with_data >= 2
ORDER BY rna_delta DESC
LIMIT 30;

#add foreign keys - claude

-- Check 1: rna_expression sample_ids not in samples
SELECT DISTINCT re.sample_id 
FROM rna_expression re
LEFT JOIN samples s ON re.sample_id = s.sample_id
WHERE s.sample_id IS NULL
LIMIT 10;  #empty good

-- Check 2: methylation sample_ids not in samples
SELECT DISTINCT m.sample_id
FROM methylation m
LEFT JOIN samples s ON s.sample_id = m.sample_id
WHERE s.sample_id IS NULL
LIMIT 10; #not empty :/ but that because we didnt have to filter! we only imported plasticity genes and they are for sure Prefrotnal cortex and 
-- AD or control => insert them

-- Fix Check 2: insert missing methylation sample_ids into samples - ok
INSERT IGNORE INTO samples (sample_id, `condition`, tissue, dataset)
SELECT DISTINCT m.sample_id, 'unknown', 'prefrontal cortex', 'GSE59685'
FROM methylation m
LEFT JOIN samples s ON s.sample_id = m.sample_id
WHERE s.sample_id IS NULL;
UPDATE samples s
JOIN metadata_import m ON m.sample_id = s.sample_id
SET s.`condition` = CASE
    WHEN LOWER(m.condition) LIKE '%alzheimer%' THEN 'Alzheimer'
    WHEN LOWER(m.condition) LIKE '%control%'
      OR LOWER(m.condition) LIKE '%normal%' THEN 'Control'
    ELSE m.condition
END
WHERE s.dataset = 'GSE59685'
  AND s.`condition` = 'unknown';

-- Check 3: atac_peaks sample_ids not in samples
SELECT DISTINCT ap.sample_id
FROM atac_peaks ap
LEFT JOIN samples s ON s.sample_id = ap.sample_id
WHERE s.sample_id IS NULL
LIMIT 10;  #empty

-- Check 4: atac_gene_links peak_ids not in atac_peaks
SELECT DISTINCT agl.peak_id
FROM atac_gene_links agl
LEFT JOIN atac_peaks ap ON ap.peak_id = agl.peak_id
WHERE ap.peak_id IS NULL
LIMIT 10; #empty

-- Check 5: methylation_gene_links cpg_ids not in methylation
SELECT DISTINCT mgl.cpg_id
FROM methylation_gene_links mgl
LEFT JOIN methylation m ON m.cpg_id = mgl.cpg_id
WHERE m.cpg_id IS NULL
LIMIT 10; #empty


-- niet gelukt
ALTER TABLE atac_gene_links
  ADD CONSTRAINT fk_agl_gene
  FOREIGN KEY (gene_symbol) REFERENCES plasticity_genes(gene_symbol)
  ON DELETE CASCADE ON UPDATE CASCADE;

-- ── methylation_gene_links → methylation ──────────────────────────────
ALTER TABLE methylation_gene_links
  ADD CONSTRAINT fk_mgl_cpg
  FOREIGN KEY (cpg_id) REFERENCES methylation(cpg_id)
  ON DELETE CASCADE ON UPDATE CASCADE;

-- wel gelukt: 
-- ── rna_expression → samples ──────────────────────────────────────────
ALTER TABLE rna_expression
  ADD CONSTRAINT fk_rna_sample
  FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
  ON DELETE CASCADE ON UPDATE CASCADE;

-- ── methylation → samples ─────────────────────────────────────────────
ALTER TABLE methylation
  ADD CONSTRAINT fk_meth_sample
  FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
  ON DELETE CASCADE ON UPDATE CASCADE;

-- ── atac_peaks → samples ──────────────────────────────────────────────
ALTER TABLE atac_peaks
  ADD CONSTRAINT fk_atac_sample
  FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
  ON DELETE CASCADE ON UPDATE CASCADE;

-- ── atac_gene_links → atac_peaks ──────────────────────────────────────
ALTER TABLE atac_gene_links
  ADD CONSTRAINT fk_agl_peak
  FOREIGN KEY (peak_id) REFERENCES atac_peaks(peak_id)
  ON DELETE CASCADE ON UPDATE CASCADE;
-- ── methylation_gene_links → plasticity_genes ─────────────────────────
ALTER TABLE methylation_gene_links
  ADD CONSTRAINT fk_mgl_gene
  FOREIGN KEY (gene_symbol) REFERENCES plasticity_genes(gene_symbol)
  ON DELETE CASCADE ON UPDATE CASCADE;

-- ── plasticity_summary → plasticity_genes ─────────────────────────────
ALTER TABLE plasticity_summary
  ADD CONSTRAINT fk_summary_gene
  FOREIGN KEY (gene_symbol) REFERENCES plasticity_genes(gene_symbol)
  ON DELETE CASCADE ON UPDATE CASCADE;


## add foreign keys - chat
ALTER TABLE rna_expression
ADD CONSTRAINT fk_rna_sample
FOREIGN KEY (sample_id)
REFERENCES samples(sample_id)
ON DELETE CASCADE
ON UPDATE CASCADE;

ALTER TABLE methylation
ADD CONSTRAINT fk_meth_sample
FOREIGN KEY (sample_id)
REFERENCES samples(sample_id)
ON DELETE CASCADE
ON UPDATE CASCADE;

ALTER TABLE atac_peaks
ADD CONSTRAINT fk_atac_sample
FOREIGN KEY (sample_id)
REFERENCES samples(sample_id)
ON DELETE CASCADE
ON UPDATE CASCADE;


-- Top epigenetically silenced plasticity genes in AD--
