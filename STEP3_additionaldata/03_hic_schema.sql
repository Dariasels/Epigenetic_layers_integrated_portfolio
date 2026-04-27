-- ============================================================================
-- HIC Integration Schema
-- ============================================================================
-- Adds HIC 3D chromatin contact data to the multi-omics integration.
-- Enables analysis of long-range chromatin interactions in AD vs control.
--
-- Tables created:
--   - hic_loops: Raw contact loops with metadata
--   - hic_gene_links: Maps HIC loop endpoints to genes
--   - hic_enhancer_links: Links HIC loops to H3K27ac enhancers
--   - hic_plasticity_summary: Query view of plasticity gene HIC contacts
--
-- Run after: 02_chip_tad_schema.sql
-- Usage: mysql -u daria -p brain_multiomics < 03_hic_schema.sql

DROP VIEW IF EXISTS hic_plasticity_contact_graph;
DROP TABLE IF EXISTS hic_enhancer_links;
DROP TABLE IF EXISTS hic_gene_links;
DROP TABLE IF EXISTS hic_loops;

-- ============================================================================
-- HIC_LOOPS: Raw 3D contact data
-- ============================================================================
-- Each row represents one significant contact between two genomic loci.
-- Contact strength normalized by library size / genomic distance.

CREATE TABLE IF NOT EXISTS hic_loops (
    hic_id INT PRIMARY KEY AUTO_INCREMENT,
    chrom1 VARCHAR(10) NOT NULL,
    start1 INT NOT NULL,
    end1 INT NOT NULL,
    chrom2 VARCHAR(10) NOT NULL,
    start2 INT NOT NULL,
    end2 INT NOT NULL,
    contact_strength FLOAT NOT NULL DEFAULT 1.0,  -- Normalized contact frequency
    source_dataset VARCHAR(100) NOT NULL,          -- e.g., "GSE105194", "AD_prefrontal"
    resolution_kb INT DEFAULT 25,                   -- Binning resolution (25kb typical)
    INDEX idx_chrom1 (chrom1, start1, end1),
    INDEX idx_chrom2 (chrom2, start2, end2),
    INDEX idx_dataset (source_dataset),
    INDEX idx_strength (contact_strength DESC)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

-- ============================================================================
-- HIC_GENE_LINKS: Maps HIC loop endpoints to genes
-- ============================================================================
-- For each HIC loop, identify genes overlapping or near endpoints.
-- Enables queries like: "Which genes physically interact in 3D space?"

CREATE TABLE IF NOT EXISTS hic_gene_links (
    hic_link_id INT PRIMARY KEY AUTO_INCREMENT,
    hic_id INT NOT NULL,
    gene_symbol VARCHAR(100) NOT NULL,
    interacting_region VARCHAR(50) NOT NULL,  -- 'loop_end1', 'loop_end2', 'within_loop'
    distance_to_bp INT DEFAULT 0,             -- Distance if not directly overlapping
    FOREIGN KEY (hic_id) REFERENCES hic_loops(hic_id) ON DELETE CASCADE,
    INDEX idx_gene (gene_symbol),
    INDEX idx_hic (hic_id),
    INDEX idx_region (interacting_region)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

-- ============================================================================
-- HIC_ENHANCER_LINKS: Maps HIC loops to H3K27ac enhancers
-- ============================================================================
-- Identifies 3D contacts between active enhancers (H3K27ac peaks).
-- Supports: "Which enhancers interact with each other in AD?"

CREATE TABLE IF NOT EXISTS hic_enhancer_links (
    hic_enh_id INT PRIMARY KEY AUTO_INCREMENT,
    hic_id INT NOT NULL,
    enhancer_id INT NOT NULL,
    loop_end_involved VARCHAR(10) NOT NULL,  -- 'end1', 'end2', 'both'
    FOREIGN KEY (hic_id) REFERENCES hic_loops(hic_id) ON DELETE CASCADE,
    FOREIGN KEY (enhancer_id) REFERENCES enhancers(enhancer_id),
    INDEX idx_hic (hic_id),
    INDEX idx_enhancer (enhancer_id),
    INDEX idx_loop_end (loop_end_involved)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

-- ============================================================================
-- HIC_PLASTICITY_CONTACT_GRAPH: View of plasticity gene 3D interactions
-- ============================================================================
-- Multi-layer view: HIC + genes + plasticity classification
-- Query: Which plasticity genes interact with each other or enhancers in 3D?

CREATE OR REPLACE VIEW hic_plasticity_contact_graph AS
SELECT
    hgl.hic_id,
    hgl.gene_symbol,
    hgl.interacting_region,
    hl.chrom1,
    hl.start1,
    hl.end1,
    hl.chrom2,
    hl.start2,
    hl.end2,
    hl.contact_strength,
    hl.source_dataset,
    pg.gene_symbol IS NOT NULL AS is_plasticity_gene,
    pg.category AS source_category,
    CASE
        WHEN hgl.interacting_region = 'loop_end1' THEN CONCAT(hl.chrom1, ':', hl.start1, '-', hl.end1)
        WHEN hgl.interacting_region = 'loop_end2' THEN CONCAT(hl.chrom2, ':', hl.start2, '-', hl.end2)
        ELSE 'within_loop'
    END AS contact_region
FROM hic_gene_links hgl
JOIN hic_loops hl ON hgl.hic_id = hl.hic_id
LEFT JOIN plasticity_genes pg ON hgl.gene_symbol = pg.gene_symbol
ORDER BY hl.contact_strength DESC;
