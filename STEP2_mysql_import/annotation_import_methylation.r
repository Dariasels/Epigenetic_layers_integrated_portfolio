#BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
write.csv(
  anno[, c("Name","chr","pos","UCSC_RefGene_Name","UCSC_RefGene_Group")],
  "GPL13534_annotation.csv",
  row.names = FALSE
)
