meta <- data.frame(
  Title = c("","",""),
  Description = c("","",""),
  BiocVersion = c("3.14","3.14","3.14"),
  Genome = c("hg38","hg19","hg19"),
  SourceType = c("","",""),
  SourceUrl = c("","",""),
  SourceVersion = c("","",""),
  Species = rep("Homo sapiens", 3),
  TaxonomyId = rep(9606, 3),
  Coordinate_1_based = TRUE,
  DataProvider = c("","",""),
  Maintainer = rep("Michael Love <michaelisaiahlove@gmail.com>", 3)
  RDataClass = c("GenomicRanges","GenomicRanges","InteractionSet")
  DispatchClass = rep("Rda", 3),
  RDataPath = c("","",""),
  Tags = c("","","")
  Notes = c("","",""))

write.csv(meta, file="metadata.csv", row.names=FALSE)
