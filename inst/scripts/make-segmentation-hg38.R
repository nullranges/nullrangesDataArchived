library(excluderanges)
suppressMessages(library(AnnotationHub))
suppressPackageStartupMessages(library(GenomicRanges))
library(nullranges)
ah <- AnnotationHub()

## Derive exclude regions from ENCODE
query_data <- query(ah, c("excluderanges", "hg38", "Exclusion regions"))
query_data
excludeGR.hg38.Kundaje.1 <- query_data[["AH95917"]]

## Derive telomere,centromere from UCSC and rCGH
query_data2 <- query(ah, c("excluderanges", "UCSC", "Homo Sapiens", "hg38"))
query_data2
telomere <- query_data2[["AH95938"]]

suppressPackageStartupMessages(library(rCGH))
# hg38 # data.frame
# Adjust chromosome names
hg38$chrom[hg38$chrom == 23] <- "X"
hg38$chrom[hg38$chrom == 24] <- "Y"
hg38$chrom <- paste0("chr", hg38$chrom)
# Make GRanges object
hg38.UCSC.centromere <- makeGRangesFromDataFrame(hg38, seqnames.field = "chrom", start.field = "centromerStart", end.field = "centromerEnd")
# Assign seqinfo data
seqlengths(hg38.UCSC.centromere) <- hg38$length
genome(hg38.UCSC.centromere)     <- "hg38"
# Resulting object
hg38.UCSC.centromere

## Combining ENCODE with centromere and telomere
excludeGR.hg38.all <- reduce(c(excludeGR.hg38.Kundaje.1, telomere, hg38.UCSC.centromere))
excludeGR.hg38.all
summary(width(excludeGR.hg38.all))
exclude <- excludeGR.hg38.all %>% plyranges::filter(width(excludeGR.hg38.all) >= 500)  ## Remove small pieces for bootstrap speed up

### Based on gene density
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
edb <- EnsDb.Hsapiens.v86
filt <- AnnotationFilterList(GeneIdFilter("ENSG", "startsWith"))
g <- genes(edb, filter = filt)

library(GenomeInfoDb)
g <- keepStandardChromosomes(g, pruning.mode = "coarse")
seqlevels(g, pruning.mode="coarse") <- setdiff(seqlevels(g), c("MT"))
# normally we would assign a new style, but for recent host issues
## seqlevelsStyle(g) <- "UCSC" 
seqlevels(g) <- paste0("chr", seqlevels(g))
genome(g) <- "hg38"

g <- sortSeqlevels(g)
g <- sort(g)

## Segmentation with tiling block length 2e6
L_s <- 2e6
seg_cbs <- segmentDensity(g, n = 3, L_s = L_s,
                      deny = exclude, type = "cbs")

seg_hmm <- segmentDensity(g, n = 3, L_s = L_s,
                          deny = exclude, type = "hmm")

save(seg_cbs, file = "data/hg38_seg2e6_cbs.rda", compress = "xz")
save(seg_hmm, file = "data/hg38_seg2e6_hmm.rda", compress = "xz")
save(excludeGR.hg38.all, file = "data/hg38_excludeall.rda", compress = "xz")

# plot <- lapply(c("ranges", "barplot", "boxplot"), function(x) plotSegment(seg_cbs,type=x,deny=exclude))
# plot[[1]]
# plot[[2]]
# plot[[3]]

