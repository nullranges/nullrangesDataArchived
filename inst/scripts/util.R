## Define helper functions for data generating scripts (from hictoolsr)

## Convert data.frames to GInteractions --------------------------------------------------
#' Convert DataFrames to GInteraction objects
#'
#' makeGInteractionsFromDataFrame takes a paired-interaction (i.e. BEDPE) formatted data-frame-like object and converts it to a GInteractions object. For convenience, \code{as_ginteractions()} can be used as an alias.
#'
#' @param df A data.table, data.frame, or DataFrame object. Assumes that the first 6   colummns are in the format chr1, start1, end1 and chr2, start2, end2, representing each pair of interactions.
#' @param keep.extra.columns TRUE or FALSE (the default). If TRUE, the columns in df that are not used to form the genomic ranges of the returned GRanges object are then returned as metadata columns on the object. Otherwise, they are ignored. If df has a width column, then it's always ignored.
#' @param starts.in.df.are.0based TRUE or FALSE (the default). If TRUE, then the start positions of the genomic ranges in df are considered to be 0-based and are converted to 1-based in the returned GRanges object. This feature is intended to make it more convenient to handle input that contains data obtained from resources using the "0-based start" convention. A notorious example of such resource is the UCSC Table Browser (http://genome.ucsc.edu/cgi-bin/hgTables).
#'
#' @return GInteraction object
#'
#' @rdname makeGInteractionsFromDataFrame
#' @aliases as_ginteractions
#'
#' @examples
#' ## data.frame
#' df <- data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
#'                  chr2 = "chr1", y1 = 30000, y2 = 40000)
#' makeGInteractionsFromDataFrame(df)
#'
#' ## data.table
#' df <- data.table(chr1 = "chr1", x1 = 10000, x2 = 20000,
#'                  chr2 = "chr1", y1 = 30000, y2 = 40000)
#' makeGInteractionsFromDataFrame(df)
#'
#' ## DataFrame
#' df <- DataFrame(chr1 = "chr1", x1 = 10000, x2 = 20000,
#'                 chr2 = "chr1", y1 = 30000, y2 = 40000)
#' makeGInteractionsFromDataFrame(df)
#'
#' ## Alias
#' df <- data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
#'                  chr2 = "chr1", y1 = 30000, y2 = 40000,
#'                  pval = 0.05, dist = 10000)
#' as_ginteractions(df)
#'
#' ## Additional metadata
#' df <- data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
#'                  chr2 = "chr1", y1 = 30000, y2 = 40000,
#'                  pval = 0.05, dist = 10000)
#' as_ginteractions(df)
#'
#' ## Remove additional metadata
#' as_ginteractions(df, keep.extra.columns = F)
#'
#' ## Add 1 to starts (for 0-based programs)
#' as_ginteractions(df, starts.in.df.are.0based = T)
#'
#' ## Bad usage
#' df <- data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
#'                  chr2 = "chr1", y1 = 30000)
#' makeGInteractionsFromDataFrame(df)
#'
makeGInteractionsFromDataFrame <- function(df,
                                           keep.extra.columns = TRUE,
                                           starts.in.df.are.0based = FALSE) {
  
  ## Convert data.table/data.frame to DataFrame
  if ("data.frame" %in% class(df)) {
    df <- DataFrame(df)
  } else if ("DFrame" %in% class(df)) {
    df <- df
  } else {
    stop("class(df) must be either 'data.frame', 'data.table', or 'DFrame'.")
  }
  
  ## Handle improper dimensions
  if(ncol(df) < 6) {
    stop("ncol(df) must be >= 6 and start with paired interactions (i.e. chr1, start1, end1 and chr2, start2, end2).")
  }
  
  ## Split into anchors
  a1 <- df[1:3] %>% `colnames<-`(c('seqnames', 'start', 'end'))
  a2 <- df[4:6] %>% `colnames<-`(c('seqnames', 'start', 'end'))
  
  ## Convert anchors to GRanges
  a1 <- makeGRangesFromDataFrame(a1, starts.in.df.are.0based = starts.in.df.are.0based)
  a2 <- makeGRangesFromDataFrame(a2, starts.in.df.are.0based = starts.in.df.are.0based)
  
  ## Create GInteractions object
  gi <- GInteractions(a1, a2)
  
  ## Add in metadata columns
  if (keep.extra.columns & ncol(df) > 6) {
    mcols(gi) <- df[7:ncol(df)]
  }
  
  ## Return
  return(gi)
  
}

#' @rdname makeGInteractionsFromDataFrame
as_ginteractions <- makeGInteractionsFromDataFrame

#' Define helper function for flexibly shifting an anchor
#'
#' @param a GRanges object
#' @param p Position within anchors to resize the bin. Can be a character or integer vector of length 1 or length(a) designating the position for each element in bedpe. Character options are "start", "end" and "center". Integers are referenced from the start position for '+' and '*' strands and from the end position for the '-' strand.
#'
#' @return GRanges object with a single position range that has been shifted appropriately.
#'
#' @rdname shiftAnchor
#'
#' @examples
#' library(GenomicRanges)
#'
#' ## Create example GRanges
#' gr1 <- GRanges(seqnames = "chr1",
#'                ranges = IRanges(start = rep(5000,3), end = rep(6000,3)),
#'                strand = c('+', '-', '*'))
#'
#' gr2 <- gr1 |> promoters(upstream = 2000, downstream = 200)
#'
#' ## Shifting anchors by keyword
#' shiftAnchor(gr1, 'start')
#' shiftAnchor(gr1, 'end')
#' shiftAnchor(gr1, 'center')
#' # shiftAnchor(gr1, 'blah') error
#'
#' ## Shifting anchors by position
#' shiftAnchor(gr1, 100)
#' shiftAnchor(gr1, c(100, 200, 300))
#' # shiftAnchor(gr1, c(100, 200, 300, 400)) error
#' # shiftAnchor(gr1, c(100, 200)) error
#'
#' ## Shifting back to TSS
#' shiftAnchor(gr2, 2000)
#'
#'
shiftAnchor <- function(a, p) {
  
  ## Shift anchors appropriately
  if (class(p) == "character") {
    
    if (p %in% c('start', 'end', 'center')) {
      
      a %<>% resize(width = 1, fix = p)
      
    } else {
      
      stop('p character must be one of "start", "end" or "center"', call. = T)
      
    }
    
  } else if (class(p) %in% c('numeric', 'integer')) {
    
    ## Convert strand Rle to vector
    sa <- as.vector(strand(a))
    
    ## Subset for strand if length(p) > 1
    if (length(p) > 1) {
      stopifnot(length(p) == length(a))
      pp <- p[which(sa %in% c('+', '*'))]
      pm <- p[which(sa == '-')] * -1
    } else {
      pp <- p
      pm <- p * -1
    }
    
    ## Shift '+' or '*' strand
    a[sa %in% c('+', '*')] %<>%
      resize(width = 1, fix = 'start') %>%
      shift(shift = pp)
    
    ## Shift '-' strand
    a[sa == '-'] %<>%
      resize(width = 1, fix = 'start') %>%
      shift(shift = pm)
    
  } else {
    
    stop('class(p) must be either character or numeric.', call. = T)
    
  }
  
  ## Return anchor
  return(a)
  
}


#' Define helper function for binning an anchor
#'
#' @param a GRanges object
#' @param p Position within anchors to resize the bin. Can be a character or integer vector of length 1 or length(a) designating the position for each element in bedpe. Character options are "start", "end" and "center". Integers are referenced from the start position for '+' and '*' strands and from the end position for the '-' strand.
#' @param res Integer - resolution in which to bin the anchor.
#'
#' @return GRanges object that has been shifted and binned into res by p.
#'
#' @rdname binAnchor
#'
#' @examples
#' library(GenomicRanges)
#'
#' ## Create example GRanges
#' gr1 <- GRanges(seqnames = "chr1",
#'                ranges = IRanges(start = rep(5000,3), end = rep(6000,3)),
#'                strand = c('+', '-', '*'))
#'
#' gr2 <- gr1 |> promoters(upstream = 2000, downstream = 200)
#'
#' ## Binning the results
#' binAnchor(gr1, 'start', 1000)
#' binAnchor(gr1, 'end', 1000)
#' binAnchor(gr1, 'center', 1000)
#'
#' ## Bin after shifting back to TSS
#' binAnchor(gr2, 2000, 1000)
#'
#'
binAnchor <- function(a, p, res) {
  
  ## Shift, bin, and trim anchors
  a %<>%
    shiftAnchor(p) %>%
    mutate(start = floor(start/res)*res,
           end = floor(start/res)*res + res) %>%
    trim() %>%
    suppressWarnings()
  
  return(a)
  
}


#' Define function to flexibly bin bedpe data by hic resolution
#'
#' @param bedpe GInteractions or data.table object with paired interactions
#' @param res Integer - resolution in which to bin bedpe anchors
#' @param a1Pos,a2Pos Position within anchors to resize the bin. Can be a character or integer vector of length 1 or length(bedpe) designating the position for each element in bedpe. Character options are "start", "end" and "center". Integers are referenced from the start position for '+' and '*' strands and from the end position for the '-' strand.
#'
#' @return GInteractions objected binned to res by a1Pos and a2Pos.
#'
#'
#' @export
#'
binBedpe <- function(bedpe, res, a1Pos, a2Pos) {
  
  if ("data.frame" %in% class(bedpe)) {
    bedpe <- try(makeGRangesFromDataFrame(bedpe))
  }
  
  ## Extract anchors
  a1 <- anchors(bedpe, type = "first")
  a2 <- anchors(bedpe, type = "second")
  
  ## Bin anchors
  a1 <- binAnchor(a = a1, p = a1Pos, res = res)
  a2 <- binAnchor(a = a2, p = a2Pos, res = res)
  
  ## Binned GInteractions object
  gi <- GInteractions(a1, a2)
  
  ## Add back metadata
  mcols(gi) <- mcols(bedpe)
  
  ## Return binned bedpe
  return(gi)
  
}


#' Define function to calculate all paired interactions within a windowSize
#'
#' @description
#' calcPairs will calculate pairs of interactions within a rolling windowSize.
#' By default, pairs of interactions are returned as StrictGInteractions where
#' anchor1 <= anchor2 for all interactions (mode = 'strict'). For more information
#' see (?InteractionSet::`GInteractions-class`). Self interactions (e.g. 1, 1 or 2, 2)
#' and order (e.g. 1, 2 is the same as 2, 1) are ignored.
#'
#' @param gr GRanges object
#' @param windowSize integer defining the window
#' @param mode If mode="strict", a StrictGInterctions object is returned with anchor
#'  indices swapped such that anchor1 <= anchor2 for all interactions. If mode="reverse",
#'  a ReverseStrictGInterctions object is returned with anchor indices swapped such that
#'  anchor1 >= anchor2.
#'
#' @return Returns a StrictGInterations object
#'
#' @examples
#' ## Load TxDb
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' ## Generate 50Kb bins across hg19
#' library(GenomicRanges)
#' bins <- tileGenome(seqinfo(txdb), tilewidth = 50e3, cut.last.tile.in.chrom = TRUE)
#'
#' ## Calculate all 50Kb bin pairs within 1Mb
#' calcPairs(gr = bins, windowSize = 1e6)
#'
calcPairs <- function(gr, windowSize, mode = 'strict') {
  
  ## Check arguments ---------------------------------------------------------------------
  
  stopifnot(isClass('GRanges', gr))
  stopifnot(isClass('integer', windowSize))
  stopifnot(length(windowSize) == 1L)
  mode <- match.arg(mode, choices = c('strict', 'reverse'))
  
  ## Begin processing --------------------------------------------------------------------
  
  ## Sort gr
  gr <- gr %>% sort()
  
  ## Define windows by windowSize
  windows <-
    gr %>%
    resize(width = windowSize, fix = 'start') %>%
    suppressWarnings() %>%
    trim()
  
  ## Group gr by windows
  ov <-
    findOverlaps(gr, windows, type = 'within') %>%
    as.data.table() %>%
    `colnames<-`(c('gr', 'windows'))
  
  ## Set up progress bar
  mx <- round(uniqueN(ov$windows)*1.01, 0)
  pb <- txtProgressBar(min = 1, max = mx, initial = NA, style = 3)
  
  ## Iterate over combinations
  ov <- ov[, {setTxtProgressBar(pb, .GRP); .(gr[1], gr[-1])}, by = windows]
  
  ## Remove NA's (last value for each chromosome)
  ov <- na.omit(ov)
  
  ## Convert coordinates to GInteractions object
  gi <- GInteractions(anchor1 = ov$V1,
                      anchor2 = ov$V2,
                      regions = gr,
                      mode = mode)
  
  ## Finish progress bar
  setTxtProgressBar(pb, value = mx)
  close(pb)
  cat('\n')
  
  return(gi)
}