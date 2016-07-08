# function to calculate binned sums from a coverage RleList and GRanges object
# adapted from http://bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesHOWTOs.pdf

binnedSum <- function(bins, numvar, mcolname)
{
  stopifnot(is(bins, "GRanges"))
  stopifnot(is(numvar, "RleList"))
  stopifnot(identical(seqlevels(bins), names(numvar)))
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  means_list <- lapply(names(numvar),
                       function(seqname) {
                         views <- Views(numvar[[seqname]],
                                        bins_per_chrom[[seqname]])
                         viewSums(views)
                       })
  new_mcol <- unsplit(means_list, as.factor(seqnames(bins)))
  mcols(bins)[[mcolname]] <- new_mcol
  bins
}

