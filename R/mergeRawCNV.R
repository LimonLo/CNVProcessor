#' @title Merge CNV region of samples
#' @description
#' Segment the CNV regions of multiple samples according to their loci, but do not populate the uncovered genomic fragments
#'
#' @param sample_list A list of GRanges objects
#' @param default_cn numeric value(The copy number of the uncovered fragment in sample B for sample A ). Default: `2`.
#'
#' @return A GRanges object containing CNV information for multiple samples in the same region
#' @export
#'
#' @examples
#' gr1_file <- system.file("extdata", "test1.txt", package = "CNVProcessor")
#' if (file.exists(gr1_file)) {
#'   gr1 <- loadCNVtoGranges(gr1_file, start = "start", end = "end", cn = "cn", seqnames = "seqnames")
#' }
#' head(gr1)
#' gr2_file <- system.file("extdata", "test2.txt", package = "CNVProcessor")
#' if (file.exists(gr2_file)) {
#'   gr2 <- loadCNVtoGranges(gr2_file, start = "start", end = "end", cn = "cn", seqnames = "seqnames")
#' }
#' head(gr2)
#' gr_merge <- mergeRawCNV(sample_list=list(gr1=gr1,gr2=gr2))
#' head(gr_merge)
mergeRawCNV <- function(sample_list, default_cn = 2) {
  # Parameter validation
  stopifnot(
    all(sapply(sample_list, function(x) is(x, "GRanges"))),
    all(sapply(sample_list, function(x) "cn" %in% colnames(mcols(x))))
  )
  # Merge GRanges objects one by one
  combined_gr <- sample_list[[1]]
  for (i in seq_along(sample_list)[-1]) {
    combined_gr <- c(combined_gr, sample_list[[i]])
  }
  combined_gr <- GenomicRanges::disjoin(combined_gr)

  # Define a function: Populate the CN value for each sample
  fill_cn_values <- function(query_gr, ref_gr, default_cn) {
    hits <- GenomicRanges::findOverlaps(query_gr, ref_gr)
    cn_vals <- rep(default_cn, length(query_gr))
    cn_vals[S4Vectors::queryHits(hits)] <- ref_gr$cn[S4Vectors::subjectHits(hits)]
    cn_vals
  }

  # Fill in CN values for each sample
  mcols_data <- base::lapply(sample_list, function(sample) {
    fill_cn_values(combined_gr, sample, default_cn)
  })
  GenomicRanges::mcols(combined_gr) <- S4Vectors::DataFrame(setNames(mcols_data, names(sample_list)))

  # Returns the sorted results
  sort(combined_gr)
}
