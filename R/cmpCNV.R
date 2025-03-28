#' @title Comparison of CNVs between two samples
#' @description
#' Comparison of CNVs between two samples
#'
#' @param sample1 A GRanges object
#' @param sample2 A GRanges object
#' @param default_cn numeric value(The copy number of the uncovered fragment in sample B for sample A ). Default: `2`.
#'
#' @return A GRanges object
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
#' gr_merge <- cmpCNV(sample1=gr1,sample2=gr2)
#' head(gr_merge)
cmpCNV <- function(sample1, sample2, default_cn = 2) {
  #Parameter validation
  stopifnot(
    is(sample1, "GRanges"),
    is(sample2, "GRanges"),
    "cn" %in% colnames(mcols(sample1)),
    "cn" %in% colnames(mcols(sample2))
  )

  # Get the name of the input parameter
  sample1_name <- deparse(substitute(sample1))
  sample2_name <- deparse(substitute(sample2))

  # The intervals of the two samples are combined and cut into non-overlapping intervals
  combined_gr <- GenomicRanges::disjoin(c(sample1, sample2))

  # Define function: Obtain the CN value of the target interval in the reference sample
  get_cn_values <- function(query_gr, ref_gr, default_cn) {
    hits <- GenomicRanges::findOverlaps(query_gr, ref_gr)
    cn_vals <- rep(default_cn, length(query_gr))
    cn_vals[S4Vectors::queryHits(hits)] <- ref_gr$cn[S4Vectors::subjectHits(hits)]
    cn_vals
  }

  # Dynamically generate column names and assign values
  GenomicRanges::mcols(combined_gr)[[paste0(sample1_name, "_cn")]] <-
    get_cn_values(combined_gr, sample1, default_cn)

  GenomicRanges::mcols(combined_gr)[[paste0(sample2_name, "_cn")]] <-
    get_cn_values(combined_gr, sample2, default_cn)

  # Determining the Type of Change
  combined_gr$FC <- GenomicRanges::mcols(combined_gr)[[paste0(sample2_name, "_cn")]]-GenomicRanges::mcols(combined_gr)[[paste0(sample1_name, "_cn")]]

  combined_gr$change <- ifelse(
    GenomicRanges::mcols(combined_gr)[[paste0(sample2_name, "_cn")]] >
      GenomicRanges::mcols(combined_gr)[[paste0(sample1_name, "_cn")]],
    "amplification",
    ifelse(
      GenomicRanges::mcols(combined_gr)[[paste0(sample2_name, "_cn")]] <
        GenomicRanges::mcols(combined_gr)[[paste0(sample1_name, "_cn")]],
      "deletion",
      "no_change"
    )
  )

  # Returns the sorted results
  sort(combined_gr)
}
