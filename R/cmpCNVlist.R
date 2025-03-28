#' @title Compare CNVs Between multiple samples and one control
#' @description
#'  CNV comparison between multiple samples and one control, resulting in a list of two samples compared.
#'      It can be thought of as multiple uses of the cmpCNV function.
#'
#' @param Treat A list of GRanges objects
#' @param Control A GRanges object
#' @param default_cn numeric value(The copy number of the uncovered fragment in sample B for sample A ). Default: `2`.
#'
#' @return A list of GRanges objects (each element is the result of a comparison)
#' @export
#'
#' @examples
#' gr1_file <- system.file("extdata", "test1.txt", package = "CNVProcessor")
#' if (file.exists(gr1_file)) {
#'   gr1 <- loadCNVtoGranges(gr1_file, start = "start", end = "end", cn = "cn", seqnames = "seqnames")
#' }
#' gr2_file <- system.file("extdata", "test2.txt", package = "CNVProcessor")
#' if (file.exists(gr2_file)) {
#'   gr2 <- loadCNVtoGranges(gr2_file, start = "start", end = "end", cn = "cn", seqnames = "seqnames")
#' }
#' gr3_file <- system.file("extdata", "test3.txt", package = "CNVProcessor")
#' if (file.exists(gr3_file)) {
#'   gr3 <- loadCNVtoGranges(gr3_file, start = "start", end = "end", cn = "cn", seqnames = "seqnames")
#' }
#' cmpCNVlist_result <- cmpCNVlist(Treat=list(gr1=gr1,gr2=gr2),Control=gr3)
#' head(cmpCNVlist_result)
cmpCNVlist <- function(
    Treat,
    Control,
    default_cn = 2
) {
  if (!is.list(Treat)) stop("Treat must be a list of GRanges objects")
  if (is.null(names(Treat))) stop("The Treat list must contain named elements")
  for (sample in Treat) {
    if (!is(sample, "GRanges")) stop("The element in the Treat list must be a GRanges object")
  }
  if (!is(Control, "GRanges")) stop("Control must be GRanges object")

  compare_cnvs <- function(sample1, sample2, default_cn = 2) {
    # Parameter validation
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

    GenomicRanges::mcols(combined_gr)$Treatment_Name <- sample_name
    GenomicRanges::mcols(combined_gr)$Control_Name <- control_name
    # Dynamically generate column names and assign values
    mcols(combined_gr)[[paste0(sample1_name, "_cn")]] <-
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
  # Get the control group name (using the variable name passed in)
  control_name <- deparse(substitute(Control))

  # Initialize the result list
  results <- list()

  # Iterate through each treatment group sample
  for (sample_name in names(Treat)) {
    sample_gr <- Treat[[sample_name]]

    # Call cmpCNV to compare
    comparison <- compare_cnvs(
      sample1 = sample_gr,
      sample2 = Control,
      default_cn = default_cn
    )

    # Save the results to a list
    results[[paste0(sample_name, "vs", control_name)]] <- comparison
  }

  return(results)
}
