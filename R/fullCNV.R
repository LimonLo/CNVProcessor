#' @title Fill Missing Chromosomal Segments
#'
#' @description
#' Fill Missing Chromosomal Segments with Default Copy Number (CN=2).
#'     This function takes a GRanges object of copy number variations (CNVs) and fills the genomic regions
#'     not covered by the input CNV segments with a default copy number (CN=2). It ensures chromosomal
#'     length consistency with a specified genome version (hg19/hg38) or user-provided chromosome lengths.
#'
#' @param gr A `GRanges` object containing CNV segments. Metadata should include a "cn" column
#'               representing copy number values.
#' @param genome_version Genome version to use for chromosome lengths. Options are `"hg19"` or `"hg38"`.
#'               Ignored if `seqlengths` is provided.
#' @param seqlengths Optional named numeric vector of chromosome lengths (e.g., `c(chr1=249250621, ...)`).
#'               If `NULL`, chromosome lengths are automatically loaded from the specified `genome_version`.
#'
#' @return A GRanges object with a continuous chromosomal segment
#' @export
#'
#' @examples
#'  # Load sample data
#' cnv_file <- system.file("extdata", "test1.txt", package = "CNVProcessor")
#' if (file.exists(cnv_file)) {
#'   gr <- loadCNVtoGranges(cnv_file, start = "start", end = "end", cn = "cn", seqnames = "seqnames")
#' }
#' head(gr)
#' # Example 1: Using genome version hg19
#' gr_full <- fullCNV(gr=gr,genome_version = "hg19")
#' head(gr_full)
#' # Example 2: Providing custom chromosome lengths
#' custom_seqlengths <- c(chr1=248956422, chr2=242193529)
#' gr_full <- fullCNV(gr=gr,seqlengths = custom_seqlengths)
fullCNV <- function(gr, genome_version = c("hg19", "hg38"),seqlengths = NULL) {

  if (!is.null(seqlengths)) {
    seqlengths_hg <- seqlengths
  } else {
    genome_version <- match.arg(genome_version)
    stopifnot(is(gr, "GRanges"))

    pkg_name <- switch(
      genome_version,
      "hg19" = "BSgenome.Hsapiens.UCSC.hg19",
      "hg38" = "BSgenome.Hsapiens.UCSC.hg38"
    )

    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      stop(
        "Requires install  ", pkg_name)
    }

    # Get genomic objects
    bs_genome <- base::getExportedValue(pkg_name, "Hsapiens")

    # Set the chromosome naming style to UCSC (chr1 format)
    GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"

    # Obtain the standard chromosome length and match it to the current GRanges object
    seqlengths_hg <- GenomeInfoDb::seqlengths(bs_genome)
    valid_chroms <- dplyr::intersect(GenomeInfoDb::seqlevels(gr), names(seqlengths_hg))
    GenomeInfoDb::seqlengths(gr) <- seqlengths_hg[valid_chroms]

    # Generate blank areas (CN=2) and merge them
    gaps_gr <- GenomicRanges::gaps(gr)
    GenomicRanges::mcols(gaps_gr)$cn <- 2
    full_gr <- sort(c(gr, gaps_gr))

    # Filter out whole chromosomal rows with strand information of +/-
    full_gr[BiocGenerics::strand(full_gr) == "*"]
  }
}
