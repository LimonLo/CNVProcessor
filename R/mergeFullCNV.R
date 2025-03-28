#' @title Merge CNV region of samples
#' @description
#' Segment the CNV regions of multiple samples according to their loci, while populating the uncovered genomic fragments
#'
#' @param sample_list A list of GRanges objects
#' @param genome_version Genome version to use for chromosome lengths. Options are `"hg19"` or `"hg38"`.Ignored if `seqlengths` is provided.
#' @param default_cn numeric value(The copy number of the uncovered fragment).  Default: `2`.
#' @param seqlengths Optional named numeric vector of chromosome lengths (e.g., `c(chr1=249250621, ...)`).
#'               If `NULL`, chromosome lengths are automatically loaded from the specified `genome_version`.
#'
#' @return A GRanges objects containing CNV information from multiple samples in the same region (but complete chromosomes)
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
#' gr_merge <- mergeFullCNV(sample_list=list(gr1=gr1,gr2=gr2),genome_version ="hg19")
#' head(gr_merge)
mergeFullCNV <- function(sample_list,
                                    genome_version = c("hg19", "hg38"),
                                    default_cn = 2,
                                    seqlengths = NULL) {
  #Parameter validation
  stopifnot(
    is.list(sample_list) && length(sample_list) > 0,
    all(vapply(sample_list, function(x) is(x, "GRanges"), logical(1))),
    all(vapply(sample_list, function(x) "cn" %in% colnames(mcols(x)), logical(1)))
  )

  #Fill the chromosomal uncovered area of each sample
  genome_version <- match.arg(genome_version)
  fill_chromosome_gaps <- function(cnv_gr) {
    # Obtain chromosome length
    if (is.null(seqlengths)) {
      pkg_name <- switch(
        genome_version,
        "hg19" = "BSgenome.Hsapiens.UCSC.hg19",
        "hg38" = "BSgenome.Hsapiens.UCSC.hg38"
      )
      if (!requireNamespace(pkg_name, quietly = TRUE)) {
        stop("Please install ", pkg_name)
      }
      bs_genome <- getExportedValue(pkg_name, "Hsapiens")
      seqlengths <- GenomeInfoDb::seqlengths(bs_genome)
    }

    # Set the chromosome length and naming style
    GenomeInfoDb::seqlevelsStyle(cnv_gr) <- "UCSC"
    valid_chroms <- dplyr::intersect(GenomeInfoDb::seqlevels(cnv_gr), names(seqlengths))
    GenomeInfoDb::seqlengths(cnv_gr) <- seqlengths[valid_chroms]

    # Generate Fill Area (Default CN = default_cn)
    gaps_gr <- GenomicRanges::gaps(cnv_gr)
    GenomicRanges::mcols(gaps_gr)$cn <- default_cn
    full_gr <- sort(c(cnv_gr, gaps_gr))

    # Filter out rows of complete chromosomes with strand information of /-
    full_gr[BiocGenerics::strand(full_gr) == "*"]
  }

  # Fill the chromosomal gaps for each sample
  filled_samples <- lapply(sample_list,fill_chromosome_gaps)
  #
  combined_gr <- filled_samples[[1]]
  for (i in seq_along(filled_samples)[-1]) {
    combined_gr <- c(combined_gr, filled_samples[[i]])
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
  mcols_data <- lapply(sample_list, function(sample) {
    fill_cn_values(combined_gr, sample, default_cn)
  })
  GenomicRanges::mcols(combined_gr) <- S4Vectors::DataFrame(setNames(mcols_data, names(sample_list)))

  # Returns the sorted results
  sort(combined_gr)

}
