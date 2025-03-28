#' @title Compare CNVs Between Two Groups
#' @description
#' Performs statistical comparison of copy number (or other genomic features) between control and treatment groups
#'     across genomic regions. Identifies regions with significant amplification/deletion based on statistical thresholds.
#'
#' @param gr A GRanges object with a record of the CN values of two groups in the same chromosomal region (available via mergeFullCNV or mergeRawCNV)
#' @param control Character vector
#' @param treat Character vector
#' @param method Statistical test method. Options: `"t.test"` or `"wilcox"`. Default: `"t.test"`.
#' @param p Adjusted p-value (FDR) cutoff for significance. Default: `0.05`.
#' @param FC Absolute fold change (FC) cutoff to define amplification/deletion. Default: `1`.
#'
#' @return A GRanges object
#' @export
#'
#' @examples
#' gr1_file <- system.file("extdata", "test1.txt", package = "CNVProcessor")
#' if (file.exists(gr1_file)) {
#'   gr1 <- loadCNVtoGranges(gr1_file, start = "start", end = "end", cn = "cn", seqnames = "seqnames")
#' }
#' gr2_file <- system.file("extdata", "test2.txt", package = "CNVProcessor")
#' if (file.exists(gr2_file)) {
#'   gr2 <- loadCNVtoGranges(gr1_file, start = "start", end = "end", cn = "cn", seqnames = "seqnames")
#' }
#' gr3_file <- system.file("extdata", "test3.txt", package = "CNVProcessor")
#' if (file.exists(gr3_file)) {
#'   gr3 <- loadCNVtoGranges(gr2_file, start = "start", end = "end", cn = "cn", seqnames = "seqnames")
#' }
#' gr4_file <- system.file("extdata", "test4.txt", package = "CNVProcessor")
#' if (file.exists(gr4_file)) {
#'   gr4 <- loadCNVtoGranges(gr3_file, start = "start", end = "end", cn = "cn", seqnames = "seqnames")
#' }
#' gr_merge <- mergeFullCNV(sample_list=list(gr1=gr1,gr2=gr2,gr3=gr3,gr4=gr4),genome_version ="hg19")
#' cmpGroups_result <- cmpGroups(gr=gr_merge,control=c("gr1","gr2"),treat=c("gr3","gr4"),method="t.test")
#' head(cmpGroups_result)
cmpGroups <- function(gr, control, treat,
                           method = c("t.test","wilcox"), p= 0.05,
                           FC = 1) {
  # Extract the data matrix
  control_mat <- as.matrix(mcols(gr)[, control, drop = FALSE])
  treat_mat <- as.matrix(mcols(gr)[, treat, drop = FALSE])

  # Calculate statistical indicators
  GenomicRanges::mcols(gr)$mean_control <- base::rowMeans(control_mat)
  GenomicRanges::mcols(gr)$mean_treat <- base::rowMeans(treat_mat)
  GenomicRanges::mcols(gr)$FC <- gr$mean_treat-gr$mean_control

  # Calculate p-values
  pvals <- vapply(1:nrow(control_mat), function(i) {
    if (method == "t.test") {
      tryCatch(t.test(treat_mat[i, ], control_mat[i, ])$p.value,
               error = function(e) NA)
    } else if (method == "wilcox") {
      suppressWarnings(wilcox.test(treat_mat[i, ], control_mat[i, ])$p.value)
    } else NA
  }, numeric(1))

  GenomicRanges::mcols(gr)$p_value <- pvals
  GenomicRanges::mcols(gr)$fdr <- p.adjust(pvals, method = "BH")

  # Define the type of difference
  GenomicRanges::mcols(gr)$change <- ifelse(
    GenomicRanges::mcols(gr)$fdr < p &
      GenomicRanges::mcols(gr)$FC >= FC,
    "amplification",
    ifelse(
      GenomicRanges::mcols(gr)$fdr < p &
        GenomicRanges::mcols(gr)$FC <= -FC,
      "deletion",
      "no_change"
    )
  )
  keep_cols <- c("mean_control", "mean_treat", "FC", "p_value", "fdr", "change")
  GenomicRanges::mcols(gr) <- GenomicRanges::mcols(gr)[, keep_cols]

  return(gr)
}
