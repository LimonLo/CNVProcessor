#' @title CNV plot
#' @description
#' It is mainly used to show the comparison results of samples CNVs.
#'
#' @param data A GRanges object
#' @param genome_version Genome version to use for chromosome lengths. Options are `"hg19"` or `"hg38"`. Ignored if `seqlengths` is provided.
#' @param chroms The chromosome that needs to be plotted,Default: chr1-chr22,chrX,chrY
#' @param fc_col The colname the CN. Usually the fold Change column of the comparison results
#' @param amp_color The color of the amplified area, Default: red
#' @param del_color The color of the deleted area, Default: blue
#' @param title The title of the diagram
#' @param normal A numerical value that defines the number of normal chromosomal copies. Default: 0 (the difference is 0 in comparison)
#' @param show_chrom_labels A logical value that displays the chromosome tag or not
#' @param xlab The X-axis title
#' @param ylab The Y-axis title
#'
#' @return Diagram of chromosomal differences
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
#' gr_merge <- cmpCNV(sample1=gr1,sample2=gr2)
#' p <- cmpCNVplot(data=gr_merge,genome_version = "hg19",fc_col = "FC")
cmpCNVplot <- function(data,
                           genome_version = c("hg19", "hg38"),
                           chroms = paste0("chr", c(1:22, "X", "Y")),
                           fc_col = "FC",
                           amp_color = "red",
                           del_color = "blue",
                           title = "",
                           normal=0,show_chrom_labels=TRUE,
                           xlab = "Chromosome", ylab = "CNV Difference (Treated - Control)") {
  # Convert data into data frames and merge global coordinates
  if (is(data, "GRanges")) data <- as.data.frame(data)
  # Parameter validation
  genome_version <- match.arg(genome_version)
  required_cols <- c("seqnames", "start", "end", fc_col)
  if (!all(required_cols %in% colnames(data))) {
    stop("The input data must contain columns: ", paste(required_cols, collapse = ", "))
  }

  # Obtain genomic chromosome length
  bs_pkg <- switch(
    genome_version,
    "hg19" = "BSgenome.Hsapiens.UCSC.hg19",
    "hg38" = "BSgenome.Hsapiens.UCSC.hg38"
  )
  if (!requireNamespace(bs_pkg, quietly = TRUE)) {
    stop("Please install ", bs_pkg)
  }
  bs_genome <- getExportedValue(bs_pkg, "Hsapiens")
  all_seqlengths <- GenomeInfoDb::seqlengths(bs_genome)

  # Verify that the chromosome specified by the user is valid
  invalid_chroms <- dplyr::setdiff(chroms, names(all_seqlengths))
  if (length(invalid_chroms) > 0) {
    stop("The following chromosomes are in ", genome_version, " does not exist：", paste(invalid_chroms, collapse = ", "))
  }
  seqlengths <- all_seqlengths[chroms]

  # Calculate cumulative coordinate offset (for user-selected chromosomes only)
  cum_offset <- data.frame(
    seqnames = chroms,
    offset = cumsum(c(0, head(seqlengths, -1))),
    stringsAsFactors = FALSE
  )

  plot_data <- data %>%
    dplyr::filter(seqnames %in% chroms) %>%
    dplyr::left_join(cum_offset, by = "seqnames") %>%
    dplyr::mutate(
      global_start = start + offset,
      global_end = end + offset,
      FC = .data[[fc_col]]
    )

  # Calculate chromosome separator coordinates (only for user-selected chromosomes)
  chrom_ends <- cum_offset %>%
    dplyr::mutate(
      chrom_end = seqlengths[seqnames] + offset,
      mid_pos = offset + seqlengths[seqnames]/2
    )

  # plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = dplyr::filter(plot_data, FC > normal),
      ggplot2::aes(xmin = global_start, xmax = global_end, ymin = 0, ymax = FC),
      fill = amp_color, color = NA, alpha = 0.7
    ) +
    ggplot2::geom_rect(
      data = dplyr::filter(plot_data, FC < normal),
      ggplot2::aes(xmin = global_start, xmax = global_end, ymin = FC, ymax = 0),
      fill = del_color, color = NA, alpha = 0.7
    ) +
    ggplot2::geom_vline(
      data = chrom_ends,
      ggplot2::aes(xintercept = chrom_end),
      color = "gray40", linetype = "dashed", linewidth = 0.3
    ) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.1)) +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey90"),
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.border = ggplot2::element_rect(color = "black", fill = NA)
    )
  # Optional: Add chromosome tag
  if (show_chrom_labels) {
    p <- p + ggplot2::geom_text(
      data = chrom_ends,
      ggplot2::aes(x = mid_pos, y = max(abs(plot_data$FC)) * 1.1,
          label = sub("chr", "", seqnames)),
      angle = 0, size = 3, vjust = 0.5
    )
  }
  return(p)
}
