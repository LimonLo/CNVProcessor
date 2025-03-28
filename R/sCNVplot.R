#' @title CNV plot
#' @description
#' Draw the CNV of a single sample.
#'
#' @param data A GRanges object
#' @param genome_version Genome version to use for chromosome lengths. Options are `"hg19"` or `"hg38"`. Ignored if `seqlengths` is provided.
#' @param chroms The chromosome that needs to be plotted,Default: chr1-chr22,chrX,chrY
#' @param cn_col The colname the CN. Usually the fold Change column of the comparison results
#' @param amp_color The color of the amplified area, Default: red
#' @param nor_color The color of the normal area, Default: yellow
#' @param del_color The color of the deleted area, Default: blue
#' @param title The title of the diagram
#' @param normal A numerical value that defines the number of normal chromosomal copies. Default: 2
#' @param xlab The X-axis title
#' @param ylab The Y-axis title
#' @param show_chrom_labels A logical value that displays the chromosome tag or not
#'
#' @return Diagram of the CNV of a single sample.
#' @export
#'
#' @examples
#' cnv_file <- system.file("extdata", "test1.txt", package = "CNVProcessor")
#' if (file.exists(cnv_file)) {
#'   gr <- loadCNVtoGranges(cnv_file, start = "start", end = "end", cn = "cn", seqnames = "seqnames")
#' }
#' p <- sCNVplot(data=gr,genome_version = "hg19",title = "",cn_col = "cn",normal = 2,chroms = paste0("chr", 1:5),ylab = "CNV",show_chrom_labels = T)
sCNVplot <- function(data,
                          genome_version = c("hg19", "hg38"),
                          chroms = paste0("chr", c(1:22, "X", "Y")),
                          cn_col = "cn",
                          amp_color = "red",
                          nor_color = "yellow",
                          del_color = "blue",
                          title = "",
                          normal = 2,
                          xlab = "Chromosome",
                          ylab = "CNV Difference (Treated - Control)",
                          show_chrom_labels = TRUE) {

  # Parameter Verification and Data Preparation
  genome_version <- match.arg(genome_version)
  if (is(data, "GRanges")) data <- as.data.frame(data)

  # Verify that the column must exist
  required_cols <- c("seqnames", "start", "end", cn_col)
  if (!all(required_cols %in% colnames(data))) {
    stop("The input data must contain columns:", paste(required_cols, collapse = ", "))
  }

  #Obtaining Genomic Information
  bs_pkg <- switch(
    genome_version,
    "hg19" = "BSgenome.Hsapiens.UCSC.hg19",
    "hg38" = "BSgenome.Hsapiens.UCSC.hg38"
  )
  if (!requireNamespace(bs_pkg, quietly = TRUE)) {
    stop("Please install", bs_pkg)
  }
  bs_genome <- getExportedValue(bs_pkg, "Hsapiens")
  all_seqlengths <- GenomeInfoDb::seqlengths(bs_genome)

  # Verify chromosomal validity
  invalid_chroms <- dplyr::setdiff(chroms, names(all_seqlengths))
  if (length(invalid_chroms) > 0) {
    stop("The following chromosomes are in ", genome_version, "does not exist in:", paste(invalid_chroms, collapse = ", "))
  }
  seqlengths <- all_seqlengths[chroms]

  #Calculate Global Coordinates
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
      FC = .data[[cn_col]],
      # Create a categorical column for legend
      region_type = dplyr::case_when(
        FC > normal ~ "Amplification",
        FC < normal ~ "Deletion",
        TRUE ~ "Normal"
      )
    )

  #Calculate chromosome separator coordinates
  chrom_ends <- cum_offset %>%
    dplyr::mutate(
      chrom_end = seqlengths[seqnames] + offset,
      mid_pos = offset + seqlengths[seqnames]/2
    )

  # plot
  p <- ggplot2::ggplot(plot_data) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = global_start, xmax = global_end,
          ymin = 0, ymax = 1,
          fill = region_type),
      alpha = 0.7
    ) +
    # Chromosome divider line (runs through all samples)
    ggplot2::geom_vline(
      data = chrom_ends,
      ggplot2::aes(xintercept = chrom_end),
      color = "gray40", linetype = "dashed", linewidth = 0.5
    ) +
    ggplot2::scale_fill_manual(
      name = NULL,
      values = c("Amplification" = amp_color,
                 "Deletion" = del_color,
                 "Normal" = nor_color)
    ) +
    # Set chromosome labels and scales
    ggplot2::scale_x_continuous(
      expand = c(0, 0),
      breaks = chrom_ends$mid_pos,  # The middle position of the chromosome serves as a scale
      labels = if(show_chrom_labels) sub("chr", "", chroms) else NULL  # tag
    ) +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, vjust = 0.5),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1),
      legend.position = "top",
      panel.spacing = ggplot2::unit(0, "lines"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA)
    )
  return(p)
}
