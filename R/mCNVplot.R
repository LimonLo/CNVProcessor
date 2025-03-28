#' @title CNV plot
#' @description
#' Draw the CNV of multiple samples.
#'
#' @param data_list A list of GRanges object
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
#' @param sample_names A string vector containing the sample names
#'
#' @return Diagram of the CNV of multiple samples.
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
#' p <- mCNVplot(data_list=list(gr1=gr1,gr2=gr2,gr3=gr3),genome_version = "hg19",cn_col = "cn")
mCNVplot <- function(data_list,
                         genome_version = c("hg19", "hg38"),
                         chroms = paste0("chr", c(1:22, "X", "Y")),
                         cn_col = "cn",
                         amp_color = "red",
                         nor_color = "yellow",
                         del_color = "blue",
                         title = "",
                         normal = 2,
                         xlab = "Chromosome",
                         ylab = "CNV Difference",
                         show_chrom_labels = TRUE,
                         sample_names = NULL) {
  #
  genome_version <- match.arg(genome_version)
  stopifnot(
    is.list(data_list),
    all(sapply(data_list, function(x) is(x, "GRanges") || is.data.frame(x))),
    length(data_list) > 0
  )

  #
  if (is.null(sample_names)) {
    sample_names <- paste0("Sample", seq_along(data_list))
  } else if (length(sample_names) != length(data_list)) {
    stop("The number of sample names must match the length of the input list")
  }

  #Merge all sample data
  combined_data <- lapply(seq_along(data_list), function(i) {
    sample_data <- data_list[[i]]
    if (is(sample_data, "GRanges")) sample_data <- as.data.frame(sample_data)

    # 添加样本ID列
    sample_data$sample_id <- sample_names[i]
    sample_data
  }) %>% dplyr::bind_rows()

  #Obtaining Genomic Information
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

  # Calculate Global Coordinates
  cum_offset <- data.frame(
    seqnames = chroms,
    offset = cumsum(c(0, head(all_seqlengths[chroms], -1))),
    stringsAsFactors = FALSE
  )

  plot_data <- combined_data %>%
    dplyr::filter(seqnames %in% chroms) %>%
    dplyr::left_join(cum_offset, by = "seqnames") %>%
    dplyr::mutate(
      global_start = start + offset,
      global_end = end + offset,
      FC = .data[[cn_col]],
      region_type = dplyr::case_when(
        FC > normal ~ "Amplification",
        FC < normal ~ "Deletion",
        TRUE ~ "Normal"
      )
    )

  #  Calculate chromosome separator coordinates
  chrom_ends <- cum_offset %>%
    dplyr::mutate(
      chrom_end = all_seqlengths[seqnames] + offset,
      mid_pos = offset + all_seqlengths[seqnames]/2
    )

  #Plot
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
    # Set facets and axes
    ggplot2::facet_grid(
      rows = vars(sample_id),
      switch = "y",
      scales = "free_y"
    ) +
    # Set chromosome labels and scales
    ggplot2::scale_x_continuous(
      expand = c(0, 0),
      breaks = chrom_ends$mid_pos,  #
      labels = if(show_chrom_labels) sub("chr", "", chroms) else NULL  #
    ) +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, vjust = 0.5),  #
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1),
      legend.position = "top",
      panel.spacing = ggplot2::unit(0, "lines")
    )

  return(p)
}
