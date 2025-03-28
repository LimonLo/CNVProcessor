#' @title Data Input
#' @description
#' Load a table file with chromosome number, start position, end position, and copy number and convert it into a GRanges object
#' @param cnvfile The path where the file is located and the name of the file
#' @param start The column name of the start position
#' @param end The column name of the start position
#' @param cn The column name of the copy number
#' @param seqnames The column name of the chromosome
#' @param header a logical value indicating whether the file contains the names of the variables as its first line.Equated \code{\link[utils]{read.table}}
#' @param sep the field separator character. Values on each line of the file are separated by this character. Equated \code{\link[utils]{read.table}}
#' @param ...
#'
#' @return A GRanges object with "cn"(copy number) as metadata
#' @export
#'
#' @examples
#'  # Load sample data
#' cnv_file <- system.file("extdata", "test1.txt", package = "CNVProcessor")
#' if (file.exists(cnv_file)) {
#'   gr <- loadCNVtoGranges(cnv_file, start = "start", end = "end", cn = "cn", seqnames = "seqnames")
#' }
#' head(gr)
loadCNVtoGranges <- function(
    cnvfile,                # file path
    start,                # start position colname
    end,                  # end position colname
    cn,                   # copy number colname
    seqnames,             # chromosome colname
    header = TRUE,            # parameter of read.table
    sep = "\t",               # parameter of read.table
    ...                       # parameter of read.table
) {
  # read.table
  data <- utils::read.table(
    file = cnvfile,
    header = header,
    sep = sep,
    stringsAsFactors = FALSE,
    ...
  )

  # Check if the necessary columns exist
  required_cols <- c(start, end, cn, seqnames)
  missing_cols <- dplyr::setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("The column does not exist: ", paste(missing_cols, collapse = ", "))
  }

  # Processing chromosome columns: Automatically add chr prefixes
  seqnames <- as.character(data[[seqnames]])
  seqnames <- ifelse(
    grepl("^chr", seqnames),
    seqnames,
    paste0("chr", seqnames)
  )

  # Creating IRanges Objects (Automatically Handling Numeric Conversions)
  ranges <- IRanges::IRanges(
    start = data[[start]],
    end = data[[end]]
  )

  # Create a GRanges object and add metadata
  gr <- GenomicRanges::GRanges(
    seqnames = seqnames,
    ranges = ranges,
    cn = data[[cn]]      # The copy number is used as metadata
  )

  return(gr)
}
