#' GRanges list into a data frame
#' @description
#' Consolidate lists with multiple GRanges objects into a single data frame
#'
#' @param list_gr A list of GRanges objects
#'
#' @return A data frame
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
#' df <- CNVlist2df(cmpCNVlist_result)
#' head(df)
#'
CNVlist2df <- function(list_gr){
  list_dataframe <- list()
  for(i in 1:length(list_gr)){
    list_dataframe[[i]] <- as.data.frame(list_gr[[i]])
  }
  combined_results <- dplyr::bind_rows(list_dataframe)
}
