#' @title Gene annotation for CNV
#' @description
#' Gene annotation for CNV
#'
#' @param gr A GRanges object
#' @param genome_version Genome version to use for chromosome lengths. Options are `"hg19"` or `"hg38"`. Ignored if `seqlengths` is provided.
#' @param gene_type The type of annotated gene
#' @param multiVals Equated to \code{\link[AnnotationDbi]{mapIds}}
#'
#' @return A GRanges object with gene annotation
#' @export
#'
#' @examples
#' gr1_file <- system.file("extdata", "test1.txt", package = "CNVProcessor")
#' if (file.exists(gr1_file)) {
#'   gr1 <- loadCNVtoGranges(gr1_file, start = "start", end = "end", cn = "cn", seqnames = "seqnames")
#' }
#' gr1_anno <- annoCNV(gr=gr1,genome_version = "hg19")
annoCNV <- function(gr, genome_version = c("hg19", "hg38"),
                    gene_type = "SYMBOL",
                    multiVals = "list") {
  genome_version <- match.arg(genome_version)
  stopifnot(
    is(gr, "GRanges"),
    gene_type %in% AnnotationDbi::columns(org.Hs.eg.db::org.Hs.eg.db)  # Verify the validity of the gene ID type
  )

  #load TxDb package -------------------#
  txdb_pkg <- switch(
    genome_version,
    "hg19" = "TxDb.Hsapiens.UCSC.hg19.knownGene",
    "hg38" = "TxDb.Hsapiens.UCSC.hg38.knownGene"
  )

  if (!requireNamespace(txdb_pkg, quietly = TRUE)) {
    stop("Please install", txdb_pkg)
  }
  txdb <- base::getExportedValue(txdb_pkg, txdb_pkg)

  #Obtaining Genomic Information
  genes_gr <- GenomicFeatures::genes(txdb)

  #annotation
  annotated_gr <- gr %>%
    plyranges::join_overlap_left(genes_gr) %>%
    as.data.frame(row.names = NULL) %>%
    dplyr::mutate(
      gene_symbol = AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = as.character(gene_id),
        column = gene_type,
        keytype = "ENTREZID",
        multiVals = multiVals
      )
    ) %>% tidyr::unnest(gene_symbol)

}
