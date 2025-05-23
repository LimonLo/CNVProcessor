% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annoCNV.R
\name{annoCNV}
\alias{annoCNV}
\title{Gene annotation for CNV}
\usage{
annoCNV(
  gr,
  genome_version = c("hg19", "hg38"),
  gene_type = "SYMBOL",
  multiVals = "list"
)
}
\arguments{
\item{gr}{A GRanges object}

\item{genome_version}{Genome version to use for chromosome lengths. Options are \code{"hg19"} or \code{"hg38"}. Ignored if \code{seqlengths} is provided.}

\item{gene_type}{The type of annotated gene}

\item{multiVals}{Equated to \code{\link[AnnotationDbi]{mapIds}}}
}
\value{
A GRanges object with gene annotation
}
\description{
Gene annotation for CNV
}
\examples{
gr1_file <- system.file("extdata", "test1.txt", package = "CNVProcessor")
if (file.exists(gr1_file)) {
  gr1 <- loadCNVtoGranges(gr1_file, start = "start", end = "end", cn = "cn", seqnames = "seqnames")
}
gr1_anno <- annoCNV(gr=gr1,genome_version = "hg19")
}
