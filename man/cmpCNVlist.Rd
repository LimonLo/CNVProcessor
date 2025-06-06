% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cmpCNVlist.R
\name{cmpCNVlist}
\alias{cmpCNVlist}
\title{Compare CNVs Between multiple samples and one control}
\usage{
cmpCNVlist(Treat, Control, default_cn = 2)
}
\arguments{
\item{Treat}{A list of GRanges objects}

\item{Control}{A GRanges object}

\item{default_cn}{numeric value(The copy number of the uncovered fragment in sample B for sample A ). Default: \code{2}.}
}
\value{
A list of GRanges objects (each element is the result of a comparison)
}
\description{
CNV comparison between multiple samples and one control, resulting in a list of two samples compared.
It can be thought of as multiple uses of the cmpCNV function.
}
\examples{
gr1_file <- system.file("extdata", "test1.txt", package = "CNVProcessor")
if (file.exists(gr1_file)) {
  gr1 <- loadCNVtoGranges(gr1_file, start = "start", end = "end", cn = "cn", seqnames = "seqnames")
}
gr2_file <- system.file("extdata", "test2.txt", package = "CNVProcessor")
if (file.exists(gr2_file)) {
  gr2 <- loadCNVtoGranges(gr2_file, start = "start", end = "end", cn = "cn", seqnames = "seqnames")
}
gr3_file <- system.file("extdata", "test3.txt", package = "CNVProcessor")
if (file.exists(gr3_file)) {
  gr3 <- loadCNVtoGranges(gr3_file, start = "start", end = "end", cn = "cn", seqnames = "seqnames")
}
cmpCNVlist_result <- cmpCNVlist(Treat=list(gr1=gr1,gr2=gr2),Control=gr3)
head(cmpCNVlist_result)
}
