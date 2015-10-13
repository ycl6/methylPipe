\name{extractBinGRanges}
\alias{extractBinGRanges}
\title{Extract genomic ranges for a given bin}
\description{
  For genomic ranges with N bins, extract the Genomic ranges for a given bin.
}
\usage{
extractBinGRanges(GenoRanges, bin, nbins)
}
\arguments{
  \item{GenoRanges}{An object of class \link{GRanges}}
  \item{bin}{numeric; the bin corresponding to which data has to be extracted}
  \item{nbins}{numeric; the number of bins in which genomic regions are divided}
}
\value{
  A \link{GRanges} Object
}
\author{Mattia Pelizzola}
\seealso{
  \code{\link{mapBSdata2GRangesBin}}
}
\examples{
gr_file <- system.file('extdata', 'GR_chr20.Rdata', package='methylPipe')
load(gr_file)
extractBinGRanges(GR_chr20, 2, 5)
}

