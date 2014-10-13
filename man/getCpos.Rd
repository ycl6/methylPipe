\name{getCpos}
\alias{getCpos}
\alias{getCposChr}
\title{Get genomic Cxx positons for a series of genomic regions}
\description{
  getCpos retrieves genomic Cxx positions, possible target of DNA
  methylation for a series of genomic regions (and bins thereof) and a
  given organism. getCposChr is a Helper function which performs the same
  task for any given DNAString sequence and is not intended for the user
  to call directly.
}
\usage{
getCpos(GenoRanges, seqContext='all', nbins, org)
getCposChr(GenoRanges, seqContext, chrseq, nbins)
}
\arguments{
  \item{GenoRanges}{An object of class \link{GRanges}}
  \item{seqContext}{character; one of: all, CG, CHG or CHH}
  \item{org}{an object of class BSgenome; typically the genome
    sequences of a given organism}
  \item{chrseq}{an object of class \link{DNAString}; typically a
    chromosome sequence of a given organism}
  \item{nbins}{numeric; the number of bins each region of genomic
    regions is divided}
}
\value{
  A list is returned with the position of the Cxx in the \link{GRanges}
  regions. The length of the list is equal to the length of the
  \link{GRanges}. For each list item a list with length equal to the
  number of bins of the \link{GRanges} is returned. For each bin the
  position of the Cxx relative to the genomic coordinates of that bin is
  returned.
}
\author{Mattia Pelizzola}
\seealso{
  \code{\link{getCposDensity}}, \code{\link{profileDNAmetBin}}
}
\examples{
require(BSgenome.Hsapiens.UCSC.hg18)
gr_file <- system.file('extdata', 'GR_chr20.Rdata', package='methylPipe')
load(gr_file)
res <- getCpos(GR_chr20, seqContext='CG', nbins=1, org=Hsapiens)
res[[1]]
}

