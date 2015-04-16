\name{getCposDensity}
\alias{getCposDensity}
\title{Determines the density of genomic Cxx positions for a series of
  genomic regions}
\description{
  After having used getCpos (or getCposChr), getCposDensity determines
  the density of Cxx sites for each bin of each genomic region.
}
\usage{
getCposDensity(GenoRanges, Cpos, nbins)
}
\arguments{
  \item{GenoRanges}{an object of class \link{GRanges} used to generate
    the Cpos list}
  \item{Cpos}{list returned by getCpos or getCposChr methods}
  \item{nbins}{numeric; the number of bins each region of genomic
    regions is divided}
}
\value{
  Returns a list with the number of Cxx sites per bp of bin size for
  each region of the \link{GRanges}.
}
\author{Mattia Pelizzola}
\seealso{
  \code{\link{getCpos}}, \code{\link{profileDNAmetBin}}
}
\examples{
require(BSgenome.Hsapiens.UCSC.hg18)
gr_file <- system.file('extdata', 'GR_chr20.Rdata', package='methylPipe')
load(gr_file)
resC <- getCposChr(GenoRanges=GR_chr20, seqContext='CG', chrseq=unmasked(Hsapiens[['chr20']]), nbins=3)
resd <- getCposDensity(GenoRanges=GR_chr20, Cpos= resC, nbins=3)
}

