\name{GEcollection-class}
\docType{class}
\alias{GEcollection-class}
\alias{GEcollection}
\alias{show,GEcollection-method}
\alias{chr,GEcollection-method}
\alias{binC,GEcollection-method}
\alias{binmC,GEcollection-method}
\alias{binrC,GEcollection-method}
\alias{nbins,GEcollection-method}
\alias{length,GEcollection-method}
\alias{binscore,GEcollection-method}
\alias{binscore<-,GEcollection-method}
\alias{chr}
\alias{binC}
\alias{binmC}
\alias{binrC}
\alias{nbins}
\alias{binscore}
\alias{length}
\alias{binscore<-}

\title{Class "GEcollection"}
\description{ This class is used in the methylPipe library to define and
  manipulate a set of genomic regions and the associated DNA methylation
  patterns}
\section{Objects from the Class}{
  This class is an extension of \link{SummarizedExperiment} class of
  GenomicRanges R package. Objects can be created using the function
  \code{profileDNAmetBin} which determines the absolute and relative
  methylation level by filling the binC, binmC and binrC slots. The
  assays slot \link{SummarizedExperiment} class here consists of four
  slots.
\itemize{
\item binC: each genomic region is divided in one or more bins and for each
  bin the density (per bp) of potential methylation sites is determined.
\item binmC:
  each genomic region is divided in one or more bins and for each bin the
  density (per bp) of methylation events is determined.
\item binrC: each genomic region is
  divided in one or more bins and for each bin the relative mC/C content is
  determined.
\item binscore: each genomic region is divided in one or more bins and
  scores can be assigned to them. In particular, it can be convenient for
storing reads count for each bin of each genomic region.
}
 The minimal set of data
  to create a \link{GEcollection} object is a set of genomic regions to be provided as a
  \link{GRanges} object and a dataset of class \link{BSdata}.
}


\section{Methods}{
  \describe{
    \item{chr}{\code{signature(object = "GEcollection")}: extracts the
      chr assignment of the genomic regions}
    \item{Strand}{\code{signature(object = "GEcollection")}: extracts
      the strand assignment of the genomic regions}
    \item{binC}{\code{signature(object = "GEcollection")}: extracts the
      binC matrix}
    \item{binmC}{\code{signature(object = "GEcollection")}: extracts the
      binmC matrix}
    \item{binrC}{\code{signature(object = "GEcollection")}: extracts the
      binrC matrix}
    \item{binscore}{\code{signature(object = "GEcollection")}: extracts
      the binscore matrix}
    \item{binscore<-}{\code{signature(object = "GEcollection")}:
      replaces the binscore matrix}
    \item{nbins}{\code{signature(object = "GEcollection")}: returns the
      number of bins}
  }
}
\author{Kamal Kishore}
\examples{
require(BSgenome.Hsapiens.UCSC.hg18)
uncov_GR <- GRanges(Rle('chr20'), IRanges(c(14350,69251,84185), c(18349,73250,88184)))
H1data <- system.file('extdata', 'H1_chr20_CG_10k_tabix_out.txt.gz', package='methylPipe')
H1.db <- BSdata(file=H1data, uncov=uncov_GR, org=Hsapiens)
gr_file <- system.file('extdata', 'GR_chr20.Rdata', package='methylPipe')
load(gr_file)
gec.H1 <- profileDNAmetBin(GenoRanges=GR_chr20, Sample=H1.db, mcCLASS='mCG')
gec.H1
}
\keyword{classes}

