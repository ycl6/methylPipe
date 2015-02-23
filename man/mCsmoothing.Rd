\name{mCsmoothing}
\alias{mCsmoothing}
\alias{mCsmoothing,methylPipe,BSdata}
\alias{mCsmoothing,methylPipe,BSdata-method}
\alias{mCsmoothing-methods}
\alias{mCsmoothing,BSdata-method}

\title{Smoothing and plotting methylation data}
\description{
  Smoothing and plotting methylation data, even chromosome wide.
}
\usage{
\S4method{mCsmoothing}{methylPipe,BSdata}(Object, refgr, Scorefun='sum', Nbins=20,
Context="CG", plot=TRUE)
}
\arguments{
  \item{Object}{An object of class \link{BSdata}}
  \item{refgr}{GRanges; Genomic Ranges to plot the data}
  \item{Scorefun}{character; either sum or mean for smoothing}
  \item{Nbins}{numeric; the number of interval each range is divided}
  \item{Context}{character; either all or a combination of CG, CHG, and CHH}
  \item{plot}{logical; whether the smoothed profile has to be plotted}
}
\details{The sum or the mean methylation level is determined on each
  window of size Binsize and smoothed with the smooth.spline function.}
\value{
  A list with three components: pos (the left most point of each
  window), score (either the sum or the mean methylation levels),
  smoothed (the smoothed methylation levels).
}
\author{Mattia Pelizzola}
\examples{
require(BSgenome.Hsapiens.UCSC.hg18)
uncov_GR <- GRanges(Rle('chr20'), IRanges(c(14350,69251,84185), c(18349,73250,88184)))
H1data <- system.file('extdata', 'H1_chr20_CG_10k_tabix_out.txt.gz', package='methylPipe')
H1.db <- BSdata(file=H1data, uncov=uncov_GR, org=Hsapiens)
gr <- GRanges("chr20",IRanges(1,5e5))
sres <- mCsmoothing(H1.db, gr, Scorefun='sum', Nbins=50, Context="CG", plot=TRUE)
}

