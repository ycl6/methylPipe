\name{methstats}
\alias{methstats}
\alias{methstats,methylPipe,BSdataSet}
\alias{methstats,methylPipe,BSdataSet-method}
\alias{methstats-methods}
\alias{methstats,BSdataSet-method}
\title{Exploratory statistics of samples in BSdataSet object}
\description{
  Exploratory methylation statistics of samples in BSdataSet object.
}
\usage{
\S4method{methstats}{methylPipe,BSdataSet}(object, chrom, mcClass='mCG', minC=1, coverage=1, pval=1, Nproc=1)
}
\arguments{
  \item{object}{An object of class \link{BSdataSet}}
  \item{chrom}{character; either NULL or an object of class character}
  \item{mcClass}{character; the mC sequence context to be considered,
    one of all, mCG, mCHG or mCHH}
  \item{minC}{numeric; the minumum number of reads with C (DNA-methylation
    events) at a given cytosine genomic position}
  \item{coverage}{numeric; the minimum number of total reads at a given
    cytosine genomic position}
  \item{pval}{numeric; the minimum binomial pValue for an mC call at a
    given cytosine genomic position}
  \item{Nproc}{numeric; the number of processors to use, one chromosome
    is ran for each processor}
}
\details{The function provides basic statistical methods which computes descriptive statistics, 
correlation matrix and clustering of samples within the BSdataSet. 
}
\value{
  A list with slots named descriptive_stats and correlation_mat.
}
\author{Kamal Kishore}
\examples{
require(BSgenome.Hsapiens.UCSC.hg18)
uncov_GR <- GRanges('chr20', IRanges(14350, 18349))
H1data <- system.file('extdata', 'H1_chr20_CG_10k_tabix_out.txt.gz', package='methylPipe')
H1.db <- BSdata(file=H1data, uncov=uncov_GR, org=Hsapiens)
IMR90data <- system.file('extdata', 'IMR90_chr20_CG_10k_tabix_out.txt.gz', package='methylPipe')
IMR90.db <- BSdata(file=IMR90data, uncov=uncov_GR, org=Hsapiens)
H1.IMR90.set <- BSdataSet(org=Hsapiens, group=c("C","C","E","E"), IMR_1=IMR90.db, 
IMR_2=IMR90.db, H1_1=H1.db,H1_2=H1.db)
stats_res <- methstats(H1.IMR90.set,chrom="chr20",mcClass='mCG', Nproc=1)
stats_res
}

