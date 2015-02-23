\name{profileDNAmetBin}
\alias{profileDNAmetBin}
\alias{profileDNAmetBinParallel}
\title{Profile DNA methylation data for a set of genomic regions}
\description{
  Profile the absolute and relative density of mC sites for each bin of
  each genomic region of a \link{GEcollection} object.
}
\usage{
profileDNAmetBin(GenoRanges, Sample, mcCLASS="mCG",
mC=1, depthThr=0, mCpv=1, minCoverage=0.75, nbins = 2)
profileDNAmetBinParallel(GenoRanges, Sample, mcCLASS="mCG", mC=1,
depthThr=0, mCpv=1, minCoverage=0.75, Nproc=1, nbins = 2)
}
\arguments{
  \item{GenoRanges}{an object of class \link{GRanges}}
  \item{Sample}{an object of class \link{BSdata}}
  \item{mcCLASS}{character; one of: mCG, mCHG, mCHH}
  \item{mC}{numeric; the minumum number of reads with C (DNA-methylation
    events) at a given cytosine genomic position}
  \item{depthThr}{numeric; the minimum number of total reads at a given
    cytosine genomic position}
  \item{mCpv}{numeric; the minimum binomial pValue for an mC call at a
    given cytosine genomic position}
  \item{minCoverage}{numeric between 0 and 1; the minimum coverage of
    for the genomic region to be profiled; currently ignored}
  \item{Nproc}{numeric; the number of processor for parallelization}
  \item{nbins}{numeric; the number of bins each genomic region is divided}
}
\details{
  For each bin of each genomic region of a \link{GRanges} object, the
  absolute and relative density of mC sites is determined and an object
  of class \link{GEcollection} is created.
}
\value{
  An object of class \link{GRanges} from which an object of class
  \link{GEcollection} is created with the binC, binmC and binrC data
  slots been filled with density matrices.
  These matrices report the density of mC sites in the sequence context
  specified by mcCLASS. They are counted for each bin in each genomic
  region and their count is divided by the bin size in bp.
  The binC data slot is filled with the density of all possible
  methylation sites in the specified sequence context.
  The binmC data slot is filled with the density of mC sites in the
  specified sequence context for the considered sample.
  The binrC data slot is filled with the ratio of binC and binmC
  matrices, representing the relative methylation for each bin in each genomic region.
}
\author{Mattia Pelizzola, Kamal Kishore}
\examples{
require(BSgenome.Hsapiens.UCSC.hg18)
H1data <- system.file('extdata', 'H1_chr20_CG_10k_tabix_out.txt.gz', package='methylPipe')
uncov_GR <- GRanges(Rle('chr20'), IRanges(c(14350,69251,84185), c(18349,73250,88184)))
H1.db <- BSdata(file=H1data, uncov= uncov_GR, org=Hsapiens)
gr_file <- system.file('extdata', 'GR_chr20.Rdata', package='methylPipe')
load(gr_file)
gec.H1 <- profileDNAmetBin(GenoRanges=GR_chr20, Sample=H1.db, mcCLASS='mCG', nbins=2)
head(binmC(gec.H1))
}

