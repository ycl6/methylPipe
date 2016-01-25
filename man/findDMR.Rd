\name{findDMR}
\alias{findDMR}
\alias{findDMR,methylPipe,BSdataSet}
\alias{findDMR,methylPipe,BSdataSet-method}
\alias{findDMR-methods}
\alias{findDMR,BSdataSet-method}
\title{Identifying Differentially Methylated Regions (DMRs)}
\description{
  Identifying differentially methylated regions for pairwise or multiple
  samples comparision.
}
\usage{
\S4method{findDMR}{methylPipe,BSdataSet}(object, Nproc=NULL, ROI=NULL,
pmdGRanges=NULL, MCClass='mCG', dmrSize=10, dmrBp=1000, binsize=0,
eprop=0.3, coverage=1, Pvalue=NULL, SNPs=NULL)
}
\arguments{
  \item{object}{An object of class \link{BSdataSet}}
  \item{Nproc}{numeric; the number of processors to use, one chromosome
    is ran for each processor}
  \item{ROI}{character; either NULL or an object of class GRanges
    consisting of genomic regions of interest for which DMRs are
    identified}
  \item{pmdGRanges}{a \link{GRanges} object containing the genomic
    coordinates of Partially Methylated Domains that will be masked}
  \item{MCClass}{character; the mC sequence context to be considered,
    one of all, mCG, mCHG or mCHH}
  \item{dmrSize}{numeric; the number of consecutive mC to be
    simulataneously considered; atleast 5}
  \item{dmrBp}{numeric; the max number of bp containing the dmrSize mC}
  \item{binsize}{numeric; the size of the bin used for smoothing the
    methylation levels, useful for nonCG methylation in human}
  \item{eprop}{numeric; the max - min methylation level is determined
    for each mC, or for each bin, and only mC (or bins) with
    difference greater than eprop are considered}
  \item{coverage}{numeric; the minimum number of total reads at a given
    cytosine genomic position}
  \item{Pvalue}{numeric; to select only those mC with significant p-value}
  \item{SNPs}{GRanges; if SNPs information is provided those cytosine
    are removed from DMR computation}
}
\details{Typically for nonCG methylation in human a dmrSize of 50, a
  dmrBp of 50000 and a binsize of 1000 are used. For CpG methylation in
  human and both CpG and nonCG methylation in plants the default
  settings are usually fine. Partially Methylated Domains are usually
  found in differentiated cells and can constitute up to one third of
  the genome (Lister R et al, Nature 2009). Usually DMRs are not
  selected within those regions especially when comparing differentiated
  and pluripotent cells. Eprop is used to speed up the
  analysis. According to the number of samples different test are used
  to compare the methylation levels (percentage of methylated reads for
  each mC). In case of two samples the non parametric wilcoxon test is
  used. In case of more than two samples the kruskal wallis non
  parametric testis used. Check \code{\link{consolidateDMRs}} to further
  process and finalize the list of DMRs.
}
\value{
  A \link{GRanges} object of DMRs with the metadata slots for pValue,
  MethDiff_Perc and log2Enrichment. When two samples are compared, MethDiff_Perc is the
  diference between percentage methylation between the conditions compared. However, log2Enrichment 
  is the log2ratio between the mean for the samples. 
}
\author{Mattia Pelizzola, Kamal Kishore}
\seealso{
  \code{\link{consolidateDMRs}}
}
\examples{
require(BSgenome.Hsapiens.UCSC.hg18)
uncov_GR <- GRanges(Rle('chr20'), IRanges(c(14350,69251,84185), c(18349,73250,88184)))
H1data <- system.file('extdata', 'H1_chr20_CG_10k_tabix_out.txt.gz', package='methylPipe')
H1.db <- BSdata(file=H1data, uncov=uncov_GR, org=Hsapiens)
IMR90data <- system.file('extdata', 'IMR90_chr20_CG_10k_tabix_out.txt.gz', package='methylPipe')
IMR90.db <- BSdata(file=IMR90data, uncov=uncov_GR, org=Hsapiens)
H1.IMR90.set <- BSdataSet(org=Hsapiens, group=c("C","E"), IMR90=IMR90.db, H1=H1.db)
gr_file <- system.file('extdata', 'GR_chr20.Rdata', package='methylPipe')
load(gr_file)
DMRs <- findDMR(object= H1.IMR90.set, Nproc=1, ROI=GR_chr20, MCClass='mCG',
dmrSize=10, dmrBp=1000, eprop=0.3)
head(DMRs)
}

