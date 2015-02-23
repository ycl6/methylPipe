\name{mapBSdata2GRanges}
\alias{mapBSdata2GRanges}
\alias{mapBSdata2GRangesBin}
\title{Retrieve mC calls for a GRanges set of genomic regions given a BSdata object for a sample}
\description{
  mapBSdata2GRanges retrieves mC calls for a \link{GRanges} given a \link{BSdata} object for a sample.
  mapBSdata2GRangesBin does the same for each bin of each genomic region.
}
\usage{
mapBSdata2GRanges(GenoRanges, Sample, context='all', mC=1, depth=0, pValue=1)
mapBSdata2GRangesBin(GenoRanges, Sample, context='all', mC=1, depth=0, pValue=1, nbins)
}
\arguments{
  \item{GenoRanges}{An object of class \link{GRanges}}
  \item{Sample}{An object of class \link{BSdata}}
  \item{context}{character; one of: all, CG, CHG or CHH}
  \item{mC}{numeric; the minumum number of reads with C (DNA-methylation
    events) at a given cytosine genomic position}
  \item{depth}{numeric; the minimum number of total reads at a given
    cytosine genomic position}
  \item{pValue}{numeric; the minimum binomial pValue for an mC call at a
    given cytosine genomic position}
  \item{nbins}{numeric; the number of bins in which Genomic regions are
    divided}
}
\details{DNA-methylation data contained for a sample within a
  \link{BSdata} object is extracted for the set of genomic regions of a
  \link{GRanges} (and in particular for each bin using the
  mapBSdata2GRangesBin method).
  It is also possible to restrict the mC sequence context, to specify the
  minimal number of reads with C events at a given cytosine genomic
  position, to specify the minimum depth of sequencing and binomial
  pValue for the mC calls. A region with no mC will be defined
  unmethylated (0 is returned for that region). However, if it is
  overlapping with at least one uncovered region then it is defined non
  evaluable (NA is retuned).}
\value{
  A list is returned. The length of the list is equal to the length of
  the \link{GRanges}. For each list item either NA, 0 or a data frame
  are returned.
  0 means that the region contains unmethylated DNA methylation sites,
  whereas NA indicates that the region or some part of region was not
  covered by the sequencing.
  If a data frame is returned, it has the following columns:
  chromosome assignment (in the form chr1, .., chr22, chrX, chrY, chrM, chrC),
  genomic position (positive integer), strand (either - or +),
  methylation sequence context (either CG, CHG or CHH),
  number (>0) of sequencing reads with C calls at that genomic position,
  number of sequencing reads with T calls at that genomic position,
  binomial pvalue (-10*log10(pvalue)) for calling a mC at that position.
}
\author{Mattia Pelizzola, Kamal Kishore}
\examples{
require(BSgenome.Hsapiens.UCSC.hg18)
uncov_GR <- GRanges(Rle('chr20'), IRanges(c(14350,69251,84185), c(18349,73250,88184)))
H1data <- system.file('extdata', 'H1_chr20_CG_10k_tabix_out.txt.gz', package='methylPipe')
H1.db <- BSdata(file=H1data, uncov=uncov_GR, org=Hsapiens)
gr_file <- system.file('extdata', 'GR_chr20.Rdata', package='methylPipe')
load(gr_file)
res <- mapBSdata2GRanges(GenoRanges=GR_chr20, Sample=H1.db, context='CG', mC=1, depth=0, pValue=1)
resbin <- mapBSdata2GRangesBin(GenoRanges=GR_chr20, Sample=H1.db, context='CG', mC=1, depth=0, pValue=1, nbins=2)
}

