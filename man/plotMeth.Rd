\name{plotMeth}
\alias{plotMeth}
\title{Plot DNA methylation together with other omics, or annotation data for a genomic region}
\description{
  Plot DNA methylation data (either high- or low-resolution) together with other omics data (ChIP-seq, RNA-seq), or annotation tracks for one genomic region (genome-browser like view based on gviz).
}
\usage{
plotMeth(grl=NULL, colors=NULL, datatype=NULL, yLim, brmeth=NULL, mcContext="CG", annodata=NULL, 
transcriptDB, chr, start, end, org)
}
\arguments{
  \item{grl}{An object of class \link{GElist} or a potentially mixed \link{list} of \link{GRanges} or \link{GEcollection} objects}
  \item{colors}{character of length equal to grl; name of colors to display data tracks from the grl object}
  \item{datatype}{character of length equal to grl; one of C, mC , rC, density or cols}
  \item{yLim}{numeric vector with the same length of grl setting
    maximum values}  
  \item{brmeth}{A list of object of class \link{BSdata}}
  \item{mcContext}{character; one of all, CG, CHG or CHH} 
  \item{annodata}{An object of class \link{GRangesList}}
  \item{transcriptDB}{An object of class TranscriptDb}
  \item{chr}{character; chromosome name}
  \item{start}{numeric; chromosome start}
  \item{end}{numeric; chromosome end}
  \item{org}{BSgenome; an object of class BSgenome}
}
\details{
This function can be used to display for one genomic region (genome-browser like) DNA methylation data together with other omics data or static annotation info. The genomic region is indicated by chr, start and end. Specifically, grl is used to display binned high- or low-resolution data, while brmeth is used to point to (unbinned) base-resolution data. Each component of grl can either be a GEcollection or a GRanges.

In case of GEcollection, binC, binmC or binrC components will be extracted as indicated in datatype (setting C, mC or rC, respectively), and if more than 1 bin is present the average value will be considered for each range. Datatype can be set to density to extract the binscore component of the GEcollection, which can be used to store low-resolution or other omics data attacched to a base-resolution dataset.

In case of a GRanges (suitable for low-resolution or other omics data independently from base-resolution data), only the 1st column of the mcols will be considered. Eventually, for both GEcolleciton and GRanges tracks, a bar with the specific value will be displayed for the ranges occurring in the considered region (if any).

Regarding unbinned base-resolution data, mcContext defines the sequence context to be considered for the methyl-cytosines for each component of the brmeth object; a bar with height equal to the methylation level of each cytosine will be displayed for each sample (track).

Annodata is an optional GRangesList that can be used to display co-occurring annotation data, such as CpG islands (presence or absence of the regions only).
transcriptDB and BSgenome are used to overlay the structure of annotated genes and chromosome ideogram, respectively.
}
\author{Kamal Kishore}
\examples{
require(TxDb.Hsapiens.UCSC.hg18.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene
require(BSgenome.Hsapiens.UCSC.hg18)
gecH1_file <- system.file('extdata', 'gec.H1.Rdata', package='methylPipe')
load(gecH1_file)
gecIMR_file <- system.file('extdata', 'gec.IMR90.Rdata', package='methylPipe')
load(gecIMR_file)
gel <- GElist(gecH1=gec.H1, gecIMR90=gec.IMR90)
uncov_GR <- GRanges(Rle('chr20'), IRanges(c(14350,69251,84185), c(18349,73250,88184)))
H1data <- system.file('extdata', 'H1_chr20_CG_10k_tabix_out.txt.gz', package='methylPipe')
H1.db <- BSdata(file=H1data, uncov=uncov_GR, org=Hsapiens)
IMR90data <- system.file('extdata', 'IMR90_chr20_CG_10k_tabix_out.txt.gz', package='methylPipe')
IMR90.db <- BSdata(file=IMR90data, uncov=uncov_GR, org=Hsapiens)
H1.IMR90.set <- list(H1=H1.db, IMR90=IMR90.db)
plotMeth(gel, colors=c("red","blue"), datatype=c("mC","mC"), yLim=c(.025, .025), brmeth=H1.IMR90.set, mcContext="CG", transcriptDB=txdb, chr="chr20", start=14350 , end=474481, org=Hsapiens)
}
