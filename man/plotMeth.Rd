\name{plotMeth}
\alias{plotMeth}
\title{Plot the methylation/omics data of a GEcollection object}
\description{
  Plot the the methylation/omics data for the genomic regions of a \link{GEcollection} object with locus specific visualization.
}
\usage{
plotMeth(grl, colors=NULL, datatype, yLim, brmeth=NULL, mcContext="CG", annodata=NULL, Datatrackname, 
                     transcriptDB, chr, start=NULL, end=NULL, org)
}
\arguments{
  \item{grl}{An object of class \link{GElist} or \link{list} of objects of class \link{GRanges} or class \link{GEcollection}}
  \item{colors}{character; name of colors to display data tracks from the grl object}
  \item{datatype}{character; one of C, mC , rC, density, cols or gr,
    indicating that you want to plot the average C* content, mC*
    absolute methylation level, mC*/C* relative methylation level, omics
    data from binscore slot of GEcollection/GElist object, data from metadata 
    column (mcols) of GEcollection/GElist object or data from GRanges columns 
    respectively}
  \item{yLim}{numeric; numeric vector equal to length of grl to define
    maximum range of y coordinates}  
  \item{brmeth}{An object of class \link{BSdataSet}}
  \item{mcContext}{character; one of: all, CG, CHG or CHH} 
  \item{annodata}{An object of class \link{GRanges}}
  \item{Datatrackname}{character; Track's name used when plotting data
    tracks}
  \item{transcriptDB}{An object of class TranscriptDb}
  \item{chr}{character; chromosome name}
  \item{start}{numeric; chromosome start}
  \item{end}{numeric; chromosome end}
  \item{org}{BSgenome; an object of class BSgenome}
}
\details{This function requires that the binC, binmC, binrC or binscore slot of
  the \link{GEcollection} were previously filled, using the \code{\link{profileDNAmetBin}} function.}
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
H1.IMR90.set <- BSdataSet(org=Hsapiens, group=c("E","C"), H1=H1.db, IMR90=IMR90.db)
plotMeth(gel, datatype=c("mC","mC"), yLim=c(.025, .025), brmeth=H1.IMR90.set, mcContext="CG", 
Datatrackname=c("Abs meth","Abs meth","H1","IMR90"), transcriptDB=txdb, chr="chr20", org=Hsapiens)
}
