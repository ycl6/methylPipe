\name{consolidateDMRs}
\alias{consolidateDMRs}
\title{Consolidating Differentially Methylated Regions (DMRs)}
\description{
  Joins differentially methylated regions according to their proximity
  to each other, statistical significance and methylation difference
}
\usage{
consolidateDMRs(DmrGR, pvThr=0.05, MethDiff_Thr=NULL, log2Er_Thr =NULL, GAP=0, type=NULL, correct=FALSE)
}
\arguments{
  \item{DmrGR}{the GRanges object resulting from the
    \code{\link{findDMR}} function}
  \item{pvThr}{numeric; the minimum pvalue for a DMR to be selected}
  \item{MethDiff_Thr}{numeric; the absolute value of minimum methylation
    difference percentage (for both hyper- and hypo-methyaltion) cutoff
    for the selection of a DMR}
  \item{log2Er_Thr}{numeric; the absolute value of minimum log2Enrichment
   (for both hyper- and hypo-methyaltion) cutoff for the selection of a DMR}
  \item{GAP}{numeric; the minimum distance between two adjacent DMRs;
    DMRs closer than that will be joined, resulting DMRs will be
    updated mean methylation difference and Pvalues combined using the
    Fisher's Method}
  \item{type}{character; one of the "hyper" or "hypo", specifies the type
    of differentially menthylated regions selected}
  \item{correct}{logical; whether to correct the pvalues using the
    Benjamini-Hochberg muliple testing correction method}
}
\details{After the DMRs are identified by \link{findDMR} method, a
  consolidation can be applied on them. This functions allows to select
  hyper/hypo differentially methylated regions based on P-value and
  absolute methylation change thresholds. Moreover, DMRs closer than a given
  distance can be joined. The final \link{GRanges} object
  with the set of final DMRs will be provided with updated mean
  methylation difference and Pvalues combined using the Fisher's Method.
}
\value{
  Either NULL or a \link{GRanges} object containing the differential
  methylated regions satisfying the criteria.
}
\author{Kamal Kishore}
\seealso{
  \code{\link{findDMR}}
}
\examples{
DMRs_file <- system.file('extdata', 'DMRs.Rdata', package='methylPipe')
load(DMRs_file)
hyper.DMRs.conso <- consolidateDMRs(DmrGR=DMRs, pvThr=0.05, GAP=100, type="hyper", correct=TRUE)
hyper.DMRs.conso
}

