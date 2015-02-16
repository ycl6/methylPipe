\name{methylPipe-package}
\alias{methylPipe-package}
\alias{methylPipe}
\docType{package}
\title{Analysis of base-pair resolution DNA methylation data.}
\description{Analysis of base-pair resolution DNA methylation data.}
\details{
\tabular{ll}{
Package: \tab methylPipe\cr
Type: \tab Package\cr
Version: \tab 0.99\cr
Date: \tab 2014-06-17\cr
License: \tab GPL\cr
Depends: \tab methods\cr
}
The package offers the following functionalities:

\itemize{
\item BSdata-class :  This class is used in to point to a TABIX compressed file containing base-resolution DNA-methylation data and reference genome sequence  
\item mCsmoothing : Smoothing and plotting methylation data, even chromosome wide	    
\item findPMDs : Find partially methylated regions for a given sample    
\item mapBSdata2GRanges : Retrieve mC calls for a GRanges from a BSdata object for a sample 
\item BSdataSet-class : This class is a set of BSdata objects
\item findDMR : Identifying differentially methylated regions for pairwise or multiple samples comparision
\item methstats : Descriptive methylation statistics of samples within BSdataSet object
\item consolidateDMRs : Joins differentially methylated regions according to their proximity to each other, statistical significance, methylation difference and DMR type
\item GEcollection-class : This class is used to define and manipulate a set of genomic regions and the associated DNA methylation patterns
\item getCpos : Get genomic Cxx positons for a series of genomic regions	       
\item getCposDensity : Determines the density of genomic Cxx positions for a series of genomic regions
\item profileDNAmetBin : Profile DNA methylation data for a set of genomic regions
\item plotMeth : Plot the methylation/omics data of a GEcollection object	    
\item BSprepare : Preparing tabular data to be used to feed a BSdata object	    
\item meth.call : Reads the methylation information from the sorted SAM files generated from BISMARK aligner	     
\item GElist-class : This class is a set of GEcollection objects	
} 
}
\author{
Computational Epigenomics Unit at the Center for Genomic Sciences of IIT@SEMM, Milan, Italy

\email{kamal.fartiyal84@gmail.com}

\url{http://genomics.iit.it/groups/computational-epigenomics.html}
}
\keyword{ package }

