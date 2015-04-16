\name{BSdataSet-class}
\docType{class}
\alias{BSdataSet-class}
\alias{BSdataSet}
\alias{length,BSdataSet-method}
\alias{show,BSdataSet-method}

% Subset methods:
\alias{[[<-,BSdataSet,ANY,ANY-method}
\alias{[[,BSdataSet,ANY,ANY-method}
\alias{[,BSdataSet,ANY,ANY-method}
\alias{$,BSdataSet-method}
\alias{$,BSdataSet}


\title{Class "BSdataSet"}
\description{ This class is used in the methylPipe library to store a
  set of BSdata objects}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("BSdataSet",
    ...)} or using the function BSdataSet(org,Objlist,names), see below.}
\section{Slots}{
  \describe{
    \item{\code{Objlist}:}{Object of class \code{"list"} : a list with
      more than one item, where each item is a \link{BSdata} object}
    \item{\code{names}:}{Object of class \code{"character"} : vector of
      the names of the objects, i.e. the sample names}
    \item{\code{group}:}{Object of class \code{"character"} : indicating
conditions and replicates; replicated samples within the same condition
have to be assigned the same group name; if object has only two groups, 
"C"(control) and "E" (Experiment) should be specified as groups name}
    \item{\code{org}:}{refrence genome of class BSgenome}
  }
}
\section{Methods}{
  \describe{
    \item{"[["}{\code{signature(x = "BSdataSet")}: subsets the
      \link{BSdataSet} returning a specific \link{BSdata} object}
    \item{"[[<-"}{\code{signature(x = "BSdataSet")}: replaces the
      specific \link{BSdata} object in the \link{BSdataSet}}
    \item{"["}{\code{signature(x = "BSdataSet")}: subsets the
      \link{BSdataSet} returning another \link{BSdataSet}}
  }
}
\author{Mattia Pelizzola, Kamal Kishore}
\seealso{
  \code{\link{BSdata-class},\link{findDMR}}
}
\examples{
require(BSgenome.Hsapiens.UCSC.hg18)
uncov_GR <- GRanges(Rle('chr20'), IRanges(c(14350,69251,84185), c(18349,73250,88184)))
H1data <- system.file('extdata', 'H1_chr20_CG_10k_tabix_out.txt.gz', package='methylPipe')
H1.db <- BSdata(file=H1data, uncov=uncov_GR, org=Hsapiens)
IMR90data <- system.file('extdata', 'IMR90_chr20_CG_10k_tabix_out.txt.gz', package='methylPipe')
IMR90.db <- BSdata(file=IMR90data, uncov=uncov_GR, org=Hsapiens)
H1.IMR90.set <- BSdataSet(org=Hsapiens, group=c("C","E"), IMR90=IMR90.db, H1=H1.db)
}
\keyword{classes}

