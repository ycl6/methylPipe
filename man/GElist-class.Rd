\name{GElist-class}
\docType{class}
\alias{GElist-class}
\alias{GElist}
\alias{show,GElist-method}
\alias{length,GElist-method}

% Subset methods:
\alias{[[<-,GElist,ANY,ANY-method}
\alias{[[,GElist,ANY,ANY-method}
\alias{[,GElist,ANY,ANY-method}
\alias{$,GElist-method}
\alias{$,GElist}


\title{Class "GElist"}
\description{ This class is used in the methylPipe library to collect a
  set of GEcollection objects}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("GElist", ...)}
  or using the function GElist(Objlist,names), see below.
  GElist are a collection of GEcollection objects (see
  \code{\link{GElist-class}}).
}
\section{Slots}{
  \describe{
    \item{\code{Objlist}:}{Object of class \code{"list"} : a list where
      each item is a \link{GEcollection} object}
    \item{\code{names}:}{Object of class \code{"character"} : vector of
      the names of the objects}
  }
}
\section{Methods}{
  \describe{
    \item{"[["}{\code{signature(x = "GElist")}: subsets the
      \link{GElist} returning a specific \link{GEcollection} object}
    \item{"[[<-"}{\code{signature(x = "GElist")}: replaces the specific
      \link{GEcollection} object in the \link{GElist}}
    \item{"["}{\code{signature(x = "GElist")}: subsets the \link{GElist}
      returning another \link{GElist}}
  }
}
\author{Mattia Pelizzola}
\seealso{
  \code{\link{GElist-class}}
}
\examples{
gecollect_file <- system.file('extdata', 'gec.H1.Rdata', package='methylPipe')
load(gecollect_file)
gec1 <- gec.H1[start(gec.H1) < 153924]
gec2 <- gec.H1[start(gec.H1) > 153924]
gel.set <- GElist(g1=gec1, g2=gec2)
}
\keyword{classes}
