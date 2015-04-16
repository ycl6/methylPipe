\name{tabixdata2GR}
\alias{tabixdata2GR}
\title{Convert the list returned by the function scanTabix into a GRanges}
\description{
  Helper function to convert the list returned by the function scanTabix into a GRanges. This
  function is not intended for the user to call directly.
}
\usage{
tabixdata2GR(x)
}
\arguments{
  \item{x}{list; the list returned by the function scanTabix}
}
\value{
  An object of class data frame.
}
\author{Mattia Pelizzola, Kamal Kishore}
\seealso{
  \link{BSdata-class}
}
