\name{chiCombP}
\alias{chiCombP}
\title{Fisher's method implementation}
\description{
  Fisher's method implementation, used to combine the results from
  several independent tests bearing upon the same overall hypothesis.
}
\usage{
chiCombP(Pvalues)
}
\arguments{
  \item{Pvalues}{an array of pvalues}
}
\details{Pvalues will not be corrected for multiple testing. The sum of
  the log of the pvalues is determined (S). -2*S has a chi-square
  distribution with 2k degrees of freedom, where k is the number of
  tests being combined. A chi-square test is then performed.}
\value{
  The chi-square final pvalue
}
\author{Mattia Pelizzola}
\examples{
chiCombP(c(1e-3,1e-5,1e-2))
}

