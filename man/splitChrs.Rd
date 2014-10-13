\name{splitChrs}
\alias{splitChrs}
\title{Partitioning genome in chunks, for parallel computation}
\description{
  Helper function to partition genome chromosome-wise for parallel computation. This
  function is not intended for the user to call directly.
}
\usage{
splitChrs(chrs, org)
}
\arguments{
  \item{chrs}{character; an aray of chromome names in the form chr1, .., chrX}
  \item{org}{an object of class BSgenome}
}
\value{
  A data frame with chromosome name, start and end position of each chunk.
}


