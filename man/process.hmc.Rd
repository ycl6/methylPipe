\name{process.hmc}
\alias{process.hmc}
\title{Processing hmC information from the MLML output}
\description{
  Processing MLML software output to generate files with hmC and CpG methylation information.
}
\usage{
process.hmc(file, output_folder, Coverage)
}
\arguments{
  \item{file}{character; the path to the file}
  \item{output_folder}{character; the path to the output files}
  \item{Coverage}{GRanges; the object containing coverage for each cytosine}
}
\details{This function allows processing of the output files from MLML software (Qu et al, Bioinformatics 2013).
  MLML read counts from BS-seq, oxBS-seq and TAB-seq to provide simultaneous estimates of 5hmC and 5mC levels. The
  input for this function is output file from this software alongwith a \link{GRanges} object consisting  of coverage
  of each cytosine. The \link{GRanges} object should contain "coverage" column. This object can be generated using the
  \link{coverage} method of R package \link{GenomicRanges}.
}
\value{
  The function will return two files one each for "CpG" and "hmC" for the given sample which can directly be used for \link{BSdata} object creation.
}
\author{Kamal Kishore}
\seealso{
  \code{\link{BSdata-class}}
}
\examples{
#process.hmc(file,"/path-to-output/", Coverage)
}
