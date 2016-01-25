\name{BSprepare}
\alias{BSprepare}
\title{Preparing tabular data to be used to feed a BSdata object}
\description{
  Appending p-values and TABIX compressing tabular data containing
  DNA-methylation data so that they can be used to create a BSdata
  object.
}
\usage{
BSprepare(files_location, output_folder, tabixPath, bc=1.5/100)
}
\arguments{
  \item{files_location}{character; the path to the files}
  \item{output_folder}{character; the path to the output files}
  \item{tabixPath}{character; the path to the Tabix application folder}
  \item{bc}{numeric; combined bisulfite conversion and sequencing error rate}
}
\details{This function can be used to convert tabular files containing
  DNA-methylation base-resolution data so that they can be used to
  create a BSdata object. Genomic coordinates in the 1-base system are
  required, i.e. the first base of each chromosome should be at position 1.
  Files have to be tab-delimited and they have to contain the
  following columns:
  chromosome assignment (in the form chr1, .., chr22, chrX, chrY, chrM, chrC),
  genomic position (positive integer),
  strand (either - or +),
  methylation sequence context (either CG, CHG or CHH),
  number of sequencing reads with C calls (>0) at that genomic position,
  number of sequencing reads with T calls at that genomic position.
  Each file has to be sorted by genomic coordinate.}
\value{
  Binomial p-values are added and a compressed file is created together with the Tabix index.
  P-values indicate for each Cytosine the probability of observing by chance the occurrence of unmethylated reads.
  The lower the p-value the higher the confidence in calling that Cytosine methylated.
}
\author{Mattia Pelizzola, Kamal Kishore}
\seealso{
  \code{\link{BSdata-class}}
}
\examples{
#BSprepare("/path-to-input/","/path-to-output/", tabix="/path-to-tabix/tabix-0.2.6")
}
