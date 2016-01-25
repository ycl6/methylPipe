\name{pool.reads}
\alias{pool.reads}
\title{Function to pool reads of replicates}
\description{
  Combine reads of replicates within a group.
}
\usage{
pool.reads(files_location)
}
\arguments{
  \item{files_location}{character; the path to the folder location consisting of tab separated text files}
}
\details{The function reads tab separated text files of methylation calls generated from \link{meth.call} or user supplied
  according to the format specified for \link{BSprepare} method. It pools all the reads of the replicates within a single
  group for each cytosine position and creates a file consisting of the cytosines with pooled reads information.
}
\value{
  A text file of methylation calls are generated in the "files_location" folder.
}
\author{Kamal Kishore}
\seealso{
  \code{\link{BSprepare}}
}
\examples{
#pool.reads(files_location)
}
