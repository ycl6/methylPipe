\name{findPMDs}
\alias{findPMDs}
\alias{findPMDs,methylPipe,BSdata}
\alias{findPMDs,methylPipe,BSdata-method}
\alias{findPMDs-methods}
\alias{findPMDs,BSdata-method}
\title{Identifying Partially Methylated Domains (PMDs)}
\description{
  This function is a wrapper function to identify partially methylated domains (PMDs) in Bis-seq data.
}
\usage{
\S4method{findPMDs}{methylPipe,BSdata}(Object, Nproc=1, Chrs=NULL)
}
\arguments{
  \item{Object}{An object of class \link{BSdataSet}}
  \item{Nproc}{numeric; the number of processors to use, one chromosome
    is ran for each processor}
  \item{Chrs}{character; Chromosome on which PMDs are identified}
}
\details{This functions is a wrapper function of segmentPMDs method of package MethylSeekR. 
  This function trains a Hidden Markov Model (HMM) to detect partially methylated domains (PMDs) in Bis-seq data.
}
\value{
      A GRangesList object containing segments that partition the genome
     into PMDs and regions outside of PMDs. The object contains two
     metadata columns indicating the type of region (PMD/notPMD) and
     the number of covered (by at least 5 reads) CpGs (nCG) in the
     region. 
     }
\author{Kamal Kishore}
\seealso{
  \code{\link{findDMR}}
}
\examples{
require(BSgenome.Hsapiens.UCSC.hg18)
uncov_GR <- GRanges(Rle('chr20'), IRanges(c(14350,69251,84185), c(18349,73250,88184)))
H1data <- system.file('extdata', 'H1_chr20_CG_10k_tabix_out.txt.gz', package='methylPipe')
H1.db <- BSdata(file=H1data, uncov=uncov_GR, org=Hsapiens)
PMDs <- findPMDs(H1.db, Nproc=1, Chrs="chr20")
}

