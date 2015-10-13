\name{meth.call}
\alias{meth.call}
\title{Function to read methylation calls}
\description{
  Reads the methylation calls from sorted SAM files generated from
  Bismark aligner.
}
\usage{
meth.call(files_location, output_folder, no_overlap, read.context, Nproc)
}
\arguments{
  \item{files_location}{character; the path(s) to the folder location
    consisting of sorted SAM files}
  \item{output_folder}{character; the path(s) to the folder location
    where the output files are written}
  \item{no_overlap}{character; if set to TRUE and the SAM file has
    paired-end reads, then one read of the overlapping paired-end read
    pair will be ignored for methylation calling}
  \item{read.context}{character; One of the 'CpG' or 'All'. Determines
    what type of methylation context will be read. If given as 'all',
    cytosine methylation information in all sequence context will be read.}
  \item{Nproc}{numeric; the number of processors to use, one sample
    is processed by each processor.}
}
\details{The function reads methylation calls from the sorted SAM file
  so that they can be used to create a \link{BSdata} object. SAM files
  must be sorted based on chr and start of reads. The user can specify
  the sequence context in which the methylation information is read from
  these files either "CpG" or "All". If "All" is specified, cytosine
  methylation in all context (CG, CHG or CHH) will be read. The
  methylation calls is saved as a text file in the output folder.
  These text files are tab-delimited and contain the following columns:
  chromosome assignment (in the form chr1, chr2..),
  genomic position (positive integer),
  strand (either - or +),
  methylation sequence context (either CG, CHG or CHH),
  number (>0) of sequencing reads with C calls at that genomic position,
  number of sequencing reads with T calls at that genomic position.
  In addition a GRanges object consisting of genomic regions with
  sequencing coverage information is generated and saved in the output 
  folder for each sample. This information is used to distinguish unmethylated 
  cytosines from those that are not covered by sequencing. This GRanges object 
  is used further to provide uncovered regions information while creating BSdata
  object by \code{\link{BSdata}} method.
}
\value{
  A text file of methylation calls and a \link{GRanges} object
  consisting of uncovered genomic regions for each sample are generated
  in the "output_folder" folder. The files are prefixed with sample name.
}
\author{Kamal Kishore}
\seealso{
  \code{\link{BSprepare}}
}
\examples{
file_loc <- system.file('extdata', 'test_methcall', package='methylPipe')
meth.call(files_location=file_loc, output_folder=tempdir(), no_overlap=TRUE, read.context="CpG", Nproc=1)
}
