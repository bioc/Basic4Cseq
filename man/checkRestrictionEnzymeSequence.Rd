\name{checkRestrictionEnzymeSequence}

\alias{checkRestrictionEnzymeSequence}
\alias{checkRestrictionEnzymeSequence,character,character-method}

\title{Remove invalid 4C-seq reads from a SAM file}

\description{
\pkg{Basic4Cseq} offers filter functions for invalid 4C-seq reads. This function removes 4C-seq reads from a provided Sequence Alignment/Map (SAM) file that show mismatches in the restriction enzyme sequence. 
}

\usage{checkRestrictionEnzymeSequence(firstCutter, inputFileName, outputFileName = "output.sam", keepOnlyUniqueReads = TRUE, writeStatistics = TRUE)}

\arguments{
  \item{firstCutter}{First restriction enzyme sequence of the 4C-seq experiment}
  \item{inputFileName}{Name of the input SAM file that contains aligned reads for the 4C-seq experiment}
  \item{outputFileName}{Name of the output SAM file that is created to store the filtered 4C-seq reads}
  \item{keepOnlyUniqueReads}{If TRUE, delete non-unique reads. Information in the SAM flag field is used to determine whether a read is unique or not.}
  \item{writeStatistics}{If TRUE, write statistics (e.g. the number of unique reads) to a text file}
}

\details{
Valid 4C-seq reads start at a primary restriction site and continue with its downstream sequence, so any mismatch in the restriction enzyme sequence of a read is an indicator for a mismatch. The mapping information of the restriction enzyme sequence bases of a read (if present) can be used for filtering purposes. \code{checkRestrictionEnzymeSequence} tests the first bases of a read (depending on the length of the first restriction enzyme either 4 or 6 bp long) for mismatches. Reads with mismatches in the restriction enzyme sequence are deleted, the filtered data is then written to a new SAM file. The function does not yet differentiate between blind and nonblind fragments, but removes potential misalignments that may overlap with valid fragment ends and distort the true 4C-seq signal.
}

\value{A SAM file containing the filtered valid 4C-seq reads}

\note{The use of the function is only possible if the restriction enzyme sequence is not trimmed or otherwise absent.}

\examples{
  if(interactive()) {
    file <- system.file("extdata", "fetalLiverCutter.sam", package="Basic4Cseq")
    checkRestrictionEnzymeSequence("aagctt", file)
  }
}

\author{Carolin Walter}

\keyword{checkRestrictionEnzymeSequence}
