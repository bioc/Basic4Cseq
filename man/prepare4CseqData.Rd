\name{prepare4CseqData}

\alias{prepare4CseqData}
\alias{prepare4CseqData,character,character,character,character-method}

\title{Alignment and filtering of raw 4C-seq data}

\description{
This function is an optional wrapper for the alignment and preliminary filtering of 4C-seq data. \code{prepare4CseqData} reads a provided 4C-seq fastq file from hard disk. Alignment of the reads is done with BWA, the function \code{checkRestrictionEnzymeSequence} is used for optional filtering. Samtools and bedtools provide the necessary functionality for intersecting the filtered reads with a given 4C-seq fragment library for visualization purposes (e.g. with the Integrative Genomics Viewer, IGV).
}

\usage{prepare4CseqData(fastqFileName, firstCutter, fragmentLibrary, referenceGenome, pathToBWA = "", pathToSam = "", pathToBED = "", controlCutterSequence = FALSE, bwaThreads = 1, minFragEndLength = 0)
}

\arguments{
  \item{fastqFileName}{The name of the fastq file that contains the 4C-seq reads}
  \item{firstCutter}{First cutting enzyme sequence for the 4C-seq experiment, e.g. "AAGCTT"}
  \item{fragmentLibrary}{Name of the fragment library to use for the current 4C-seq experiment; has to correspond to the chosen cutters and chosen genome}
  \item{referenceGenome}{Name (plus path) of the reference genome to use}
  \item{pathToBWA}{Path to BWA}
  \item{pathToSam}{Path to samtools}
  \item{pathToBED}{Path to bedtools}
  \item{controlCutterSequence}{If TRUE, the function \code{checkRestrictionEnzymeSequence} is used to filter non-valid 4C-seq reads}
  \item{bwaThreads}{Number of BWA threads}
  \item{minFragEndLength}{Minimum fragment end length to use for BED export}
}

\seealso{\code{\link{checkRestrictionEnzymeSequence}} 
}

\references{
Li, H. and Durbin, R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform, Bioinformatics, 25, 1754-60.

Helga Thorvaldsdottir, James T. Robinson, Jill P. Mesirov. Integrative Genomics Viewer (IGV): high-performance genomics data visualization and exploration.  Briefings in Bioinformatics 2012.
}

\examples{
  if(interactive()) {
    # BWA, samtools and bedtools must be installed
    # It is assumed that the example data files (from the package) are in the active directory
    prepare4CseqData("veryShortExample.fastq", "CATG", "veryShortLib.csv", referenceGenome = "veryShortReference.fasta")
  }
}

\author{Carolin Walter}

\keyword{prepare4CseqData}
