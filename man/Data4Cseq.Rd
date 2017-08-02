\name{Data4Cseq}

\alias{Data4Cseq}
\alias{Data4Cseq,character,numeric,numeric,data.frame,GAlignments-method}
\alias{Data4Cseq,character,numeric,numeric,missing,missing-method}


\title{Creating a Data4Cseq object}

\description{
    This function creates a Data4Cseq object. Data on the 4C-seq experiment, e.g. the chromosome of the viewpoint, is stored and checked for consistency.
}

\usage{Data4Cseq(viewpointChromosome, viewpointInterval, readLength, pointsOfInterest, rawReads)}

\arguments{
    \item{viewpointChromosome}{The experiment's viewpoint chromosome}
    \item{viewpointInterval}{The interval of the experiment's viewpoint, consisting of a start and end coordinate}
    \item{readLength}{The experiment's read length (in base pairs)}
    \item{pointsOfInterest}{Points of interest to be marked in a near-cis visualization}
    \item{rawReads}{Reads of the 4C-seq experiment, aligned and stored as an \code{GAlignments} object}
}

\details{
    A \code{Data4Cseq} object contains basic information for the corresponding 4C-seq experiment, including the viewpoint chromosome, the viewpoint region and reads from the experiment. See \code{\link{Data4Cseq-class}} for more details. The constructor collects the basic data; fragment data or normalized read counts are added later.   

Fragments at the experiment's viewpoint are usually vastly overrepresented due to self-ligation; \code{chooseNearCisFragments} offers the option to discard all fragments in the specified viewpoint region. The specified viewpoint interval of a \code{Data4Cseq} object is supposed to correspond to the positions of the biological primers on the genome, but can also be increased in size if more fragments around the viewpoint should be removed.
}

\value{
    An instance of the Data4Cseq class.
}

\seealso{
    \code{\link{Data4Cseq-class}}
}

\examples{
    # create a Data4Cseq object with a minimum of data
    liverData = Data4Cseq(viewpointChromosome = "10", viewpointInterval = c(20879870, 20882209), readLength = 54)
    liverData

    # create a Data4Cseq object, including possible points of interest and raw reads 
    bamFile <- system.file("extdata", "fetalLiverShort.bam", package="Basic4Cseq")
    liverReads <- readGAlignments(bamFile)
    pointsOfInterestFile <- system.file("extdata", "fetalLiverVP.bed", package="Basic4Cseq")
    liverData = Data4Cseq(viewpointChromosome = "10", viewpointInterval = c(20879870, 20882209), readLength = 54, pointsOfInterest = readPointsOfInterestFile(pointsOfInterestFile), rawReads = liverReads)
    liverData
}

\author{Carolin Walter}

\keyword{Data4Cseq}
