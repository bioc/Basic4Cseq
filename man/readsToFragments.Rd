\name{readsToFragments}

\alias{readsToFragments}
\alias{readsToFragments,Data4Cseq,character-method}

\title{Determine fragment coverage of a 4C-seq fragment library}

\description{
    This function maps aligned reads to fragment ends of the virtual fragment library to calculate the coverage of the fragments. The number of reads at the start and end of a fragment is provided, as well as the average of both fragment ends.
}

\usage{readsToFragments(expData, fragmentLib)}

\arguments{
    \item{expData}{Experiment data of class \code{Data4Cseq} with information on the 4C-seq experiment, including raw 4C-seq read data}
    \item{fragmentLib}{Fragment library for the given genome and cutting enzyme combination}
}

\value{
    Data frame containing fragment-based data, i.e. a fragment's position and read coverage 
}

\examples{
    data(liverData)
    file <- system.file("extdata", "vfl_aagctt_catg_mm9_54_vp.csv", package="Basic4Cseq")
    rawFragments(liverData) = readsToFragments(liverData, file)
    head(rawFragments(liverData))
}

\author{Carolin Walter}

\keyword{readsToFragments}
