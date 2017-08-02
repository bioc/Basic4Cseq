\name{normalizeFragmentData}

\alias{normalizeFragmentData}
\alias{normalizeFragmentData,Data4Cseq-method}

\title{Normalize near-cis fragment data read count}

\description{
This function provides a simple RPM (reads per million) normalization for near-cis fragment data read counts of a \code{Data4Cseq} object. A form of normalization is especially important for the comparison of samples with a different read count.
}

\usage{normalizeFragmentData(expData)}

\arguments{
    \item{expData}{Experiment data of class \code{Data4Cseq} with information on the 4C-seq experiment}
}

\value{Data frame with RPM-normalized data}

\examples{
    data(liverDataRaw)
    normalizedFragments<-normalizeFragmentData(liverDataRaw)
    head(normalizedFragments)
}

\author{Carolin Walter}

\keyword{normalizeFragmentData}

