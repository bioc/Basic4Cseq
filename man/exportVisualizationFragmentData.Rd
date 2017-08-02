\name{exportVisualizationFragmentData}

\alias{exportVisualizationFragmentData}
\alias{exportVisualizationFragmentData,Data4Cseq,character-method}

\title{Export near-cis fragment data of a \code{Data4Cseq} object}

\description{
This function is a simple helper function that writes the near-cis data of a \code{Data4Cseq} object as tab-separated file to hard disk.
}

\usage{exportVisualizationFragmentData(expData, fileName, fullData = FALSE)}

\arguments{
    \item{expData}{Experiment data of class \code{Data4Cseq} information on the 4C-seq experiment, including visualization data}
    \item{fileName}{Name for the tab-separated file}
    \item{fullData}{If TRUE, the function exports the full fragment data (including fragment end length etc). If FALSE, only the minimum fragment information is exported, i.e. chromosome, start, end and (normalized) read count.}
}

\value{A tab-separated file containing near-cis framgent data of a \code{Data4Cseq} object}

\examples{
    if(interactive()) {
        data(liverData)
        exportVisualizationFragmentData(liverData, "fetalLiverData.csv")
    }
}

\author{Carolin Walter}

\keyword{exportVisualizationFragmentData}

