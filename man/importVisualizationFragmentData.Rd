\name{importVisualizationFragmentData}

\alias{importVisualizationFragmentData}
\alias{importVisualizationFragmentData,character-method}

\title{Import visualization data from a file}

\description{
This function is a simple helper function that can import near-cis data which was previously exported and stored as tab-separated file.
}

\usage{importVisualizationFragmentData(fileName)}

\arguments{
    \item{fileName}{Name for the tab-separated file with near-cis fragment data}
}

\value{Data frame containing the near-cis fragment data}

\examples{
    file <- system.file("extdata", "fetalLiver_finalFragments.csv", package="Basic4Cseq")
    importVisualizationFragmentData(file)
    head(file)
}

\author{Carolin Walter}

\keyword{importVisualizationFragmentData}

