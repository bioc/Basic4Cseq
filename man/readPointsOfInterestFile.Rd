\name{readPointsOfInterestFile}

\alias{readPointsOfInterestFile}
\alias{readPointsOfInterestFile,character-method}

\title{Read a file with coordinates of marker points}

\description{
This small helper function reads a tab-separated file with points of interest information stored in a BED-like format. The file has to provide the columns "chromosome", "start", "end", "name" and "colour" of the regions. The data can then be used for marking the points in near-cis visualization plots, as provided by \code{visualizeViewpoint} and \code{drawHeatmap}.
}

\usage{readPointsOfInterestFile(poiFile)}

\arguments{
  \item{poiFile}{Name of the input file (tab-separated)}
}

\value{Data frame with information on points of interest for the near-cis visualizations}

\examples{
  file <- system.file("extdata", "fetalLiverVP.bed", package="Basic4Cseq")
  pointsOfInterests = readPointsOfInterestFile(file)
  pointsOfInterests
}

\author{Carolin Walter}

\keyword{readPointsOfInterestFile}

