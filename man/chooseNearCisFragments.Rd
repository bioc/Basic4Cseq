\name{chooseNearCisFragments}

\alias{chooseNearCisFragments}
\alias{chooseNearCisFragments,Data4Cseq,numeric-method}

\title{Choose fragments in a provided region around the viewpoint}

\description{
    This function extracts fragment data from a \code{Data4Cseq} object's \code{rawFragments} slot for visualization with the functions \code{visualizeViewpoint} and \code{drawHeatmap} . Relevant fragments are located within the chosen visualization range; the viewpoint itself can be excluded or included.
}

\usage{chooseNearCisFragments(expData, regionCoordinates, deleteViewpoint = TRUE)}

\arguments{
  \item{expData}{Experiment data of class \code{Data4Cseq} with information on the 4C-seq experiment, including fragment data for the viewpoint chromosome}
  \item{regionCoordinates}{Interval on the viewpoint chromosome for the intended visualization}
  \item{deleteViewpoint}{If TRUE, delete all fragments that intersect with the experiment's viewpoint interval}
}

\value{A data frame containing the chosen near-cis fragments}

\note{Viewpoint fragments are removed per default to prevent bias through overrepresented sequences caused by self-ligation. These fragments can be included, but should be interpreted with caution.}

\examples{
  # read example data
  data(liverData)
  fragments<-chooseNearCisFragments(liverData, regionCoordinates = c(20800000, 21000000))
  head(fragments)
}

\author{Carolin Walter}

\keyword{chooseNearCisFragments}
