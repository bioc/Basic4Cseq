\name{drawDigestionFragmentHistogram}

\alias{drawDigestionFragmentHistogram}
\alias{drawDigestionFragmentHistogram,data.frame-method}

\title{Visualize digestion fragments with a histogram}

\description{
This function is a small convenience function to plot the results of \code{simulateDigestion} as a histogram. Minimum and maximum fragment lengths can be specified to visualize a specified interval of the fragment data.
}

\usage{
  drawDigestionFragmentHistogram(fragments, minLength = 0, maxLength = 10000)
}

\arguments{
  \item{fragments}{Fragment data to visualize (data frame with lengths and corresponding frequencies)}
  \item{minLength}{Minimum fragment length to visualize}
  \item{maxLength}{Maximum fragment length to visualize}
}

\value{Histogram plot of the fragment data}

\examples{
  shortTestGenome = "ATCCATGTAGGCTAAGTACACATGTTAAGGTACAGTACAATTGCACGATCAT"
  fragments = simulateDigestion("catg", "gtac", shortTestGenome)
  drawDigestionFragmentHistogram(fragments)
}

\author{Carolin Walter}

\keyword{drawDigestionFragmentHistogram}

