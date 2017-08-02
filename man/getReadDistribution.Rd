\name{getReadDistribution}

\alias{getReadDistribution}
\alias{getReadDistribution,Data4Cseq-method}

\title{Calculate the read distribution for a 4C-seq experiment}

\description{
This function provides some 4C-seq quality statistics based on the experiment's read distribution throughout the genome. \code{getReadDistribution} calculates the number of total reads, cis to overall ratio of reads, and the percentage of covered fragment ends within a certain distance around the experiment's viewpoint. Reference values for high-quality experiments, as provided by van de Werken et al, 2012, are more than one million reads total, a cis to overall ratio of more than 40\% and a large fraction of covered fragment ends in the viewpoint's vicinity. 
}

\usage{getReadDistribution(expData, distanceFromVP = 100000, useFragEnds = TRUE, outputName = "")
}

\arguments{
    \item{expData}{Experiment data of class \code{Data4Cseq} with information on the 4C-seq experiment}
    \item{distanceFromVP}{Distance from the viewpoint that is checked for covered fragments}
    \item{useFragEnds}{If TRUE, the function uses fragment end data; if FALSE, an average value for whole fragments is used.}
    \item{outputName}{An optional name for an output text file containing the statistics data}
}

\value{Text with statistics data on the 4C-seq experiment}

\references{van de Werken, H., de Vree, P., Splinter, E., et al. (2012): 4C technology: protocols and data analysis, Methods Enzymology, 513, 89-112}

\note{
Text export is supported; if no file name is provided, the results are printed on screen.
}

\examples{
    data(liverData)
    getReadDistribution(liverData)
}

\author{Carolin Walter}

\keyword{getReadDistribution}

