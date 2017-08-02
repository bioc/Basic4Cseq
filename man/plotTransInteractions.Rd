\name{plotTransInteractions}

\alias{plotTransInteractions}
\alias{plotTransInteractions,character,character,numeric,character-method}
\alias{plotTransInteractions,data.frame,character,numeric,data.frame-method}

\title{Visualize trans interaction intervals}

\description{
This function visualizes trans interaction intervals of a 4C-seq experiment with the help of the \pkg{RCircos} package. Significant interactions can be obtained by use of Splinter et al's \code{significant_interactions} code or similar algorithms.
}

\usage{
plotTransInteractions(interactionFile, chromosomeViewpoint, coordViewpoint, ideogramData, PlotColor = "default", expandBands = FALSE, expansionValue = 0, plotFileName = "", picDim = c(8, 8))
}

\arguments{
    \item{interactionFile}{Interaction interval data; either a file name or a data frame}
    \item{chromosomeViewpoint}{Viewpoint chromosome of the 4C-seq experiment}
    \item{coordViewpoint}{Viewpoint coordinates of the 4C-seq experiment}
    \item{ideogramData}{Ideogram data to be visualized in the RCirco-plot; either a file name or a data frame}
    \item{PlotColor}{Plot colours for the visualized interactions}
    \item{expandBands}{If TRUE, add a specified value to the size of the interaction intervals to increase the visibility of very small interactions}
    \item{expansionValue}{Value that is added to each interaction interval end}
    \item{plotFileName}{Optional name for an output file}
    \item{picDim}{Dimensions of the plot}
}

\value{An RCircos-plot of trans interaction intervals}

\details{The code of Splinter et al to determine significant interactions provides chromosome, start and end of interaction intervals and a forth column with information on far-cis or trans data. This column is ignored by \code{plotTransInteractions}; it is assumed that all interactions for trans visualization are indeed trans interactions. Otherwise, far-cis interactions are visualized as well. While not a mistake per se, the (usually more numerous) far-cis interactions are easier to interpret if visualized with Splinter et al's spider-plot functions.}

\note{
PDF export and output as TIFF format are supported. The export format is chosen depending on the plot file name's ending. If no plot file name is provided, the result is plotted on screen.
}

\references{
Zhang, H., Meltzer, P. and Davis, S. (2013) RCircos: an R package for Circos 2D track plots, BMC Bioinformatics, 14, 244

Splinter, E., de Wit, E., van de Werken, H., et al. (2012) Determining long-range chromatin interactions for selected genomic sites using 4C-seq technology: From fixation to computation, Methods, 58, 221-230.
}

\examples{
    if(interactive()) {
        library(RCircos)
        interactions <- system.file("extdata", "transInteractionData.txt", package="Basic4Cseq")
        ideograms <- system.file("extdata", "RCircos_GRCm38_ideogram.csv", package="Basic4Cseq")
        plotTransInteractions(interactions, "10", c(20000042, 20001000), ideograms, PlotColor = "blue", expandBands = TRUE, expansionValue = 1000000, plotFileName = "")
    }
}

\author{Carolin Walter}

\keyword{plotTransInteractions}

