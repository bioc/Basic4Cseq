\name{visualizeViewpoint}

\alias{visualizeViewpoint}
\alias{visualizeViewpoint,data.frame-method}
\alias{visualizeViewpoint,Data4Cseq-method}

\title{Draw a near-cis coverage plot for 4C-seq data}

\description{
    This method creates a plot of near-cis 4C-seq fragment data around the experiment's viewpoint. Fragment-based raw data is visualized as grey dots, interpolated data (running median / running mean) as coloured dots. Trend line and quantiles are loess-smoothed; the trend line is shown as colored line whereas the quantiles are depicted as light-grey bands. A corresponding quantile legend is added in an extra plot.
}

\usage{
visualizeViewpoint(expData, poi = data.frame(chr = character(), start = character(), end = character(), name = character(), colour = character()), plotFileName = "", windowLength = 5, interpolationType = "median", picDim = c(9, 5), maxY = -1, minQuantile = 0.2, maxQuantile = 0.8, mainColour = "blue", plotTitle = "4C-seq plot", loessSpan = 0.1, xAxisIntervalLength = 50000, yAxisIntervalLength = 500, useFragEnds = TRUE)
}


\arguments{
    \item{expData}{Experiment data of class \code{Data4Cseq} with information on the 4C-seq experiment, including normalized near-cis fragment data for visualization}
    \item{poi}{Points of interest that will be marked in the plot}
    \item{plotFileName}{Name for the 4C-seq plot file}
    \item{windowLength}{Length of the window for running median / running mean that is used to smooth the trend line}
    \item{interpolationType}{Type of interpolation, either running median or running mean}
    \item{picDim}{Dimensions of the plot}
    \item{maxY}{Maximum y-value to plot. If no maximum is given, the maximum running median / mean value is used}  
    \item{minQuantile}{Minimum quantile to draw}
    \item{maxQuantile}{Maximum quantile to draw}
    \item{mainColour}{Main colour of the plot}
    \item{plotTitle}{Title of the 4C-seq plot, depicted above the main plot}
    \item{loessSpan}{Span value for the loess curve; smaller values mean a tighter fit to the data points, but a value that is too small may produce errors}
    \item{xAxisIntervalLength}{Length of the x axis intervals in the plot}
    \item{yAxisIntervalLength}{Length of the y axis intervals in the plot}
    \item{useFragEnds}{Indicates whether fragment end data is used directly or interpolated on fragment level}
}

\value{A near-cis coverage plot and a corresponding quantile legend}

\note{
PDF export and output as TIFF format are supported. The export format is chosen depending on the plot file name's ending. If no plot file name is provided, the result is plotted on screen.
}

\examples{
    data(liverData)
    file <- system.file("extdata", "fetalLiverVP.bed", package="Basic4Cseq")
        visualizeViewpoint(liverData, readPointsOfInterestFile(file), plotFileName = "", mainColour = "red", plotTitle = "Fetal Liver Near-Cis Plot", loessSpan = 0.1, maxY = 6000, xAxisIntervalLength = 50000, yAxisIntervalLength = 1000)
}

\author{Carolin Walter}

\keyword{visualizeViewpoint}
