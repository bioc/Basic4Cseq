\name{drawHeatmap}

\alias{drawHeatmap}
\alias{drawHeatmap,Data4Cseq-method}
\alias{drawHeatmap,data.frame-method}

\title{Draw a heatmap-like multi-scale contact profile}

\description{
This method draws a fragment-based heatmap-like plot for 4C-seq data around a given viewpoint. For a given number of bands, color-coded running medians or running means of signal intensity (normalized and log-scaled) in different fragments are displayed; the window size of the running medians or running means increases from top to bottom. A corresponding colour legend is added in an extra plot.
}

\usage{
  \S4method{drawHeatmap}{Data4Cseq}(expData, plotFileName = "", smoothingType = "median", picDim = c(9, 2.2), bands = 5, cutoffLog = -7.0, xAxisIntervalLength = 50000, legendLabels = expression(2^-7, 2^0), useFragEnds = TRUE)
  \S4method{drawHeatmap}{data.frame}(expData, plotFileName = "", smoothingType = "median", picDim = c(9, 2.2), bands = 5, cutoffLog = -7.0, xAxisIntervalLength = 50000, legendLabels = expression(2^-7, 2^0), useFragEnds = TRUE)
}

\arguments{
  \item{expData}{Experiment data from a given 4C-seq experiment for visualization; can be a \code{Data4Cseq} object or a data frame}
  \item{plotFileName}{Name for the heatmap plot}
  \item{smoothingType}{Type of interpolation (running mean or running median). Default value is "median" (i.e. running median)}
  \item{picDim}{Dimensions of the plot. Default value is c(9, 2.2), to fit a small heatmap plot below the main 4C-seq plot that is created by \code{visualizeViewpoint}}
  \item{bands}{Number of coloured "bands" (rows) to visualize. The first band contains the raw data (running median or running mean with window size 1), the following bands increase in window size (+2 per band)}
  \item{cutoffLog}{Cut off value for the logarithmic scale}
  \item{xAxisIntervalLength}{Length of the x axis intervals in the plot}
  \item{legendLabels}{Labels for a heat colour legend plot; labels should correspond to the logarithmic cut offs}
  \item{useFragEnds}{Indicates whether fragment end data is used directly or interpolated on fragment level}
}

\value{A multiscale intensity contact profile plot and a corresponding colour legend)}

\note{
PDF export and output as TIFF format are supported. The export format is chosen depending on the plot file name's ending. If no plot file name is provided, the result is plotted on screen.
}

\examples{
  if(interactive()) {
    data(liverData)
    drawHeatmap(liverData)
  }
}

\author{Carolin Walter}

\keyword{drawHeatmap}
