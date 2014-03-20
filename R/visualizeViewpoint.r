.visualizeViewpoint_Data4Cseq <- function(expData, poi = data.frame(chr = character(), start = character(), end = character(), name = character(), colour = character()), plotFileName = "", windowLength = 5, interpolationType = "median", picDim = c(9, 5), maxY = -1, minQuantile = 0.2, maxQuantile = 0.8, mainColour = "blue", plotTitle = "4C-seq plot", loessSpan = 0.1, xAxisIntervalLength = 50000, yAxisIntervalLength = 500, useFragEnds = TRUE) {
 
  fragsToVisualize = formatFragmentData(expData, useFragEnds)
  poi = pointsOfInterest(expData)

  visualizeViewpointMain(fragsToVisualize, poi, plotFileName, windowLength, interpolationType, picDim, maxY, minQuantile, maxQuantile, mainColour, plotTitle, loessSpan, xAxisIntervalLength, yAxisIntervalLength, useFragEnds)
}


.visualizeViewpoint_DataFrame <- function(expData, poi = data.frame(chr = character(), start = character(), end = character(), name = character(), colour = character()), plotFileName = "", windowLength = 5, interpolationType = "median", picDim = c(9, 5), maxY = -1, minQuantile = 0.2, maxQuantile = 0.8, mainColour = "blue", plotTitle = "4C-seq plot", loessSpan = 0.1, xAxisIntervalLength = 50000, yAxisIntervalLength = 500, useFragEnds = TRUE) {
 
  fragsToVisualize = expData
  
  if (!useFragEnds) {
    print("visualizeViewpoint: use of fragment-wise interpolated data not possible for data frame input")
  }

  visualizeViewpointMain(fragsToVisualize, poi, plotFileName, windowLength, interpolationType, picDim, maxY, minQuantile, maxQuantile, mainColour, plotTitle, loessSpan, xAxisIntervalLength, yAxisIntervalLength, useFragEnds)
}


.visualizeViewpointMain <- function(fragsToVisualize, poi = data.frame(chr = character(), start = character(), end = character(), name = character(), colour = character()), plotFileName = "", windowLength = 5, interpolationType = "median", picDim = c(9, 5), maxY = -1, minQuantile = 0.2, maxQuantile = 0.8, mainColour = "blue", plotTitle = "4C-seq plot", loessSpan = 0.1, xAxisIntervalLength = 50000, yAxisIntervalLength = 500, useFragEnds = TRUE) {

  # prepare different color shades  
  colours = c(paste(mainColour, "", sep = ""), paste(mainColour, "2", sep = ""), paste(mainColour, "3", sep = ""))
  
  # prepare data
  yVal = 0
  
  position = round((fragsToVisualize$start + fragsToVisualize$end) / 2)
  averageReads = fragsToVisualize$reads 
  
  if (interpolationType == "median") {
    # smoothed data (median with defined window length)
    yVal = runmed(averageReads, windowLength, endrule = "constant")
    
  } else if (interpolationType == "mean") {
    # smoothed data (mean with defined window length)
    yVal = runmean(averageReads, windowLength, endrule = "constant")
    
  } else {
    # raw data
    yVal = averageReads    
  }
  
  # if no (reasonable) maximum y-value is provided, use the maximum of yValA and yValB, rounded up to the next y-axis marker
  if (maxY == -1) {
    maxY = ceiling(max(yVal) / 500) * 500
  }
  
  # prepare plot
  if (plotFileName != "") {
    if (grepl(".pdf", plotFileName)) {
      pdf(file = plotFileName, title = plotTitle, width = picDim[1], height = picDim[2], useDingbats = FALSE)  
    } else {
      tiff(filename = plotFileName, width = picDim[1], height = picDim[2], compression = "none", bg = "white", pointsize = 20)
    }
  }
  
  # basic plot
  plot(position, yVal,  ylim = c(0,maxY), ylab = "fragment read count (RPM)", xlab = "fragment position", type = "n", main = plotTitle, lty = 1, axes = FALSE)
  
  # plot sample data
  basic4CPlot(position, yVal, averageReads, minQuantile, maxQuantile, windowLength, colours, loessSpan) 
  
  drawMetaData(fragsToVisualize$start[1], fragsToVisualize$end[nrow(fragsToVisualize)], maxY, poi, xAxisIntervalLength, yAxisIntervalLength)   
  
  if (plotFileName != "") {
    dev.off()
  }
  
  # prepare legend plot
  scale = 100
  
  if (plotFileName != "") {
    if (grepl(".pdf", plotFileName)) {
      pdf(file = "quantile_legend.pdf", title = "Main Trend", width = 2, height = 5, useDingbats = FALSE)  
    } else {
      tiff(filename = "quantile_legend.tiff", width = 200, height = 500, compression = "none", bg = "white", pointsize = 20)
    }
    
    # basic plot
    plot(c(0,10), c(0,1), type='n', main="Main trend", axes = FALSE, xlab = '', ylab = "loess-smoothed quantiles")
  
    # add axis...
    axis(2, c(0, 0.5, 1), labels = c(paste(minQuantile, "%", sep = ""), "50%", paste(maxQuantile, "%", sep = "")), las = 1)
  
    # ... quantile representation...
    for (i in 1:100) {    
      y = (i-1) / scale
      rect(0, 0, 10, 1, col = "grey90", border=NA)
    }
  
    # ... and trend line
    lines(c(0, 10), c(0.5, 0.5), col = "black", lwd = 3)
  
    dev.off()
  }
}


.basic4CPlot <- function(position, yVal, averageReads, minQuantile, maxQuantile, windowLength, colours, loessSpan = 0.1) {

  # calculate quantiles and smooth the curve with R's loess function
  quanMin = runquantile(averageReads, windowLength, probs = minQuantile)
  yMin <- loess(quanMin ~ position, span=loessSpan, data.frame(position = position, quanMin = quanMin))
  yMinPredict <- predict(yMin, data.frame(position = position))
    
  quanMax = runquantile(averageReads, windowLength, probs = maxQuantile)
  yMax <- loess(quanMax ~ position, span=loessSpan, data.frame(position = position, quanMax = quanMax))
  yMaxPredict <- predict(yMax, data.frame(position = position))

  # visualize quantiles as polygon (light-grey)
  polygon(c(position, rev(position)), c(yMinPredict, rev(yMaxPredict)), col = "grey90", border = NA)
 
  # plot interpolated data points (median or mean)
  lines(position, yVal, type="p", col= colours[1], lwd = 0.2, lty = 1, pch = 20, cex = 0.5)
  
  # smooth main trend (loess)
  yLoess <- loess(yVal ~ position, span = loessSpan, data.frame(position = position, yVal = yVal))
  yPredict <- predict(yLoess, data.frame(position = position))
  lines(position, yPredict, col = colours[3], lwd = 1)

  # raw data points (fragment-based)
  points(position, averageReads, col = "grey40", pch = ".")  
}


.drawMetaData <- function(minX, maxX, maxY, poi = data.frame(chr = character(), start = character(), end = character(), name = character(), colour = character()), xAxisIntervalLength = 50000, yAxisIntervalLength = 500) {
  
  # if present, print points of interest as defined by the user
  if (nrow(poi) > 0) {
    pointsX = c(poi$start, poi$end)
    pointsY = round(maxY * 0.97)

    points(round((poi$start + poi$end) / 2), rep(round(maxY * 0.95), times = nrow(poi)), cex = 1.0, pch = 25, col = as.vector(poi$colour), bg = as.vector(poi$colour))
    text(round((poi$start + poi$end) / 2), round(maxY * 0.99), poi$name, cex = 1.0)
  }

  # print axes with custom intervals
  minXAxis = floor(minX / xAxisIntervalLength) * xAxisIntervalLength
  maxXAxis = ceiling(maxX / xAxisIntervalLength) * xAxisIntervalLength
  seqXAxis = seq(minXAxis, maxXAxis, by = xAxisIntervalLength)
  axis(side = 1, at = seqXAxis, labels = seqXAxis, las = 0)
  
  seqYAxis = seq(0, ceiling(maxY / yAxisIntervalLength) * yAxisIntervalLength, by = yAxisIntervalLength)
  axis(side = 2, at = seqYAxis, labels = seqYAxis, las = 2)  
}




setMethod("visualizeViewpoint",
    signature=signature(expData="Data4Cseq"),
    .visualizeViewpoint_Data4Cseq)

setMethod("visualizeViewpoint",
    signature=signature(expData="data.frame"),
    .visualizeViewpoint_DataFrame)

setMethod("visualizeViewpointMain",
    signature=signature(fragsToVisualize="data.frame"),
    .visualizeViewpointMain)


setMethod("basic4CPlot",
    signature=signature(position="numeric", yVal="numeric", averageReads="numeric", minQuantile="numeric", maxQuantile="numeric", windowLength="numeric", colours="character"),
    .basic4CPlot)


setMethod("drawMetaData",
    signature=signature(minX="numeric", maxX="numeric", maxY="numeric"),
    .drawMetaData)
