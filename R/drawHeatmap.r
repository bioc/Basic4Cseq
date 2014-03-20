.drawHeatmap_Data4Cseq <- function(expData, plotFileName = "", smoothingType = "median", picDim = c(9, 2.2), bands = 5, cutoffLog = -7.0, xAxisIntervalLength = 50000, legendLabels = expression(2^-7, 2^0), useFragEnds = TRUE) {

  fragsToVisualize = formatFragmentData(expData, useFragEnds)

  drawHeatmapMain(fragsToVisualize, plotFileName, smoothingType, picDim, bands, cutoffLog, xAxisIntervalLength, legendLabels, useFragEnds)
}


.drawHeatmap_DataFrame <- function(expData, plotFileName = "", smoothingType = "median", picDim = c(9, 2.2), bands = 5, cutoffLog = -7.0, xAxisIntervalLength = 50000, legendLabels = expression(2^-7, 2^0), useFragEnds = TRUE) {

  fragsToVisualize = expData
  
  if (!useFragEnds) {
    print("drawHeatmap: use of fragment-wise interpolated data not possible for data frame input")
  }

  drawHeatmapMain(fragsToVisualize, plotFileName, smoothingType, picDim, bands, cutoffLog, xAxisIntervalLength, legendLabels, useFragEnds)
}
  

.drawHeatmapMain <- function(fragsToVisualize, plotFileName = "", smoothingType = "median", picDim = c(9, 2.2), bands = 5, cutoffLog = -7.0, xAxisIntervalLength = 50000, legendLabels = expression(2^-7, 2^0), useFragEnds = TRUE) {

  position = round((fragsToVisualize$start + fragsToVisualize$end) / 2)
  averageReads = fragsToVisualize$reads
  
  # prepare heatmap-like data (varied window lengths for running mean / running median)
  signal = matrix(0, bands, nrow(fragsToVisualize))
  
  for (i in 1:bands) {
    
    wl = i * 2 - 1
    
    if (smoothingType == "median") {
      # smoothed data (median with defined window length)
      signal[i,] = runmed(averageReads, wl, endrule = "keep")
      
    } else {
      # smoothed data (mean with defined window length)
      signal[i,] = runmean(averageReads, wl, endrule = "keep")
    }     
  }
  
  # prepare heatmap plot
  if (plotFileName != "") {
    if (grepl(".pdf", plotFileName)) {
      pdf(file = plotFileName, title = plotFileName, width = picDim[1], height = picDim[2], useDingbats = FALSE)  
    } else {
      tiff(filename = plotFileName, width = picDim[1], height = picDim[2], compression = "none", bg = "white", pointsize = 20)
    }
  }
  
  plot(position, averageReads,  ylim = c(0,1), ylab = "wl", xlab = "fragment position", type = "n", lty = 1, axes = FALSE)
  
  scale = 1 / (bands + 1.5)
  maxSignal = max(signal)
  
  # log-scale signal
  signal = log2(signal/maxSignal)
  
  # set cut-offs for log-scaled data
  capMin = cutoffLog
  capMax = 0.0
  
  # define axes with custom intervals
  minXAxis = floor(fragsToVisualize$start[1] / xAxisIntervalLength) * xAxisIntervalLength
  maxXAxis = ceiling(fragsToVisualize$end[nrow(fragsToVisualize)] / xAxisIntervalLength) * xAxisIntervalLength
  seqXAxis = seq(minXAxis, maxXAxis, by = xAxisIntervalLength)
  axis(side = 1, at = seqXAxis, labels = seqXAxis, las = 0)
  
  seqYAxis = c(0.5, bands - 0.5) * scale
  labelsYAxis = c(2*bands-1,1)    
  
  axis(side = 2, at = seqYAxis, labels = labelsYAxis, las = 2)
  
  # mark unused fragments (blind / repeated)
  rect(minXAxis, 0, maxXAxis, bands * scale, col = "black", border = NA)
  
  # prepare color palette 
  colorNumber = capMax*100-capMin*100 + 1
  colors = heat.colors(colorNumber)
  
  # pick colours for the fragments
  for (i in 1:bands) {
    
    for (j in 1:nrow(fragsToVisualize)) {
      
      if (signal[i,j] > capMax) {
        signal[i,j] = capMax
      }
      if (signal[i,j] < capMin) {
        signal[i,j] = capMin
      }
      
      index = round((signal[i,j] - capMin) * 99) + 1
      
      pickedColor = colors[index]
      
      # mark 'missing' fragments at the sides (only fully present running median / mean windows are used)
      if ((j < i) | (j > (nrow(fragsToVisualize) - i))) {
        pickedColor = "black"
      }
      
      # print fragments in chosen colours
      rect(fragsToVisualize[,2][j], (bands-i+1)*scale, fragsToVisualize[,3][j], (bands-i)*scale, col = pickedColor, border = NA)
      
    }
  }
  
  if (plotFileName != "") {
    dev.off()
  }
  
  # prepare heatmap legend
  chosenPalette = heat.colors(100)
  
  # prepare color legend plot
  if (plotFileName != "") {
    if (grepl(".pdf", plotFileName)) {
      pdf(file = "color_legend.pdf", title = "Color Legend", width = 2, height = 5, useDingbats=FALSE)  
    } else {
      tiff(filename = "color_legend.tiff", width = 200, height = 500, compression = "none", bg = "white", pointsize = 20)
    }
    
    par(mar=c(5.1,4.1,4.1,4.1))
  
    # basic plot
    plot(c(0,10), c(0,1), type='n', main="Signal intensity", axes = FALSE, xlab = "", ylab = "relative window coverage (log-scaled)")
  
    # add axis...
    legendLength = length(legendLabels) - 1
    axis(2, c((0:legendLength)/legendLength), labels = legendLabels, las = 1)
  
    # ... and color bands
    for (i in 1:(length(chosenPalette)-1)) {    
      y = (i-1) / 100
      rect(0, y, 10, y+1/100, col = chosenPalette[i], border=NA)
    }
  
    dev.off()
  }
}




setMethod("drawHeatmap",
    signature=signature(expData="Data4Cseq"),
    .drawHeatmap_Data4Cseq)

setMethod("drawHeatmap",
    signature=signature(expData="data.frame"),
    .drawHeatmap_DataFrame)

setMethod("drawHeatmapMain",
    signature=signature(fragsToVisualize="data.frame"),
    .drawHeatmapMain)
