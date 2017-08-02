.plotTransInteractions_Character <- function(interactionFile, chromosomeViewpoint, coordViewpoint, ideogramData, PlotColor = "default", expandBands = FALSE, expansionValue = 0, plotFileName = "", picDim = c(8, 8)) {

    interactionData = read.table(interactionFile)
    cyto.info = read.table(ideogramData, header = TRUE)

    plotTransInteractionsMain(interactionData, chromosomeViewpoint, coordViewpoint, cyto.info, PlotColor, expandBands, expansionValue, plotFileName, picDim)
}


.plotTransInteractions_DataFrame <- function(interactionFile, chromosomeViewpoint, coordViewpoint, ideogramData, PlotColor = "default", expandBands = FALSE, expansionValue = 0, plotFileName = "", picDim = c(8, 8)) {

    interactionData = interactionFile
    cyto.info = ideogramData

    plotTransInteractionsMain(interactionData, chromosomeViewpoint, coordViewpoint, cyto.info, PlotColor, expandBands, expansionValue, plotFileName, picDim)
}


.plotTransInteractionsMain <- function(interactionData, chromosomeViewpoint, coordViewpoint, cyto.info, PlotColor = "default", expandBands = FALSE, expansionValue = 0, plotFileName = "", picDim = c(8, 8)) {

    chromA = chromosomeViewpoint
    chromStartA = coordViewpoint[1]
    chromEndA = coordViewpoint[2]
    chromB = as.vector(interactionData[,1])
    chromStartB = interactionData[,2]
    chromEndB = interactionData[,3]

    if (expandBands) {
        chromStartA = max(0, chromStartA - expansionValue)
        chromEndA = chromEndA + expansionValue
        for (i in 1:length(chromB)) {
            chromStartB[i] = max(0, chromStartB[i] - expansionValue)
            chromEndB[i] = chromEndB[i] + expansionValue
        }
    }

    if (PlotColor != "default") {
        circosData = data.frame(chromA, chromStartA, chromEndA, chromB, chromStartB, chromEndB, PlotColor)
    } else {
        circosData = data.frame(chromA, chromStartA, chromEndA, chromB, chromStartB, chromEndB)
    }
    chr.exclude = NULL

    tracks.inside = 1
    tracks.outside = 0
    RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

    if (plotFileName != "") {
        if (grepl(".pdf", plotFileName)) {
            pdf(file = plotFileName, title = plotFileName, width = picDim[1], height = picDim[2], useDingbats = FALSE)  
        } else {
            tiff(filename = plotFileName, width = picDim[1], height = picDim[2], compression = "none", bg = "white", pointsize = 20)
        }
    }

    RCircos.Set.Plot.Area()

    plot.window(c(-1.5,1.5), c(-1.5,1.5))

    RCircos.Chromosome.Ideogram.Plot()

    RCircos.Ribbon.Plot(ribbon.data=circosData, track.num=1, by.chromosome=FALSE, twist=FALSE)

    if (plotFileName != "") {
        dev.off()
    }
}




setMethod("plotTransInteractions",
    signature=signature(interactionFile="character", chromosomeViewpoint="character", coordViewpoint="numeric", ideogramData="character"),
    .plotTransInteractions_Character)

setMethod("plotTransInteractions",
    signature=signature(interactionFile="data.frame", chromosomeViewpoint="character", coordViewpoint="numeric", ideogramData="data.frame"),
    .plotTransInteractions_DataFrame)

setMethod("plotTransInteractionsMain",
    signature=signature(interactionData="data.frame", chromosomeViewpoint="character", coordViewpoint="numeric", cyto.info="data.frame"),
    .plotTransInteractionsMain)


