# generics for data preprocessing
setGeneric("prepare4CseqData", signature=c("fastqFileName", "firstCutter", "fragmentLibrary", "referenceGenome"),
    function(fastqFileName, firstCutter, fragmentLibrary, referenceGenome, pathToBWA = "", pathToSam = "", pathToBED = "", controlCutterSequence = FALSE, bwaThreads = 1, minFragEndLength = 0)
        standardGeneric("prepare4CseqData"))


# generics for quality controls
setGeneric("getReadDistribution", signature=c("expData"),
    function(expData, distanceFromVP = 100000, useFragEnds = TRUE, outputName = "")
        standardGeneric("getReadDistribution"))


# generics for export / import functions
setGeneric("printBEDFragmentLibrary", signature=c("fragmentLibrary", "BEDLibraryName"),
    function(fragmentLibrary, BEDLibraryName, minFragEndLength = 0, zeroBased = FALSE)
        standardGeneric("printBEDFragmentLibrary"))

setGeneric("printWigFile", signature=c("expData"),
    function(expData, wigFileName = "output.wig", fixedSpan = TRUE, headerUCSC = "", useOnlyIndex = FALSE)
        standardGeneric("printWigFile"))

setGeneric("giveWigDataChromosome", signature=c("fragmentDataChromosome", "readLength", "chromosomeID"),
    function(fragmentDataChromosome, readLength, chromosomeID, fixedSpan = TRUE)
        standardGeneric("giveWigDataChromosome"))

setGeneric("readPointsOfInterestFile", signature=c("poiFile"),
    function(poiFile)
        standardGeneric("readPointsOfInterestFile"))

setGeneric("exportVisualizationFragmentData", signature=c("expData", "fileName"),
    function(expData, fileName, fullData = FALSE)
        standardGeneric("exportVisualizationFragmentData"))

setGeneric("importVisualizationFragmentData", signature=c("fileName"),
    function(fileName)
        standardGeneric("importVisualizationFragmentData"))


# generics for 4C-seq fragment analysis and helper functions
setGeneric("readsToFragments", signature=c("expData", "fragmentLib"),
    function(expData, fragmentLib)
        standardGeneric("readsToFragments"))

setGeneric("chooseNearCisFragments", signature=c("expData", "regionCoordinates"),
    function(expData, regionCoordinates, deleteViewpoint = TRUE)
        standardGeneric("chooseNearCisFragments"))

setGeneric("normalizeFragmentData", signature=c("expData"),
    function(expData)
        standardGeneric("normalizeFragmentData"))

setGeneric("formatFragmentData", signature=c("expData"),
    function(expData, useFragEnds = TRUE)
        standardGeneric("formatFragmentData"))

setGeneric("formatFragmentDataMain", signature=c("normalizedFragData"),
    function(normalizedFragData, useFragEnds = TRUE)
        standardGeneric("formatFragmentDataMain"))

setGeneric("giveEnzymeSequence", signature=c("fileNameDatabase", "enzymeName"),
    function(fileNameDatabase, enzymeName)
        standardGeneric("giveEnzymeSequence"))


# generics for 4C-seq visualization
setGeneric("visualizeViewpoint", signature=c("expData"),
    function(expData, poi = data.frame(chr = character(), start = character(), end = character(), name = character(), colour = character()), plotFileName = "", windowLength = 5, interpolationType = "median", picDim = c(9, 5), maxY = -1, minQuantile = 0.2, maxQuantile = 0.8, mainColour = "blue", plotTitle = "4C-seq plot", loessSpan = 0.1, xAxisIntervalLength = 50000, yAxisIntervalLength = 500, useFragEnds = TRUE)
        standardGeneric("visualizeViewpoint"))

setGeneric("visualizeViewpointMain", signature=c("fragsToVisualize"),
    function(fragsToVisualize, poi = data.frame(chr = character(), start = character(), end = character(), name = character(), colour = character()), plotFileName = "", windowLength = 5, interpolationType = "median", picDim = c(9, 5), maxY = -1, minQuantile = 0.2, maxQuantile = 0.8, mainColour = "blue", plotTitle = "4C-seq plot", loessSpan = 0.1, xAxisIntervalLength = 50000, yAxisIntervalLength = 500, useFragEnds = TRUE)
        standardGeneric("visualizeViewpointMain"))

setGeneric("basic4CPlot", signature=c("position", "yVal", "averageReads", "minQuantile", "maxQuantile", "windowLength", "colours"),
    function(position, yVal, averageReads, minQuantile, maxQuantile, windowLength, colours, loessSpan = 0.1)
        standardGeneric("basic4CPlot"))

setGeneric("drawMetaData", signature=c("minX", "maxX", "maxY"),
    function(minX, maxX, maxY, poi = data.frame(chr = character(), start = character(), end = character(), name = character(), colour = character()), xAxisIntervalLength = 50000, yAxisIntervalLength = 500)
        standardGeneric("drawMetaData"))

setGeneric("drawHeatmap", signature=c("expData"),
    function(expData, plotFileName = "", smoothingType = "median", picDim = c(9, 2.2), bands = 5, cutoffLog = -7.0, xAxisIntervalLength = 50000, legendLabels = expression(2^-7, 2^0), useFragEnds = TRUE)
        standardGeneric("drawHeatmap"))

setGeneric("drawHeatmapMain", signature=c("fragsToVisualize"),
    function(fragsToVisualize, plotFileName = "", smoothingType = "median", picDim = c(9, 2.2), bands = 5, cutoffLog = -7.0, xAxisIntervalLength = 50000, legendLabels = expression(2^-7, 2^0), useFragEnds = TRUE)
        standardGeneric("drawHeatmapMain"))

setGeneric("plotTransInteractions",
    function(interactionFile, chromosomeViewpoint, coordViewpoint, ideogramData, PlotColor = "default", expandBands = FALSE, expansionValue = 0, plotFileName = "", picDim = c(8, 8))
        standardGeneric("plotTransInteractions"))

setGeneric("plotTransInteractionsMain",
    function(interactionData, chromosomeViewpoint, coordViewpoint, cyto.info, PlotColor = "default", expandBands = FALSE, expansionValue = 0, plotFileName = "", picDim = c(8, 8))
        standardGeneric("plotTransInteractionsMain"))


# generics for first cutter check
setGeneric("checkRestrictionEnzymeSequence", signature=c("firstCutter", "inputFileName"),
    function(firstCutter, inputFileName, outputFileName = "output.sam", keepOnlyUniqueReads = TRUE, writeStatistics = TRUE)
        standardGeneric("checkRestrictionEnzymeSequence"))


# generics for fragment library creation
setGeneric("createVirtualFragmentLibrary", signature=c("chosenGenome", "firstCutter", "secondCutter", "readLength"),
    function(chosenGenome, firstCutter, secondCutter, readLength, onlyNonBlind = TRUE, useOnlyIndex = FALSE, minSize = 0, maxSize = -1, minFragEndSize = 0, maxFragEndSize = 10000000, useAllData = TRUE, chromosomeName = "chr1", libraryName = "default")
        standardGeneric("createVirtualFragmentLibrary"))

setGeneric("splitChromosome", signature=c("firstCutter", "secondCutter", "chromosomeToSplit", "chromosomeName"),
    function(firstCutter, secondCutter, chromosomeToSplit, chromosomeName, onlyNonBlind = TRUE)
        standardGeneric("splitChromosome"))

setGeneric("createVirtualFragmentLibraryMain", signature=c("totalFragments", "totalFragmentsRev", "firstCutter", "secondCutter", "readLength"),
    function(totalFragments, totalFragmentsRev, firstCutter, secondCutter, readLength, onlyNonBlind = TRUE, useOnlyIndex = FALSE, minSize = 0, maxSize = -1, minFragEndSize = 0, maxFragEndSize = 10000000, chromosomeName = "chr1", libraryName = "default")
        standardGeneric("createVirtualFragmentLibraryMain"))


# generics for in silico digestion ("real" fragment library control -> expected fragment sizes)
setGeneric("simulateDigestion", signature=c("firstCutter", "secondCutter", "dnaSequence"),
    function(firstCutter, secondCutter, dnaSequence)
        standardGeneric("simulateDigestion"))

setGeneric("simulateDigestionChromosome", signature=c("firstCutter", "secondCutter", "dnaSequence"),
    function(firstCutter, secondCutter, dnaSequence)
        standardGeneric("simulateDigestionChromosome"))

setGeneric("drawDigestionFragmentHistogram", signature=c("fragments"),
    function(fragments, minLength = 0, maxLength = 10000)
        standardGeneric("drawDigestionFragmentHistogram"))


# class generics
setGeneric("Data4Cseq", function(viewpointChromosome, viewpointInterval, readLength, pointsOfInterest, rawReads) {
    standardGeneric("Data4Cseq")})

setGeneric("viewpointChromosome", function(object) {
    standardGeneric("viewpointChromosome")})

setGeneric("viewpointChromosome<-", function(object, value) {
    standardGeneric("viewpointChromosome<-")})

setGeneric("viewpointInterval", function(object) {
    standardGeneric("viewpointInterval")})

setGeneric("viewpointInterval<-", function(object, value) {
    standardGeneric("viewpointInterval<-")})

setGeneric("readLength", function(object) {
    standardGeneric("readLength")})

setGeneric("readLength<-", function(object, value) {
    standardGeneric("readLength<-")})

setGeneric("pointsOfInterest", function(object) {
    standardGeneric("pointsOfInterest")})

setGeneric("pointsOfInterest<-", function(object, value) {
    standardGeneric("pointsOfInterest<-")})

setGeneric("rawReads", function(object) {
    standardGeneric("rawReads")})

setGeneric("rawReads<-", function(object, value) {
    standardGeneric("rawReads<-")})

setGeneric("rawFragments", function(object) {
    standardGeneric("rawFragments")})

setGeneric("rawFragments<-", function(object, value) {
    standardGeneric("rawFragments<-")})

setGeneric("nearCisFragments", function(object) {
    standardGeneric("nearCisFragments")})

setGeneric("nearCisFragments<-", function(object, value) {
    standardGeneric("nearCisFragments<-")})
