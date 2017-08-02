.readPointsOfInterestFile <- function(poiFile) {
    
    if (poiFile != "") {
        poiTable = read.table(poiFile, sep = "\t", header = FALSE)
        colnames(poiTable) = c("chr", "start", "end", "name", "colour")

        chrNames = unique(poiTable$chr)
        if (length(chrNames) > 1) {
            message("Please provide only points of interest on the viewpoint chromosome for the near-cis visualizations")
        }
    } else {
        poiTable = data.frame(chr = character(), start = character(), end = character(), name = character(), colour = character())
    }

    return(poiTable)
}


.exportVisualizationFragmentData <- function(expData, fileName, fullData = FALSE) {
    
    if (fullData) {
        fragmentData = nearCisFragments(expData)    
    } else {
        fragmentData = formatFragmentData(expData)
    } 
    write.table(fragmentData, file = fileName, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}


.importVisualizationFragmentData <- function(fileName) {

    fragmentData = read.table(file = fileName, header = TRUE)
    return(fragmentData)
}




setMethod("readPointsOfInterestFile",
    signature=signature(poiFile="character"),
    .readPointsOfInterestFile)


setMethod("exportVisualizationFragmentData",
    signature=signature(expData="Data4Cseq", fileName="character"),
    .exportVisualizationFragmentData)


setMethod("importVisualizationFragmentData",
    signature=signature(fileName="character"),
    .importVisualizationFragmentData)
