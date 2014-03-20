.chooseNearCisFragments <- function(expData, regionCoordinates, deleteViewpoint = TRUE) {

  fragmentData = rawFragments(expData)
  
  # for cis-interactions: only consider viewpoint chromosome
  fragmentData <- subset(fragmentData, fragmentData$chromosomeName == viewpointChromosome(expData))
  
  # delete viewpoint fragments, if chosen (default: TRUE)
  if (deleteViewpoint) {
    fragmentData = subset(fragmentData, (fragmentData$fragmentEnd < viewpointInterval(expData)[1] | fragmentData$fragmentStart > viewpointInterval(expData)[2]))
  }

  # pick relevant fragments for visualization range
  relevantFragments = (subset(fragmentData, (fragmentData$fragmentStart >= regionCoordinates[1] & fragmentData$fragmentEnd <= regionCoordinates[2])))

  return(relevantFragments)  
}




setMethod("chooseNearCisFragments",
    signature=signature(expData="Data4Cseq", regionCoordinates="numeric"),
    .chooseNearCisFragments)
