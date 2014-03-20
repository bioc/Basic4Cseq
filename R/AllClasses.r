setClass("Data4Cseq",
  representation=list(           
    viewpointChromosome="character",
    viewpointInterval="numeric",
    readLength="numeric",
    pointsOfInterest="data.frame",
    rawReads="AlignedRead",
    rawFragments="data.frame",
    nearCisFragments="data.frame"
  ),         
         
  prototype=prototype(
    viewpointChromosome="chr1",
    viewpointInterval=c(1,2),
    readLength=30,
    pointsOfInterest=data.frame(),
    rawReads=AlignedRead(),
    rawFragments=data.frame(),
    nearCisFragments=data.frame()
  ),

  validity=function(object) {
    if (length(readLength(object)) != 1)
      return("Read length has to be a single number")
    if (length(viewpointInterval(object)) != 2)
      return("The viewpoint interval should consist of a start and end coordinate of the form c(A, B)")
    return(TRUE)
  }
)





