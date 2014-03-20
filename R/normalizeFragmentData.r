.normalizeFragmentData <- function(expData) {
  
  readCount = length(rawReads(expData))
  relevantFragData = nearCisFragments(expData)

  # simple RPM normalization (reads per million)
  RPMFrags <- relevantFragData
  RPMFrags$leftFragEndReads <- RPMFrags$leftFragEndReads / readCount * 1000000
  RPMFrags$rightFragEndReads <- RPMFrags$rightFragEndReads / readCount * 1000000
  RPMFrags$fragEndReadsAverage <- RPMFrags$fragEndReadsAverage / readCount * 1000000

  return(RPMFrags)
}




setMethod("normalizeFragmentData",
    signature=signature(expData="Data4Cseq"),
    .normalizeFragmentData)
