.readsToFragments <- function(expData, fragmentLib) {

  readsTotal = rawReads(expData)
  # check all chromosomes where reads are present
  chromosomeNames = unique(as.vector(seqnames(readsTotal)))

  # read fragment data (with  meta data, like uniqueness of frag-ends + second cutter site)
  fragmentTableTotal = read.table(fragmentLib, header = TRUE)
  
  fragmentDataTotal = NULL
  
  for (i in 1:length(chromosomeNames)) {

    reads = subset(readsTotal, seqnames(readsTotal) == chromosomeNames[i])

    readsPlus = subset(reads, strand(reads) == "+")
    readsMinus = subset(reads, strand(reads) == "-")

    if (length(readsPlus) > 0 & length(readsMinus) > 0) {

      baseCovPlus = as.numeric(coverage(readsPlus)[[chromosomeNames[i]]])
      baseCovMinus = as.numeric(coverage(readsMinus)[[chromosomeNames[i]]])

      fragmentTable = subset(fragmentTableTotal, fragmentTableTotal$chromosomeName == chromosomeNames[i])

      # add zeroes if last fragment is not covered
      tempPlus = fragmentTable[nrow(fragmentTable),3] - length(baseCovPlus)
      if (tempPlus < 0) {
        tempPlus = 0
      }
      baseCovPlus = c(baseCovPlus,  rep(0, times=tempPlus))
      tempMinus = fragmentTable[nrow(fragmentTable),3] - length(baseCovMinus)
      if (tempMinus < 0) {
        tempMinus = 0
      }
      baseCovMinus = c(baseCovMinus,  rep(0, times=tempMinus))

      # calculate read count for fragment start, end, and average
      leftFragEndReads = baseCovPlus[fragmentTable$fragmentStart]
      rightFragEndReads = baseCovMinus[fragmentTable$fragmentEnd]
      fragEndReadsAverage = round((leftFragEndReads + rightFragEndReads) / 2)

      fragmentData = data.frame(fragmentTable, leftFragEndReads, rightFragEndReads, fragEndReadsAverage)
      fragmentDataTotal = rbind(fragmentDataTotal, fragmentData)
    }
  }

  return(fragmentDataTotal)  
}




setMethod("readsToFragments",
    signature=signature(expData="Data4Cseq", fragmentLib="character"),
    .readsToFragments)
