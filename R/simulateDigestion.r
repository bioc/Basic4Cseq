.simulateDigestion_Character <- function(firstCutter, secondCutter, dnaSequence) {

    fragmentSequencesFirstDigestion = unlist(strsplit(dnaSequence, split=toupper(firstCutter)))
    fragmentSequencesSecondDigestion = unlist(strsplit(fragmentSequencesFirstDigestion, split=toupper(secondCutter)))
    fragmentLengths = as.data.frame(table(nchar(fragmentSequencesSecondDigestion)))
    colnames(fragmentLengths) = c("length", "frequency")
    return(fragmentLengths)
}


.simulateDigestion_BSgenome <- function(firstCutter, secondCutter, dnaSequence) {

    chromosomeNames = seqnames(dnaSequence)
    totalFragments = NULL

    for (i in 1:length(chromosomeNames)) {
        chromosomeToSplit = dnaSequence[[i]]
        currentChromosomeFragments = simulateDigestionChromosome(firstCutter, secondCutter, chromosomeToSplit)
        totalFragments = c(totalFragments, currentChromosomeFragments)
    }

    fragmentLengths = as.data.frame(table(nchar(totalFragments)))
    colnames(fragmentLengths) = c("length", "frequency")
    return(fragmentLengths)
}


.simulateDigestionChromosome_MaskedDNAString <- function(firstCutter, secondCutter, dnaSequence) {

    dnaSequence = toString(as(unmasked(dnaSequence), "Views"))
    fragmentSequencesFirstDigestion = unlist(strsplit(dnaSequence, split=toupper(firstCutter)))
    fragmentSequencesSecondDigestion = unlist(strsplit(fragmentSequencesFirstDigestion, split=toupper(secondCutter)))

    return(fragmentSequencesSecondDigestion)
}


.simulateDigestionChromosome_DNAString <- function(firstCutter, secondCutter, dnaSequence) {

    dnaSequence = toString(dnaSequence)
    fragmentSequencesFirstDigestion = unlist(strsplit(dnaSequence, split=toupper(firstCutter)))
    fragmentSequencesSecondDigestion = unlist(strsplit(fragmentSequencesFirstDigestion, split=toupper(secondCutter)))

    return(fragmentSequencesSecondDigestion)
}


.drawDigestionFragmentHistogram <- function(fragments, minLength = 0, maxLength = 10000) {

    fragmentData = subset(fragments, as.numeric(fragments$length) >= minLength & as.numeric(fragments$length) <= maxLength)
    plot(as.vector(fragmentData$length), as.vector(fragmentData$frequency), type = "h", main = "Fragment frequencies", lty = 1, xlab = "length", ylab = "frequency")
}


setMethod("simulateDigestion",
    signature=signature(firstCutter="character", secondCutter="character", dnaSequence="character"),
    .simulateDigestion_Character)

setMethod("simulateDigestion",
    signature=signature(firstCutter="character", secondCutter="character", dnaSequence="BSgenome"),
    .simulateDigestion_BSgenome)

setMethod("simulateDigestionChromosome",
    signature=signature(firstCutter="character", secondCutter="character", dnaSequence="MaskedDNAString"),
    .simulateDigestionChromosome_MaskedDNAString)

setMethod("simulateDigestionChromosome",
    signature=signature(firstCutter="character", secondCutter="character", dnaSequence="DNAString"),
    .simulateDigestionChromosome_DNAString)

setMethod("drawDigestionFragmentHistogram",
    signature=signature(fragments="data.frame"),
    .drawDigestionFragmentHistogram)
