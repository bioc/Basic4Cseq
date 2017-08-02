.printWigFile <- function(expData, wigFileName = "output.wig", fixedSpan = TRUE, headerUCSC = "", useOnlyIndex = FALSE) {
    
    fragmentData = rawFragments(expData)
    chromosomes = unique(as.vector(fragmentData$chromosomeName))

    # prepare chromosome names for wig --> even if "1"..."Y" have been used for the library or in the reference genome, 
    # one may wish to print the wig with "chr1"..."chrY" notations
    chromosomeNamesWig = gsub("chr", "", as.character(chromosomes))
    if (!useOnlyIndex) {
        chromosomeNamesWig = paste("chr", chromosomeNamesWig, sep = "")
    }

    wigFileContent = NULL

    if (headerUCSC != "") {
        wigFileContent = headerUCSC
    }
    
    for (i in 1:length(chromosomes)) {

        fragsChromosome = subset(fragmentData, fragmentData$chromosomeName == chromosomes[i])
        chromosomeWig = giveWigDataChromosome(fragsChromosome, readLength(expData), chromosomeNamesWig[i], fixedSpan)

        wigFileContent = c(wigFileContent, chromosomeWig)
    }
    
    fileConn <- file(wigFileName)

    writeLines(wigFileContent, fileConn)

    close(fileConn)
}


.giveWigDataChromosome <- function(fragmentDataChromosome, readLength, chromosomeID, fixedSpan = TRUE) {

    if (!fixedSpan) {
        # variable step
        wigHeader = paste("variableStep chrom=", chromosomeID, sep = "")
    } else {
        # variable step with fixed span (read length)
        wigHeader = paste("variableStep chrom=", chromosomeID, " span=", readLength, sep = "")
    }

    # get reads on valid unique fragments
    fragmentDataLeftValid = subset(fragmentDataChromosome, fragmentDataChromosome$leftFragEndValid == TRUE)
    fragmentDataRightValid = subset(fragmentDataChromosome, fragmentDataChromosome$rightFragEndValid == TRUE)

    fragsLeft = cbind(fragmentDataLeftValid$fragmentStart, fragmentDataLeftValid$leftFragEndReads)
    fragsRight = cbind(fragmentDataRightValid$fragmentEnd - readLength + 1, fragmentDataRightValid$rightFragEndReads)

    if (!fixedSpan) {
        fragsLeftPos = NULL
        fragsLeftReads = NULL
        fragsRightPos = NULL
        fragsRightReads = NULL
        for (i in 1:readLength) {
            fragsLeftPos = c(fragsLeftPos, fragsLeft[,1]+i-1)
            fragsLeftReads = c(fragsLeftReads, fragsLeft[,2])
            fragsRightPos = c(fragsRightPos, fragsRight[,1]+i-1)
            fragsRightReads = c(fragsRightReads, fragsRight[,2])
        }
        fragsLeft = cbind(fragsLeftPos, fragsLeftReads)
        fragsRight = cbind(fragsRightPos, fragsRightReads)
    }

    fragEnds = rbind(fragsLeft, fragsRight)
    fragEnds = fragEnds[order(fragEnds[,1]),]

    wigData = c(wigHeader, paste(fragEnds[,1], fragEnds[,2], sep = " "))

    return(wigData)
}




setMethod("printWigFile",
    signature=signature(expData="Data4Cseq"),
    .printWigFile)


setMethod("giveWigDataChromosome",
    signature=signature(fragmentDataChromosome="data.frame", readLength="numeric", chromosomeID="character"),
    .giveWigDataChromosome)
