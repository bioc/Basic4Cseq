.getReadDistribution <- function(expData, distanceFromVP = 100000, useFragEnds = TRUE, outputName = "") {

    coordinatesVP = viewpointInterval(expData)

    fragmentData = subset(rawFragments(expData), rawFragments(expData)$chromosomeName == viewpointChromosome(expData) & (rawFragments(expData)$fragmentEnd < coordinatesVP[1] | rawFragments(expData)$fragmentStart > coordinatesVP[2]))
    vpAreaFragments = subset(fragmentData, (fragmentData$fragmentStart >= coordinatesVP[1]-distanceFromVP & fragmentData$fragmentEnd <= coordinatesVP[2]+distanceFromVP))

    vpAreaFragments <- formatFragmentData(vpAreaFragments, useFragEnds)

    # filter invalid reads
    readsTotal = subset(rawReads(expData), !(is.na(start(rawReads(expData)))))

    # get cis reads
    readsVPChrom = subset(readsTotal, as.character(seqnames(readsTotal)) == viewpointChromosome(expData))

    # prepare region for viewpoint region read control
    maxLeft = coordinatesVP[1] - distanceFromVP
    maxRight = coordinatesVP[2] + distanceFromVP

    readsVPRegion = subset(readsVPChrom, start(readsVPChrom) >= maxLeft & end(readsVPChrom) <= maxRight)

    percentageChrom = round(length(readsVPChrom) / length(readsTotal) * 100, digits = 2)
    percentageVP = round(length(readsVPRegion) / length(readsTotal) * 100, digits = 2)

    coveredReadFragments = subset(vpAreaFragments, vpAreaFragments$reads != 0)

    percentageCov = round(nrow(coveredReadFragments) / nrow(vpAreaFragments) * 100, digits = 2)

    readsTotalText = paste("total reads:", length(readsTotal), sep = " ")
    readsVPChromText = paste("reads on the viewpoint chromosome: ", length(readsVPChrom), " (", percentageChrom, "% of total reads)", sep = "")
    readsVPRegionText = paste("reads in the viewpoint region: ", length(readsVPRegion), " (", percentageVP, "% of total reads)", sep = "")
    coverageVPRegionText = paste("covered fragment ends in the viewpoint region: ", percentageCov, "%", sep = "")

    statistics = paste(readsTotalText, readsVPChromText, readsVPRegionText, coverageVPRegionText, sep = "\n")

    if (outputName == "") {
        print(readsTotalText)
        print(readsVPChromText)
        print(readsVPRegionText)
        print(coverageVPRegionText)
    } else {
        fileConn = file(outputName)
        writeLines(statistics, fileConn)
        close(fileConn)
    }
}




setMethod("getReadDistribution",
    signature=signature(expData="Data4Cseq"),
    .getReadDistribution)
