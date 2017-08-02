.formatFragmentData_Data4Cseq <- function(expData, useFragEnds = TRUE) {

    normalizedFragData = nearCisFragments(expData)

    formatFragmentDataMain(normalizedFragData, useFragEnds)
}


.formatFragmentData_DataFrame <- function(expData, useFragEnds = TRUE) {

    normalizedFragData = expData

    formatFragmentDataMain(normalizedFragData, useFragEnds)
}


.formatFragmentDataMain <- function(normalizedFragData, useFragEnds = TRUE) {

    if (useFragEnds) {
        fragmentDataLeftUnique = subset(normalizedFragData, normalizedFragData$leftFragEndValid == TRUE)
        fragmentDataRightUnique = subset(normalizedFragData, normalizedFragData$rightFragEndValid == TRUE)

        fragsLeft = fragmentDataLeftUnique[,c(1,2,4,11)]
        fragsRight = fragmentDataRightUnique[,c(1,4,3,12)]
        fragsRight[,2] = fragsRight[,2] + 1

        colnames(fragsLeft) = c("chrom", "start", "end", "reads")
        colnames(fragsRight) = c("chrom", "start", "end", "reads")

        fragEnds = rbind(fragsLeft, fragsRight)
        fragmentDataFinal = fragEnds[order(fragEnds[,1], fragEnds[,2]),]

    } else {
        fragments = subset(normalizedFragData, normalizedFragData$leftFragEndValid == TRUE & normalizedFragData$rightFragEndValid == TRUE)

        fragmentDataFinal = fragments[,c(1:3,13)]
        colnames(fragmentDataFinal) = c("chrom", "start", "end", "reads")
    }

    row.names(fragmentDataFinal) = NULL

    return(fragmentDataFinal)
}




setMethod("formatFragmentData",
    signature=signature(expData="Data4Cseq"),
    .formatFragmentData_Data4Cseq)

setMethod("formatFragmentData",
    signature=signature(expData="data.frame"),
    .formatFragmentData_DataFrame)

setMethod("formatFragmentDataMain",
    signature=signature(normalizedFragData="data.frame"),
    .formatFragmentDataMain)
