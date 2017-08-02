.printBEDFragmentLibrary <- function(fragmentLibrary, BEDLibraryName, minFragEndLength = 0, zeroBased = FALSE) {

    allFragments = read.table(fragmentLibrary, header = TRUE)
    
    # relevant fragments: all fragments with at least one relevant (i.e. unique and of sufficient length) frag end 
    relevantFragments = subset(allFragments, ((allFragments$leftFragEndLength >= minFragEndLength) & (allFragments$leftFragEndValid == TRUE)) | ((allFragments$rightFragEndLength >= minFragEndLength) & (allFragments$rightFragEndValid == TRUE)) )

    if (zeroBased) {
        relevantFragments$fragmentStart = relevantFragments$fragmentStart - 1
    }

    write.table(relevantFragments[,1:3], file = BEDLibraryName, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}




setMethod("printBEDFragmentLibrary",
    signature=signature(fragmentLibrary="character", BEDLibraryName="character"),
    .printBEDFragmentLibrary)
