.giveEnzymeSequence <- function(fileNameDatabase, enzymeName) {

    enzymeDB = read.table(fileNameDatabase, header = TRUE)
    chosenEnzyme = subset(enzymeDB, enzymeDB$name == enzymeName)
    enzymeSequence = as.character(chosenEnzyme[1,2])  

    return(enzymeSequence)  
}




setMethod("giveEnzymeSequence",
    signature=signature(fileNameDatabase="character", enzymeName="character"),
    .giveEnzymeSequence)
