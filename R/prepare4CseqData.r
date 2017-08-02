.prepare4CseqData <- function(fastqFileName, firstCutter, fragmentLibrary, referenceGenome, pathToBWA = "", pathToSam = "", pathToBED = "", controlCutterSequence = FALSE, bwaThreads = 1, minFragEndLength = 0) {

    tempAnalysisID = unlist(strsplit(fastqFileName, "/"))
    analysisID = sub(".fastq", "", tempAnalysisID[length(tempAnalysisID)])

    # align the raw 4C-seq reads with BWA
    system(paste(pathToBWA, "bwa aln -n 2 -t ", bwaThreads, " ", referenceGenome, " ", fastqFileName, " > ", analysisID, ".sai", sep = ""))
    system(paste(pathToBWA, "bwa samse ", referenceGenome, " ", analysisID, ".sai ", fastqFileName, " > ", analysisID, "_raw.sam", sep = ""))

    system(paste("rm ", analysisID, ".sai", sep = ""))

    if (controlCutterSequence) {
        # filter: check cutter sequence for mismatches
        checkRestrictionEnzymeSequence(firstCutter, paste(analysisID, "_raw.sam", sep = ""), paste(analysisID, "_filtered.sam", sep = ""))
        longAnalysisID = paste(analysisID, "_filtered", sep = "")
        #system(paste("rm ", analysisID, "_raw.sam", sep = ""))
    } else {
        longAnalysisID = paste(analysisID, "_raw", sep = "")
    }

    # prepare sorted bam-file and index for possible IGV visualization 
    system(paste(pathToSam, "samtools view -b -S ", longAnalysisID, ".sam > ", longAnalysisID, ".bam", sep = ""))
    system(paste(pathToSam, "samtools sort ", longAnalysisID, ".bam ", longAnalysisID, sep = ""))
    system(paste(pathToSam, "samtools index ", longAnalysisID, ".bam", sep = ""))

    # extract relevant rows from [virtual] fragment library for intersectBed
    printBEDFragmentLibrary(fragmentLibrary, "fragmentLibBED.bed", minFragEndLength = minFragEndLength, zeroBased = FALSE)

    # intersect filtered 4C-seq reads with fragment library
    system(paste(pathToBED, "intersectBed -abam ", longAnalysisID, ".bam -b fragmentLibBED.bed > ", analysisID, "_4Cseq.bam", sep = ""))

    # remove bed file (full information stored in fragment library)
    system("rm fragmentLibBED.bed")

    # sort and index remaining 4C-seq reads for possible IGV visualization
    system(paste(pathToSam, "samtools sort ", analysisID, "_4Cseq.bam ", analysisID, "_4Cseq", sep = ""))
    system(paste(pathToSam, "samtools index ", analysisID, "_4Cseq.bam", sep = ""))
}




setMethod("prepare4CseqData",
    signature=signature(fastqFileName="character", firstCutter="character", fragmentLibrary="character", referenceGenome="character"),
    .prepare4CseqData)
