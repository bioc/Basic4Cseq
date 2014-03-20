.checkRestrictionEnzymeSequence <- function(firstCutter, inputFileName, outputFileName = "output.sam", keepOnlyUniqueReads = TRUE, writeStatistics = TRUE) {

  cutterLength = nchar(firstCutter)
  
  # read lines from specified sam-file
  sam = readLines(inputFileName)

  # initialize statistic variables
  unmapped = 0
  mismatchStart = 0
  mismatchEnd = 0
  incompleteCutter = 0
  totalReads = length(sam)
  mapped = 0
  notUnique = 0

  revCompFirstCutter = toString(reverseComplement(DNAString(firstCutter)))

  for (i in 1:length(sam)) {

    temp = sam[i]

    # check each line if it is a part of the header and use only non-header lines
    if (substr(temp, start = 1, stop = 1) == "@") {
      totalReads = totalReads - 1
    } else {

      # for each line in the sam file that does not belong to the header, split at "TAB"
      split = unlist(strsplit(temp, "\\\t"))

      # get read length and save flag entry
      readLength = nchar(split[10])
      flag = split[2]
 
      # for now, leave reads without "XT:A:" tag and mark only those with XT:A != U for deletion
      if (keepOnlyUniqueReads) {
    
        splitUnique = unlist(strsplit(temp, "XT:A:"))
        isNotUnique = FALSE
    
        if (length(splitUnique) == 2) {
          if (substr(splitUnique[2], start = 1, stop = 1) != "U") {
            # not a unique read
            isNotUnique = TRUE
          }
        }
      }
  
      # identify unmapped reads via flag-entry and mark for deletion
      if (intToBits(as.numeric(flag))[3] == 01) {
        unmapped = unmapped + 1
        sam[i] = "toDelete"    
      } else if (keepOnlyUniqueReads == TRUE & isNotUnique == TRUE) {
        notUnique = notUnique + 1
        sam[i] = "toDelete"    
      } else { 
        # mapped read (both strands):
        # 1. check if cutter sequence is complete
        # 2. check if cigar string shows correct mapping (absense of insertion / deletions)  of cutter sequence
        # 3. check if MD string shows correct mapping (absense of SNPs) of cutter sequence

        seq = split[10]
        seqStart = substr(seq, start=1, stop=cutterLength)
        seqEnd = substr(seq, start=readLength-cutterLength+1, stop=readLength)
   
        if (intToBits(as.numeric(flag))[5] == 00)  {
          # case a): read maps to forward strand

          # check if cutter sequence is intact
          if (seqStart == toupper(firstCutter)) {    
      
            startCigar = unlist(strsplit(split[6], "M"))

            if (nchar(startCigar[1]) < (cutterLength-1)) {
              # start of string has form XXM or XM -->  no I, D or extended cigar symbols possible
              if (as.numeric(startCigar[1]) < cutterLength) {
                # less than (lengthCutter) matches --> cutter sequence not properly mapped --> delete read

                mismatchStart = mismatchStart + 1
                sam[i] = "toDelete"
              } else {
                # enough matches and complete cutter sequence at start of read (cigar string) --> check MD string for possible SNPs 

                MDpre = unlist(strsplit(temp, "MD:Z:"))
                MD = unlist(strsplit(MDpre[2], "\\\t"))[1]
  
                MDA = substr(MD[1], start = 1, stop = 1)
                MDB = substr(MD[1], start = 2, stop = 2)

                if ((MDB == "A") || (MDB == "G") || (MDB == "C") || (MDB == "T")) {
                  # start of MD has form XA, XG, XC or XD (X in 0:9) --> check X (number of matches) if at least (cutterLength) matches are present 
 
                  if (as.numeric(MDA) < cutterLength) {

                    mismatchStart = mismatchStart + 1
                    sam[i] = "toDelete"
                  } else {
                    # cutter sequence complete
                    mapped = mapped + 1
                  }
                } else {
                  # start of MD has form XX (X in 0:9) --> number of matches >= 10 --> enough matches
                  mapped = mapped + 1
                }
                
              }
            } else {
              # insertions, deletions or extended cigar symbols at start of string --> cutter sequence not properly mapped --> delete read

              mismatchStart = mismatchStart + 1
              sam[i] = "toDelete"          
            }

          } else {
            # cutter sequence is not intact
        
            incompleteCutter = incompleteCutter + 1
            sam[i] = "toDelete"
          }

        } else {
          # case b): read maps to reverse strand

          if (seqEnd == revCompFirstCutter) {
    
            # check end of cigar string: at least (cutterLength; 4 or 6) matches required
            lCig = nchar(split[6])
            cigA = substr(split[6], start = lCig-2, stop = lCig-2)
            cigB = substr(split[6], start = lCig-1, stop = lCig-1)
            cigC = substr(split[6], start = lCig, stop = lCig)
  
            if (cigC != "M") {
              # end of string has insertions or deletions or extended cigar symbols --> cutter sequence not properly mapped --> delete read
              sam[i] = "toDelete"

              mismatchEnd = mismatchEnd + 1
        
            } else if ((cigA == "I" || cigA == "D" || cigA == "N" || cigA == "S" || cigA == "H" || cigA == "P") && (as.numeric(cigB) < cutterLength) ){
              # number of matches at the end of the string is < cutterLength --> cutter sequence not properly mapped --> delete read
              sam[i] = "toDelete"

              mismatchEnd = mismatchEnd + 1
            } else {
              # >= (cutterLength) mapped bases at the end of the read, complete cutter sequence --> check MD string for possible SNPs

              MDpre = unlist(strsplit(temp, "MD:Z:"))
              MD = unlist(strsplit(MDpre[2], "\\\t"))[1]

              MDA = substr(MD[1], start = readLength-1, stop = readLength-1)
              MDB = substr(MD[1], start = readLength, stop = readLength)
              
              if ((MDA == "A") || (MDA == "G") || (MDA == "C") || (MDA == "T")) {
                # end of MD has form AX, GX, CX or DX (X in 0:9) --> check X (number of matches) if at least four matches are present

                if (as.numeric(MDB) < cutterLength) {
                  mismatchEnd = mismatchEnd + 1
                  sam[i] = "toDelete"
                } else {
                  # cutter sequence complete
                  mapped = mapped + 1
                }
              } else {
                # end of MD has form XX (X in 0:9) --> number of matches >= 10 --> enough matches
                mapped = mapped + 1            
              }
              
            }      
            
          } else {
            # cutter sequence not complete: delete read

            incompleteCutter = incompleteCutter + 1
            sam[i] = "toDelete"        
          }      
          
        }
      }
    }
  }
 
  # optional: write statistics to additional output file
  if (writeStatistics) {
    totalStats = paste("total reads:", totalReads, sep = " ")
    mappedStats = paste("mapped reads:", mapped, sep = " ")
    unmappedStats = paste("unmapped reads:", unmapped, sep = " ")
    incompleteCutterStats = paste("incomplete cutter:", incompleteCutter, sep = " ")
    mismatchStats = paste("mismatched cutter:", mismatchStart + mismatchEnd, sep = " ")
    notUniqueStats = paste("not unique read:", notUnique, sep = " ")
  
    statistics = paste(inputFileName, totalStats, mappedStats, unmappedStats, incompleteCutterStats, mismatchStats, notUniqueStats, sep = "\n")
  
    fileConn = file(sub(".sam", "_statistics.txt", outputFileName))
    writeLines(statistics, fileConn)
    close(fileConn)
  }  
  
  # filter reads that are marked for deletion
  samNew =  sam[!sam == "toDelete"]
  
  # write remaining data to specified file
  fileConn = file(outputFileName)
  writeLines(samNew, fileConn)
  close(fileConn)
}




setMethod("checkRestrictionEnzymeSequence",
    signature=signature(firstCutter="character", inputFileName="character"),
    .checkRestrictionEnzymeSequence)
