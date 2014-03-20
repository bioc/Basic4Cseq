# constructor for class "Data4Cseq"

setMethod("Data4Cseq", signature(viewpointChromosome="character", viewpointInterval="numeric", readLength="numeric", pointsOfInterest="data.frame", rawReads="AlignedRead"),
    function(viewpointChromosome, viewpointInterval, readLength, pointsOfInterest, rawReads) {
        
      for (i in 1:nrow(pointsOfInterest)) {
        if (pointsOfInterest[i,1] != viewpointChromosome) {
          message("One point of interest is not on the viewpoint chromosome and therefore not visible in a near-cis visualization. Alternatively, check for 'chrX / X' notation variety.")
        }        
      }
      pointsOfInterest = subset(pointsOfInterest, pointsOfInterest[,1] == viewpointChromosome)

      if (length(readLength) == 1) {
        if (readLength < 10) {
          message("The specified read length seems to be very short, please check if the number is correct")
        }
        if (readLength > 1000) {
          message("The specified read length seems to be very long, please check if the number is correct")
        }
      }

      if (viewpointInterval[1] > viewpointInterval[2]) {
        temp = viewpointInterval[1]
        viewpointInterval[1] = viewpointInterval[2]
        viewpointInterval[2] = temp
        message("The viewpoint coordinates (a,b) are expected to have the form a < b; switching coordinates...")
      }

      newData4CseqObject = new("Data4Cseq", viewpointChromosome = viewpointChromosome, viewpointInterval = viewpointInterval, readLength = readLength, pointsOfInterest = pointsOfInterest, rawReads = rawReads)
      
      return(newData4CseqObject)
    }
)


setMethod("Data4Cseq", signature(viewpointChromosome="character", viewpointInterval="numeric", readLength="numeric", pointsOfInterest="missing", rawReads="missing"),
    function(viewpointChromosome, viewpointInterval, readLength, pointsOfInterest, rawReads) {
        
      if (length(readLength) == 1) {
        if (readLength < 10) {
          message("The specified read length seems to be very short, please check if the number is correct")
        }
        if (readLength > 1000) {
          message("The specified read length seems to be very long, please check if the number is correct")
        }
      }

      if (viewpointInterval[1] > viewpointInterval[2]) {
        temp = viewpointInterval[1]
        viewpointInterval[1] = viewpointInterval[2]
        viewpointInterval[2] = temp
        message("The viewpoint coordinates (a,b) are expected to have the form a < b; switching coordinates...")
      }

      newData4CseqObject = new("Data4Cseq", viewpointChromosome = viewpointChromosome, viewpointInterval = viewpointInterval)
      
      return(newData4CseqObject)
    }
)
     

# "set" and "replace" methods for class "Data4Cseq"

setMethod("viewpointChromosome", signature(object="Data4Cseq"),
    function(object) {
        return(object@viewpointChromosome)
    }
)

setReplaceMethod("viewpointChromosome",
    signature=signature(object="Data4Cseq", value="character"),
    function(object, value) {
        object@viewpointChromosome = value
        return(object)
    }
)


setMethod("viewpointInterval", signature(object="Data4Cseq"),
    function(object) {
        return(object@viewpointInterval)
    }
)

setReplaceMethod("viewpointInterval",
    signature=signature(object="Data4Cseq", value="numeric"),
    function(object, value) {
        object@viewpointInterval = value
        return(object)
    }
)


setMethod("readLength", signature(object="Data4Cseq"),
    function(object) {
        return(object@readLength)
    }
)

setReplaceMethod("readLength",
    signature=signature(object="Data4Cseq", value="numeric"),
    function(object, value) {
        object@readLength = value
        return(object)
    }
)


setMethod("pointsOfInterest", signature(object="Data4Cseq"),
    function(object) {
        return(object@pointsOfInterest)
    }
)

setReplaceMethod("pointsOfInterest",
    signature=signature(object="Data4Cseq", value="data.frame"),
    function(object, value) {
        object@pointsOfInterest = value
        return(object)
    }
)


setMethod("rawReads", signature(object="Data4Cseq"),
    function(object) {
        return(object@rawReads)
    }
)

setReplaceMethod("rawReads",
    signature=signature(object="Data4Cseq", value="AlignedRead"),
    function(object, value) {
        object@rawReads = value
        return(object)
    }
)


setMethod("rawFragments", signature(object="Data4Cseq"),
    function(object) {
        return(object@rawFragments)
    }
)

setReplaceMethod("rawFragments",
    signature=signature(object="Data4Cseq", value="data.frame"),
    function(object, value) {
        object@rawFragments= value
        return(object)
    }
)

setMethod("nearCisFragments", signature(object="Data4Cseq"),
    function(object) {
        return(object@nearCisFragments)
    }
)

setReplaceMethod("nearCisFragments",
    signature=signature(object="Data4Cseq", value="data.frame"),
    function(object, value) {
        object@nearCisFragments= value
        return(object)
    }
)


setMethod("show", "Data4Cseq",
    function(object){
        cat("4C-seq experiment data\n")
        cat("Type:", class(object), "\n")
        cat("Viewpoint: ", viewpointChromosome(object), ":", viewpointInterval(object)[1], "-", viewpointInterval(object)[2], "\n")
        cat("Read length: ", readLength(object), "\n")
        cat("Number of reads: ", length(rawReads(object)), "\n")
        cat("Number of total fragments: ", nrow(rawFragments(object)), "\n")
        cat("Number of near-cis fragments: ", nrow(nearCisFragments(object)), "\n")
        cat("Points of interest: ", nrow(pointsOfInterest(object)), "\n")
    }
)
