\name{createVirtualFragmentLibrary}

\alias{createVirtualFragmentLibrary}
\alias{createVirtualFragmentLibrary,BSgenome,character,character,numeric-method}
\alias{createVirtualFragmentLibrary,DNAString,character,character,numeric-method}

\title{Create a virtual fragment library from a provided genome and two restriction enzymes}

\description{\pkg{Basic4Cseq} can create virtual fragment libraries from any BSgenome package or DNAString object. Two restriction enzymes have to be specified to cut the DNA, the read length is needed to check the fragment ends of corresponding length for uniqueness. Filter options (minimum and maximum size) are provided on fragment level and on fragment end level.}

\usage{
createVirtualFragmentLibrary(chosenGenome, firstCutter, secondCutter, readLength, onlyNonBlind = TRUE, useOnlyIndex = TRUE, minSize = 0, maxSize = -1, minFragEndSize = 0, maxFragEndSize = 10000000, useAllData = TRUE, chromosomeName = "chr1", libraryName = "default")
}

\arguments{
 \item{chosenGenome}{The genome that is to be digested in silico with the provided enzymes; can be an instance of BSgenome or DNAString}
 \item{firstCutter}{First of two restriction enzymes}
 \item{secondCutter}{Second of two restriction enzymes}
 \item{readLength}{Read length for the experiment}
 \item{onlyNonBlind}{Variable that is TRUE (default) if only non-blind fragments are considered (i.e. all blind fragments are removed)}
 \item{useOnlyIndex}{Convenience function to adapt the annotation style of the chromosomes ("chr1", \dots "chrY" or "1", \dots, "Y"); parameter depends on personal preferences}
 \item{minSize}{Filter option that allows to delete fragments below a certain size (in bp)}
 \item{maxSize}{Filter option that allows to delete fragments above a certain size (in bp)}
 \item{minFragEndSize}{Filter option that allows to delete fragment ends below a certain size (in bp)}
 \item{maxFragEndSize}{Filter option that allows to delete fragment ends above a certain size (in bp)}
 \item{useAllData}{Variable that indicates if all data of a BSgenome package is to be used. If FALSE, chromosome names including a "_" are removed, reducing the set of chromosomes to (1 \dots 19, X, Y, MT) for the mouse genome or (1 \dots 22, X, Y, MT) for the human genome}
 \item{chromosomeName}{Chromosome name for the virtual fragment library if a \code{DNAString} object is used instead of a \code{BSgenome} object.}
 \item{libraryName}{Name of the file the created virtual fragment library is written to. Per default the file is called "fragments_firstCutter_secondCutter.csv". The fragment data is returned as a data frame if and only if an empty character string is chosen as \code{libraryName}.}
}

\details{
  \itemize{
    \item \code{readLength} is relevant for the creation of the virtual fragment library to differenciate between unique and non-unique fragment ends. While two fragments can be unique, their respective ends may be repetitive if only the first few bases are considered. For 4C-seq data, reads can only map to the start (or end, respectively) of a 4C-seq fragment, the remaining fragment part is not covered. The length of a fragment end that has to be checked for uniqueness therefore depends on the read length of the experiment.
    \item \code{useAllData} uses the lengths of the chromosomes to identify relevant ones, based on the current BSgenome packages for mm10 or hg19, and may therefore provide undesirable results for smaller genomes with different lengths (i.e. discard all chromosomes).
     \item The length of a fragment influences the expected read count of a 4C-seq fragment. Per default, \pkg{Basic4Cseq} uses the experiment's read length as minimum fragment end size and places virtually no limit on the maximum fragment end size.
  }
}

\value{A tab-separated file with the specified virtual fragment library (containing fragment position, length, presence of second restriction enzyme and uniqueness of the fragment ends)}

\note{
  \itemize{
    \item It is strongly recommended to preprocess and store the virtual fragment library if a number of experiments with the same restriction enzyme combination, read length and underlying genome have to be analyzed.
    \item Processing one of the larger BSgenome packages takes some time and computer data storage.
    \item If no library name for the virtual fragment library is specified, the fragment data is returned as a data frame. If the library name "default" is chosen, the tab-separated file is named "fragments_firstCutter_secondCutter" (with variable cutter sequences).
  }
}

\examples{
  if(interactive()) {
    library(BSgenome.Ecoli.NCBI.20080805)
    fragmentData = createVirtualFragmentLibrary(chosenGenome = Ecoli$NC_002655, firstCutter = "catg", secondCutter = "gtac", readLength = 30,  onlyNonBlind = TRUE, chromosomeName = "NC_002655", libraryName = "fragments_Ecoli.csv")
  }
}

\author{Carolin Walter}

\keyword{createVirtualFragmentLibrary}
