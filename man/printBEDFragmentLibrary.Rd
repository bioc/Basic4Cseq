\name{printBEDFragmentLibrary}

\alias{printBEDFragmentLibrary}
\alias{printBEDFragmentLibrary,character,character-method}

\title{Print a BED-file fragment library}

\description{
This function extracts the first columns of a virtual fragment library file and exports them as a BED-file for use with other tools (e.g. visualization in the Integrative Genomics Viewer (IGV).)
}

\usage{printBEDFragmentLibrary(fragmentLibrary, BEDLibraryName, minFragEndLength = 0, zeroBased = FALSE)
}

\arguments{
    \item{fragmentLibrary}{Virtual fragment library file name}
    \item{BEDLibraryName}{File name for the exported BED file}
    \item{minFragEndLength}{Minimum fragment end length to be considered}
    \item{zeroBased}{If TRUE, adapt the start of the BED-file fragments}
}

\value{writes BED-file containing the virtual fragment library position data}

\examples{
    if(interactive()) {
        file <- system.file("extdata", "vfl_aagctt_catg_mm9_54_vp.csv", package="Basic4Cseq")
        printBEDFragmentLibrary(file, "BEDLibrary_FL_vp.bed")
    }
}

\author{Carolin Walter}

\keyword{printBEDFragmentLibrary}
