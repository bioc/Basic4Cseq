\name{simulateDigestion}

\alias{simulateDigestion}
\alias{simulateDigestion,character,character,character-method}
\alias{simulateDigestion,character,character,BSgenome-method}

\title{Simulate the digestion of a genome}

\description{
This function simulates the digestion process with two restriction enzymes for a dna sequence or a BSgenome package. The information can then be used for quality controls of the biological fragment library.
}

\usage{
    simulateDigestion(firstCutter, secondCutter, dnaSequence)
}

\arguments{
    \item{firstCutter}{First restriction enzyme sequence for the digestion process}
    \item{secondCutter}{Second restriction enzyme sequence for the digestion process}
    \item{dnaSequence}{DNA sequence that is digested}
}

\value{Data frame with lengths and corresponding frequences of fragments}

\details{The resulting virtual library of fragment parts does not provide information on blind or non-blind fragments, but provides information on the fragment length distribution of the real (i.e. biological) 4C-seq library. In contrast to the regular virtual fragment library for 4C-seq data, fragments between two adjacent secondary restriction enzyme sites are counted as well.}

\note{
    \itemize{
        \item The resulting fragment lengths and corresponding frequencies can easily be visualized with R's plot function or the small convenience function \code{drawDigestionFragmentHistogram}
        \item The resulting table of fragment frequencies can easily be exported with R's \code{write.table} function
    }
}

\examples{
    shortTestGenome = "ATCCATGTAGGCTAAGTACACATGTTAAGGTACAGTACAATTGCACGATCAT"
    fragments = simulateDigestion("catg", "gtac", shortTestGenome)
    head(fragments)
}

\author{Carolin Walter}

\keyword{simulateDigestion}

