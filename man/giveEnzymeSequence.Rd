\name{giveEnzymeSequence}

\alias{giveEnzymeSequence}
\alias{giveEnzymeSequence,character,character-method}

\title{Provide the corresponding enzyme sequence for an enzyme name}

\description{
This function is a small convenience function that reads in a prepared file with restriction enzyme sequence names and sequences. \code{giveEnzymeSequence} then provides restriction enzyme sequences for the example enzymes listed in van de Werken et al's 4Cseqpipe data base. 
}

\usage{giveEnzymeSequence(fileNameDatabase, enzymeName)}

\arguments{
  \item{fileNameDatabase}{File name of the prepared enzyme database}
  \item{enzymeName}{Name of the enzyme for which the sequence is to be returned}
}

\value{Character string with the restriction enzyme sequence}

\references{van de Werken, H., Landan, G., Holwerda, S., et al. (2012): Robust 4C-seq data analysis to screen for regulatory DNA interactions, Nature Methods, 9, 969-971.}

\note{For any custom-made enzyme list it is assumed that there are no duplicate enzyme names in the database.}

\examples{
  file <- system.file("extdata", "enzymeData.csv", package="Basic4Cseq")  
  giveEnzymeSequence(file, "NlaIII")
}

\author{Carolin Walter}

\keyword{giveEnzymeSequence}

