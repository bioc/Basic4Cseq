\name{printWigFile}

\alias{printWigFile}
\alias{printWigFile,Data4Cseq-method}

\title{Print a wig file from 4C-seq read data}

\description{
This function provides wig files from filtered fragment data. Only reads on unique frag-ends are considered for the export. Export of wig files with a fixed span length requires a uniform read length throughout the data.

While some tools (e.g. the Integrative Genomics Viewer, IGV) accept 'raw' wig data, the UCSC browser needs a header line for correct visualizations. A basic header line has the form 'track type=wiggle_0', but may also contain information on the track's name and a short description. Since the header line may complicate possible downstream analysis of the wig files, no header is included per default.   
}

\usage{printWigFile(expData, wigFileName = "output.wig", fixedSpan = TRUE, headerUCSC = "", useOnlyIndex = FALSE)
}

\arguments{
    \item{expData}{Experiment data of class \code{Data4Cseq} with information on the 4C-seq experiment}
    \item{wigFileName}{Name of the wig file that is written to hard disk}
    \item{fixedSpan}{If TRUE, use a fixed span for the wig file}
    \item{headerUCSC}{A header line for the UCSC browser}
    \item{useOnlyIndex}{If TRUE, use only '1,2,...Y' as chromosome names, if FALSE, use 'chr1,chr2...chrY'.}
}

\value{A wig file containing the experiment's reads}

\references{
UCSC Genome Browser: Kent WJ, Sugnet CW, Furey TS, Roskin KM, Pringle TH, Zahler AM, Haussler D. The human genome browser at UCSC. Genome Res. 2002 Jun;12(6):996-1006. 

http://genome.ucsc.edu/
}

\examples{
    if(interactive()) {
        data(liverData)
        printWigFile(liverData, wigFileName = "fetalLiver.wig")
    }
}

\author{Carolin Walter}

\keyword{printWigFile}

