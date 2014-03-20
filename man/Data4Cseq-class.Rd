\name{Data4Cseq-class}

\docType{class}
\alias{Data4Cseq-class}
\alias{viewpointChromosome}
\alias{viewpointChromosome,Data4Cseq-method}
\alias{viewpointChromosome<-}
\alias{viewpointChromosome<-,Data4Cseq,character-method}
\alias{viewpointInterval}
\alias{viewpointInterval,Data4Cseq-method}
\alias{viewpointInterval<-}
\alias{viewpointInterval<-,Data4Cseq,numeric-method}
\alias{readLength}
\alias{readLength,Data4Cseq-method}
\alias{readLength<-}
\alias{readLength<-,Data4Cseq,numeric-method}
\alias{pointsOfInterest}
\alias{pointsOfInterest,Data4Cseq-method}
\alias{pointsOfInterest<-}
\alias{pointsOfInterest<-,Data4Cseq,data.frame-method}
\alias{rawReads}
\alias{rawReads,Data4Cseq-method}
\alias{rawReads<-}
\alias{rawReads<-,Data4Cseq,AlignedRead-method}
\alias{rawFragments}
\alias{rawFragments,Data4Cseq-method}
\alias{rawFragments<-}
\alias{rawFragments<-,Data4Cseq,data.frame-method}
\alias{nearCisFragments}
\alias{nearCisFragments,Data4Cseq-method}
\alias{nearCisFragments<-}
\alias{nearCisFragments<-,Data4Cseq,data.frame-method}

\title{Class \code{"Data4Cseq"}}
\description{
This class is a container for information on a specific 4C-seq experiment. Stored information includes raw reads, fragment data, and the experiment's viewpoint location.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Data4Cseq", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{viewpointChromosome}:}{Object of class \code{"character"} representing the viewpoint chromosome's name}
    \item{\code{viewpointInterval}:}{Object of class \code{"numeric"} representing the viewpoint interval's location}
    \item{\code{readLength}:}{Object of class \code{"numeric"} representing the experiment's read length}
    \item{\code{pointsOfInterest}:}{Object of class \code{"data.frame"} representing any points of interest to be marked in the near-cis visualizations}
    \item{\code{rawReads}:}{Object of class \code{"AlignedRead"} representing the raw 4C-seq reads of the experiment}
    \item{\code{rawFragments}:}{Object of class \code{"data.frame"} representing the experiment's corresponding virtual fragment library}
    \item{\code{nearCisFragments}:}{Object of class \code{"data.frame"} representing near-cis data in fragment form}
  }
}

\section{Methods}{
  \describe{
    \item{viewpointChromosome<-}{\code{signature(object = "Data4Cseq", value = "character")}: 
      Setter-method for the viewpointChromosome slot.}
    \item{viewpointChromosome}{\code{signature(object = "Data4Cseq")}: 
      Getter-method for the viewpointChromosome slot.}
    \item{viewpointInterval<-}{\code{signature(object = "Data4Cseq", value = "numeric")}: 
      Setter-method for the viewpointInterval slot.}
    \item{viewpointInterval}{\code{signature(object = "Data4Cseq")}: 
      Getter-method for the viewpointInterval slot.}
    \item{readLength<-}{\code{signature(object = "Data4Cseq", value = "numeric")}: 
      Setter-method for the readLength slot.}
    \item{readLength}{\code{signature(object = "Data4Cseq")}: 
      Getter-method for the readLength slot.}
    \item{pointsOfInterest<-}{\code{signature(object = "Data4Cseq", value = "data.frame")}: 
      Setter-method for the pointsOfInterest slot.}
    \item{pointsOfInterest}{\code{signature(object = "Data4Cseq")}: 
      Getter-method for the pointsOfInterest slot.}
    \item{rawReads<-}{\code{signature(object = "Data4Cseq", value = "AlignedRead")}: 
      Setter-method for the rawReads slot.}
    \item{rawReads}{\code{signature(object = "Data4Cseq")}: 
      Getter-method for the rawReads slot.}
    \item{rawFragments<-}{\code{signature(object = "Data4Cseq", value = "data.frame")}: 
      Setter-method for the rawFragments slot.}
    \item{rawFragments}{\code{signature(object = "Data4Cseq")}: 
      Getter-method for the rawFragments slot.}
    \item{nearCisFragments<-}{\code{signature(object = "Data4Cseq", value = "data.frame")}: 
      Setter-method for the nearCisFragments slot.}
    \item{nearCisFragments}{\code{signature(object = "Data4Cseq")}: 
      Getter-method for the nearCisFragments slot.}
  }
}

\author{
Carolin Walter
}

\examples{
showClass("Data4Cseq")
}

\keyword{classes}
