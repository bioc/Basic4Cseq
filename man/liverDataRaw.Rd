\name{liverDataRaw}

\alias{liverDataRaw}

\docType{data}

\title{
  Example 4C-seq data set of fetal liver data 
}

\description{
  This data set contains an instance of a \code{Data4Cseq} object; 2175 reads on 453 fragments are included. Raw reads are mapped to fragments, but the read count has not yet been normalized.
  
  The 4C-seq data was taken from Stadhouders et al's fetal liver data set.
}

\usage{
  data("liverDataRaw")
}

\format{
  Formal class 'Data4Cseq'
}

\source{
  Shortened version of Stadhouders et al's fetal liver data:
  
  Stadhouders, R., Thongjuea, S., et al. (2012): Dynamic long-range chromatin interactions control Myb proto-oncogene transcription during erythroid development. EMBO, 31, 986-999.
}

\examples{
  data("liverDataRaw")
  liverDataRaw
}

\keyword{datasets}
