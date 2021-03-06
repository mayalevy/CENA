% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BinsConstructor.R
\docType{class}
\name{LeafBinNode-class}
\alias{LeafBinNode-class}
\alias{LeafBinNode}
\title{Reference Class to represent an internal bin in the bins tree}
\description{
Reference Class to represent an internal bin in the bins tree
}
\section{Fields}{

\describe{
\item{\code{left_son}}{the left son of the internal bin}

\item{\code{right_son}}{the right son of the internal bin}

\item{\code{midpoint}}{the value that seperate the bin to left and right bins}

\item{\code{discriminator_dim}}{the dimension by which the bin was separated}

\item{\code{count_diff}}{the different points number between the left and the right son}

\item{\code{points}}{a matrix of the points in this bin}
}}

\section{Methods}{

\describe{
\item{\code{split(dim)}}{split the bin into 2 sub-bins in the simension in which the count_diff is the largest}
}}

\keyword{Class}
\keyword{Reference}
\keyword{an}
\keyword{bin}
\keyword{bins}
\keyword{in}
\keyword{internal}
\keyword{leaf}
\keyword{represent}
\keyword{the}
\keyword{to}
\keyword{tree}
