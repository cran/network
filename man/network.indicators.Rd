\name{network.indicators}
\alias{network.indicators}
\alias{has.loops}
\alias{is.directed}
\alias{is.hyper}
\alias{is.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Indicator Functions for Network Properties }
\description{
  Various indicators for properties of \code{network} class objects.
}
\usage{
has.loops(x)
is.directed(x)
is.hyper(x)
is.multiplex(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class \code{network} }
}
\details{
  \code{has.loops} returns \code{TRUE} iff \code{x} is allowed to contain loops (or loop-like edges, in the hypergraphic case).
  
  \code{is.directed} returns \code{TRUE} iff the edges of \code{x} are to be interpreted as directed.

  \code{is.hyper} returns \code{TRUE} iff \code{x} is allowed to contain hypergraphic edges.

  \code{is.multiplex} returns \code{TRUE} iff \code{x} is allowed to contain multiplex edges.
}
\value{
  \code{TRUE} or \code{FALSE}
}
\references{ Butts, C.T.  2002.  ``Memory Structures for Relational Data in R: Classes and Interfaces''  Working Paper. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ ~~further notes~~ }
%
% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{network}}, \code{\link{get.network.attribute}}, \code{set.network.attribute} }
\examples{
g<-network.initialize(5)    #Initialize the network
is.directed(g)
has.loops(g)
}
\keyword{ classes }% at least one, from doc/KEYWORDS
\keyword{ graphs }% __ONLY ONE__ keyword per line
