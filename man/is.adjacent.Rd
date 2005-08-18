\name{is.adjacent}
\alias{is.adjacent}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Determine Whether Two Vertices Are Adjacent }
\description{
  \code{is.adjacent} returns \code{TRUE} iff \code{vi} is adjacent to \code{vj} in \code{x}.  Missing edges may be omitted or not, as per \code{na.omit}.
}
\usage{
is.adjacent(x, vi, vj, na.omit = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class \code{network} }
  \item{vi}{ a vertex ID }
  \item{vj}{ a second vertex ID }
  \item{na.omit}{ logical; should missing edges count when assessing adjacency? }
}
\details{
  
}
\value{
  A logical, giving the status of the (i,j) edge
}
\references{ Wasserman, S. and Faust, K.  1994.  \emph{Social Network Analysis: Methods and Applications}.  Cambridge: Cambridge University Press. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ ~~further notes~~ }
%
% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{get.neighborhood}} }
\examples{
#Create a very simple graph
g<-network.initialize(3)
g<-add.edge(g,1,2)
is.adjacent(g,1,2)  #TRUE
is.adjacent(g,2,1)  #FALSE
}
\keyword{ graphs }% at least one, from doc/KEYWORDS
