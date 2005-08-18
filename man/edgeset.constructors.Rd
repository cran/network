\name{edgeset.constructors}
\alias{edgeset.constructors}
\alias{network.adjacency}
\alias{network.incidence}
\alias{network.edgelist}
\alias{network.bipartite}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Edgeset Constructors for Network Objects }
\description{
  These functions convert relational data in matrix form to network edge sets.
}
\usage{
network.adjacency(x, g, ignore.eval = TRUE, names.eval = NULL, ...)
network.edgelist(x, g, ignore.eval = TRUE, names.eval = NULL, ...)
network.incidence(x, g, ignore.eval = TRUE, names.eval = NULL, ...)
network.bipartite(x, g, ignore.eval = TRUE, names.eval = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a matrix containing edge information }
  \item{g}{ an object of class \code{network} }
  \item{ignore.eval}{ logical; ignore edge values? }
  \item{names.eval}{ the edge attribute under which to store edge values, if any }
  \item{\dots}{ additional arguments to \code{\link{add.edge}} }
}
\details{
  Each of the above functions takes a \code{network} and a matrix as input, and returns a \code{network} object whose edgeset reflects the given information.  \code{network.adjacency} takes \code{x} to be an adjacency matrix; \code{code.edgelist} takes \code{x} to be an edgelist matrix; and \code{network.incidence} takes \code{x} to be an incidence matrix.  \code{network.bipartite} takes \code{x} to be an adjacency matrix where the rows are thge actors and the columns are the events; If \code{ignore.eval==FALSE}, (non-zero) edge values are stored as edgewise attributes with name \code{names.eval}.  Any additional command line parameters are passed to \code{\link{add.edge}}.
}
\value{
 An object of class \code{network}
}
\references{ Butts, C.T.  2002.  ``Memory Structures for Relational Data in R: Classes and Interfaces''  Working Paper. }
\author{ Carter T. Butts \email{buttsc@uci.edu} and David Hunter \email{dhunter@stat.psu.edu} }
\note{ Handling of missing data is not yet fully implemented. }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{network}}, \code{\link{network.initialize}}, \code{\link{add.edge}} }
\examples{
#Create an arbitrary adjacency matrix
m<-matrix(rbinom(25,1,0.5),5,5)
diag(m)<-0

g<-network.initialize(5)    #Initialize the network
g<-network.adjacency(m,g)   #Import the edge data
}
\keyword{ classes }% at least one, from doc/KEYWORDS
\keyword{ graphs }% __ONLY ONE__ keyword per line
