\name{add.edges}
\alias{add.edges}
\alias{add.edge}
\title{Add Edges to a Network Object}
\description{
  Add one or more edges to an existing network object.
}
\usage{
add.edge(x, tail, head, names.eval=NULL, vals.eval=NULL, 
    edge.check=FALSE, ...)
add.edges(x, tail, head, names.eval=NULL, vals.eval=NULL, ...)
}
\arguments{
  \item{x}{an object of class \code{network}}
  \item{tail}{for \code{add.edge}, a vector of vertex IDs reflecting the tail set for the edge to be added; for \code{add.edges}, a list of such vectors}
  \item{head}{for \code{add.edge}, a vector of vertex IDs reflecting the head set for the edge to be added; for \code{add.edges}, a list of such vectors}
  \item{names.eval}{for \code{add.edge}, an optional list of names for edge attributes; for \code{add.edges}, a list of such lists}
  \item{vals.eval}{for \code{add.edge}, an optional list of edge attribute values (matching \code{names.eval}); for \code{add.edges}, a list of such lists}
  \item{edge.check}{ logical; should we perform (computationally expensive) tests to check for the legality of submitted edges? }
  \item{...}{ additional arguments }
}
\details{
The edge checking procedure is very slow, but should always be employed when debugging; without it, one cannot guarantee that the network state is consistent with network level variables (see \code{\link{network.indicators}}).
}
\value{
\code{add.edges} returns the network object with the edges added.
}
\author{ Carter T. Butts \email{buttsc@uci.edu}}
\references{ Butts, C.T.  2002.  ``Memory Structures for Relational Data in R: Classes and Interfaces''  Working Paper. }
\note{ As with many routines in the network class, the iterative nature of these procedures makes them very slow within R.  Eventually, such things will be handled via C backends. }
\seealso{\code{\link{network}}, \code{\link{add.vertices}}}
\examples{
#Initialize a small, empty network
g<-network.initialize(3)

#Add an edge
g<-add.edge(g,1,2)
g
}
\keyword{classes}
\keyword{graphs}
