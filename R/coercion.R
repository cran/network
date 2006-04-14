######################################################################
#
# coercion.R
#
# Written by Carter T. Butts <buttsc@uci.edu>; portions contributed by
# David Hunter <dhunter@stat.psu.edu> and Mark S. Handcock
# <handcock@u.washington.edu>.
#
# Last Modified 4/10/06
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/network package
#
# This file contains various routines for coercion to/from network
# class objects.
#
# Contents:
#
#   as.matrix.network
#   as.matrix.network.adjacency
#   as.matrix.network.edgelist
#   as.matrix.network.incidence
#   as.network
#   as.network.default
#   as.network.network
#   as.network.matrix
#   as.sociomatrix
#
######################################################################


# Method for general coercion of network class objects into matrices.
# Matrix type is indicated by the eponymous argument; note that some
# types may not be supported for certain networks.  Where
# attrname!=NULL, an edge attribute of name attrname is used to supply
# edge values.  Otherwise, edges are assumed to be unvalued.
#
as.matrix.network<-function(x,matrix.type=NULL,attrname=NULL){
  #Get the matrix type
  if(is.null(matrix.type))
    matrix.type<-which.matrix.type(x)
  else
    matrix.type<-match.arg(matrix.type,c("adjacency","incidence","edgelist"))
  #Dispatch as needed
  switch(matrix.type,
    adjacency=as.matrix.network.adjacency(x=x,attrname=attrname),
    incidence=as.matrix.network.incidence(x=x,attrname=attrname),
    edgelist=as.matrix.network.edgelist(x=x,attrname=attrname)
  )
}


# Coerce a network object to an adjacency matrix (where possible).  If
# provided, attrname is used to identify an attribute to use for edge
# values.
#
as.matrix.network.adjacency<-function(x,attrname=NULL){
  #Check to make sure this is a supported network type
  if(is.hyper(x))
    stop("Hypergraphs not currently supported in as.matrix.network.adjacency.  Exiting.\n")
  if(is.multiplex(x))
    stop("Multigraphs not currently supported in as.matrix.network.adjacency.  Exiting.\n")
  #Generate the adjacency matrix 
  m<-matrix(0,nr=network.size(x),nc=network.size(x))
  tl<-unlist(sapply(x$mel,"[[","outl")) #Can unlist b/c no hyperedges
  hl<-unlist(sapply(x$mel,"[[","inl"))
  nal<-as.logical(get.edge.attribute(x$mel,"na",unlist=TRUE))
  if(!is.null(attrname)){
    val<-unlist(get.edge.attribute(x$mel,attrname))
    if(is.null(val)){
     warning(paste("There is no edge attribute named", attrname))
     val<-rep(1,length(tl))
    }
  }else{
    val<-rep(1,length(tl))
  }
  if(length(hl[!nal])>0){
    m[tl[!nal]+(hl[!nal]-1)*network.size(x)]<-val[!nal]
  }
  if(length(hl[ nal])>0){
   m[tl[ nal]+(hl[ nal]-1)*network.size(x)]<-NA
  }
  #If undirected, symmetrize
  if(!is.directed(x)){
# changed by MSH to allow non binary values
#   m<-pmax(m,t(m))
    m[m==0] <- t(m)[m==0]
  }
  #Set row/colnames to vertex names
  xnames <- network.vertex.names(x)
  dimnames(m) <- list(xnames, xnames)
  #If bipartite extract only it
  if(is.bipartite(x)){
    nactors <- get.network.attribute(x, "bipartite")
    nevents <- network.size(x) - nactors
    m <- m[1:nactors, nactors+(1:nevents)]
  }
  #Return the result
  m
}


# Coerce a network object to an edgelist matrix.  If provided, attrname is 
# used to identify an attribute to use for edge values.
#
as.matrix.network.edgelist<-function(x,attrname=NULL){
  #Check to make sure this is a supported network type
  if(is.hyper(x))
    stop("Hypergraphs not currently supported in as.matrix.network.edgelist.  Exiting.\n")
  #Find the missing edges
  nal<-as.logical(get.edge.attribute(x$mel,"na"))
  #Generate the edgelist matrix
  m<-cbind(unlist(sapply(x$mel,"[[","outl")), unlist(sapply(x$mel,"[[","inl")))
  #Add edge values, if needed
  if(!is.null(attrname))
    m<-cbind(m,unlist(get.edge.attribute(x$mel,attrname)))
  #Return the result
  m[!nal,]
}


# Coerce a network object to an incidence matrix (where possible).  If
# provided, attrname is used to identify an attribute to use for edge
# values.
#
as.matrix.network.incidence<-function(x,attrname=NULL){
  #Perform preprocessing
  n<-network.size(x)
  inl<-lapply(x$mel,"[[","inl")
  outl<-lapply(x$mel,"[[","outl")
  if(!is.null(attrname))
    evals<-unlist(get.edge.attribute(x$mel,attrname))
  else
    evals<-rep(1,length(x$mel))
  ena<-as.logical(get.edge.attribute(x$mel,"na"))
  #Generate the incidence matrix
  dir<-is.directed(x)
  f<-function(a,m,k){y<-rep(0,m); y[a]<-k; y}
  im<-sapply(inl,f,n,1)+sapply(outl,f,n,ifelse(dir,-1,1))
  if(!dir)
    im<-pmin(im,1)
  im<-sweep(im,2,evals,"*")              #Fill in edge values
  im[sapply(ena,rep,n)*(im!=0)]<-NA      #Add NAs, if needed
  #Return the result
  im
}


as.network<-function(x,...)
  UseMethod("as.network")


as.network.default<-function(x,...)
  as.network.matrix(x,...)


as.network.network<-function(x,...)
  x


#
# MSH modified for bipartite
#
as.network.matrix<-function(x, matrix.type=NULL,
        directed=TRUE, hyper=FALSE, loops=FALSE, multiple=FALSE,
        bipartite=FALSE,
        ignore.eval=TRUE, names.eval=NULL, na.rm=FALSE, edge.check=FALSE, ...){
  if(is.logical(x)){x <- 1*x}
  #Get the matrix type
  if(is.null(matrix.type))
    matrix.type<-which.matrix.type(x)
  else
    matrix.type<-match.arg(matrix.type,c("adjacency","incidence","edgelist",
                                         "bipartite"))
  # Add names if available
  unames <- NULL
  if(matrix.type=="edgelist"){
    if(dim(x)[2]>2)
      vals<-x[,-(1:2)]
    else
      vals<-NULL
    if(is.character(x<-as.matrix(x[,1:2,drop=FALSE]))){
      unames <- sort(unique(as.vector(x)))
      x <- cbind(match(x[,1],unames),match(x[,2],unames))
    }
    if(!is.null(vals))
      x<-cbind(x,vals)
  }
  if(matrix.type=="adjacency" && !is.null(colnames(x))){
    unames <- colnames(x)
  }
  if(matrix.type=="bipartite"){
   directed <- FALSE
   bipartite <- dim(x)[1]
   unames <- 1:sum(dim(x))
   if(!is.null(rownames(x))){
     unames[1:(dim(x)[1])] <- rownames(x)
   }
   if(!is.null(colnames(x))){
     unames[(dim(x)[1])+(1:(dim(x)[2]))] <- colnames(x)
   }
  }
  #Initialize the network object
  n<-switch(matrix.type,	#Extract n based on matrix type
    adjacency=dim(x)[1],
    incidence=dim(x)[1],
    bipartite=sum(dim(x)),
    edgelist=max(x[,1:2]),
  )
  g<-network.initialize(n,directed=directed, hyper=hyper, loops=loops, multiple=multiple,bipartite=bipartite)
  #Call the specific coercion routine, depending on matrix type
  g<-switch(matrix.type,
    adjacency=network.adjacency(x,g,
     ignore.eval,names.eval,na.rm,edge.check),
    incidence=network.incidence(x,g,
     ignore.eval,names.eval,na.rm,edge.check),
    bipartite=network.bipartite(x,g,
     ignore.eval,names.eval,na.rm,edge.check),
    edgelist=network.edgelist(x,g, 
     ignore.eval,names.eval,na.rm,edge.check)
  )

  if(!is.null(unames)){
   g <- set.vertex.attribute(g,"vertex.names", unames)
  }
  #Return the result
  g
}


#Force the input into sociomatrix form.  This is a shortcut to 
#as.matrix.network.adjacency, which ensures that a raw matrix is
#passed through as-is.
as.sociomatrix<-function(x, attrname=NULL, simplify=TRUE){
  if(is.network(x)){ #If network, coerce to adjacency matrix
    g<-as.matrix.network.adjacency(x,attrname=attrname)
  }else if(is.matrix(x)||is.array(x)){ #If an array/matrix, use as-is
    g<-x
  }else if(is.list(x)){  #If a list, recurse on list elements
    g<-lapply(x,as.sociomatrix,attrname=attrname,simplify=simplify)
  }else{
    stop("as.sociomatrix input must be an adjacency matrix/array, network, or list.")
  }
  #Convert into the appropriate return format
  if(is.list(g)){   #Collapse if needed
    if(length(g)==1){
      g<-g[[1]]
      if((!simplify)&&(length(dim(g))==3)){  #Coerce to a list of matrices?
        out<-list()
        for(i in 1:dim(g)[1])
          out[[i]]<-g[i,,]
      }else{
        out<-g
      }
    }else{
      #Coerce to array form?
      if(simplify){
        dims<-sapply(g,dim)
        if(is.list(dims)){      #Dims must not be of equal length
          mats<-sapply(dims,length)
          mats[mats==1]<-0
          mats[mats==2]<-1
          mats[mats==3]<-sapply(dims[mats==3],"[[",1)
          mats<-cumsum(mats)
          dims<-sapply(dims,"[",2)
        }else{                  #Dims are of equal length
          if(NROW(dims)==3)      #Determine number of matrices per entry
            mats<-cumsum(dims[1,])
          else
            mats<-1:NCOL(dims)
          dims<-dims[2,]         #Get ncols
        }
        if((!any(is.null(dims)))&&(length(unique(dims))==1)&&(all(mats>0))){
          out<-array(dim=c(mats[length(mats)],dims[1],dims[1]))
          for(i in 1:length(mats))
            out[(c(0,mats)[i]+1):(mats[i]),,]<-g[[i]]
        }else
          out<-g
      }else
        out<-g
    }
  }else{
    if((!simplify)&&(length(dim(g))==3)){  #Coerce to a list of matrices?
      out<-list()
      for(i in 1:dim(g)[1])
        out[[i]]<-g[i,,]
    }else
      out<-g
  }
  #Return the result
  out
}