######################################################################
#
# access.R
#
# Written by Carter T. Butts <buttsc@uci.edu>; portions contributed by
# David Hunter <dhunter@stat.psu.edu> and Mark S. Handcock
# <handcock@u.washington.edu>.
#
# Last Modified 8/18/05
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/network package
#
# This file contains various routines for accessing network class objects.
#
# Contents:
#
#   add.edge
#   add.edges
#   add.vertices
#   delete.edge.attribute
#   delete.edges
#   delete.network.attribute
#   delete.vertex.attribute
#   delete.vertices
#   get.edge.attribute
#   get.edge.value
#   get.edgeIDs
#   get.edges
#   get.network.attribute
#   get.neighborhood
#   get.vertex.attribute
#   has.loops
#   is.adjacent
#   is.directed
#   is.hyper
#   is.multiplex
#   is.network
#   list.edge.attributes
#   list.network.attributes
#   list.vertex.attributes
#   network.edgecount
#   network.size
#   network.vertex.names
#   permute.vertexIDs
#   set.edge.attribute
#   set.edge.value
#   set.network.attribute
#   set.vertex.attribute
#
######################################################################


#Add a single edge to a network object.
#
add.edge<-function(x, tail, head, names.eval=NULL, vals.eval=NULL, edge.check=FALSE, ...){
  #Define a useful list equality function, for edge checks
  leq<-function(a,b){
    if(length(a)!=length(b))
      FALSE
    else
      all(a==b)
  }
  #Define some constants, for convenience
  n<-network.size(x)
  #Verify the integrity of the edge information (at some 
  #performance cost)
  if(edge.check){
    if(length(tail)*length(head)==0)	#Both sets must be nonempty
      stop("Empty head/tail list in add.edge.  Exiting.\n")
    if(any(head<1,tail<1,head>n,tail>n))  #Only allow legal vertices 
      stop("Illegal vertex reference in add.edge.  Exiting.\n")
    if(!is.hyper(x))			#Hyperedge check
      if((max(length(tail),length(head))>1))
        stop("Attempted to add hyperedge where hyper==FALSE in add.edge.  Exiting.")
    if((!is.multiplex(x))&&(length(x$mel)>0))
      if(is.directed(x)){
        if(any(sapply(lapply(x$mel,"[[","outl"),leq,tail)&
	sapply(lapply(x$mel,"[[","inl"),leq,head)))
          stop("Attempted to add multiplex edge where multiple==FALSE in add.edge.  Exiting.")
      }else{
        if(any((sapply(lapply(x$mel,"[[","outl"),leq,tail)&
	sapply(lapply(x$mel,"[[","inl"),leq,head))|
	(sapply(lapply(x$mel,"[[","outl"),leq,head)&
	sapply(lapply(x$mel,"[[","inl"),leq,tail))))
          stop("Attempted to add multiplex edge where multiple==FALSE in add.edge.  Exiting.")
      }
  }
  #Set up the edge attribute list
  atl<-as.list(vals.eval)
  if(!is.null(vals.eval)){
    if(length(names.eval)>length(vals.eval)){              #Too many labels!
      warning(paste("Too many labels in add.edge: wanted ",
      length(vals.eval)," got ",length(names.eval),". Truncating name list.\n",collapse=""))
      names(atl)<-names.eval[1:length(atl)]
    }else if(length(names.eval)<length(vals.eval)){	    #Too few labels!
      warning(paste("Too few labels in add.edge: wanted ",
      length(vals.eval)," got ",length(names.eval),".  Naming numerically.\n",collapse=""))
      names(atl)<-c(names.eval,
      (length(names.eval)+1):length(vals.eval))
    }else{                                           #Just the right number.
      names(atl)<-names.eval
    }
  }
  if(!("na"%in%names(atl)))		#Add the na attrib, if not present
    atl$na<-FALSE
  #Add the edge to the master edge list
  x$mel[[x$gal$mnext]]<-list(inl=as.vector(head), outl=as.vector(tail),atl=atl)
  #Add the edge pointer to the incoming/outgoing edge lists
  tf<-tail[[1]]
  hf<-head[[1]]
  for(i in tail){
    if(length(x$oel[[i]])>0){
      ihf<-sapply(lapply(x$mel[x$oel[[i]]],"[[","inl"),"[",1)
      x$oel[[i]]<-c(x$oel[[i]][ihf<hf],x$gal$mnext, x$oel[[i]][ihf>=hf])
    }else
      x$oel[[i]]<-x$gal$mnext
  }
  for(i in head){
    if(length(x$iel[[i]])>0){
      itf<-sapply(lapply(x$mel[x$iel[[i]]],"[[","outl"),"[",1)
      x$iel[[i]]<-c(x$iel[[i]][itf<tf],x$gal$mnext, x$iel[[i]][itf>=tf])
    }else
      x$iel[[i]]<-x$gal$mnext
  }
  #Increment the edge counter
  x$gal$mnext<-x$gal$mnext+1
  #Return the network
  x
}


# Add multiple edges to network x.  Tail must be a list, each element of
# which is the tail set for a given edge (ditto for head).  If edge values
# are provided, they must be given similarly as lists of lists.
#
add.edges<-function(x, tail, head, names.eval=NULL, vals.eval=NULL, ...){
  for(i in 1:length(tail))
    x<-add.edge(x,tail[[i]],head[[i]],names.eval[[i]],vals.eval[[i]],...)
  #Return the updated network
  x
}


# Add nv vertices to network x.  Vertex attributes (in addition to those which
# are required) are to be provided in vattr; vattr must be a list containing
# nv elements, each of which is equal to the desired val[i] entry.
#
add.vertices<-function(x, nv, vattr=NULL){
  #Set the size attribute
  x<-set.network.attribute(x,"n",network.size(x)+nv)
  #Add entries to the outgoing/incoming edge lists
  x$oel<-c(x$oel,replicate(nv,vector(mode="numeric")))
  x$iel<-c(x$iel,replicate(nv,vector(mode="numeric")))
  #Set up the vertex attributes
  if(is.null(vattr))
    vattr<-replicate(nv,list())
  nona<-sapply(lapply(vattr,"[[","na"),is.null)
  if(any(nona))
    vattr[nona]<-lapply(vattr[nona],function(a){a$na<-FALSE})
  #Add vertex attributes to the vertex attribute list
  x$val<-c(x$val,vattr)
  #Return the updated network
  x
}


# Remove all instances of the specified attribute from the edge set
#
delete.edge.attribute<-function(x,attrname){
  #Set the attribute to NULL, thereby deleting it.
  set.edge.attribute(x,attrname,NULL)
}
 
 
# Remove specified edges from the network.
#
delete.edges<-function(x,eid){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("delete.edges requires an argument of class network.")
  #Step through the edges, one by one
  for(i in eid){
    #Identify those in the head and tail
    head<-x$mel[[i]]$inl
    tail<-x$mel[[i]]$outl
    #Remove the edge from the respective lists
    for(j in head)
      x$iel[[j]]<-x$iel[[j]][x$iel[[j]]!=i]
    for(j in tail)
      x$oel[[j]]<-x$oel[[j]][x$oel[[j]]!=i]
    #Remove the edge itself
    x$mel[i]<-list(NULL)
  }
  #Return the network
  x
}


# Remove the specified network-level attribute
#
delete.network.attribute<-function(x,attrname){
  #Set the attribute to NULL, thereby deleting it.
  set.network.attribute(x,attrname,NULL)
}


# Remove all instances of the specified attribute from the vertex set
#
delete.vertex.attribute<-function(x,attrname){
  #Set the attribute to NULL, thereby deleting it.
  set.vertex.attribute(x,attrname,NULL)
}


# Remove specified vertices (and associated edges) from the network.
#
delete.vertices<-function(x,vid){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("delete.vertices requires an argument of class network.")
  n<-network.size(x)
  #Perform some sanity checks
  if(any((vid>n)|(vid<1)))
    stop("Vertex ID does not correspond to actual vertex in delete.vertices.")
  vid<-unique(vid)  #Get rid of duplicates, if any
  #Collect the list of edges having any vid vertex as an endpoint
  eid<-vector()
  for(i in vid){
    eid<-c(eid,get.edgeIDs(x,i,neighborhood="combined",na.omit=FALSE))
  }
  eid<-unique(eid)
  #Remove the edges
  x<-delete.edges(x,eid)
  #Permute the vid vertices to the end of the graph
  nord<-c((1:n)[-vid],vid)
  x<-permute.vertexIDs(x,nord)
  newn<-n-length(vid)
  #Now, get rid of the vertices
  x$val<-x$val[1:newn]
  x$iel<-x$iel[1:newn]
  x$oel<-x$oel[1:newn]
  x<-set.network.attribute(x,"n",newn)
  #Return the result
  x
}


# Retrieve a specified edge attribute from edge list el.  The attribute
# is returned as a list, regardless of type.
#
get.edge.attribute<-function(el, attrname, unlist=TRUE){
  x <- lapply(lapply(el,"[[","atl"),"[[",attrname)
  if(unlist){unlist(x)}else{x}
}


# Retrieve a specified edge attribute from all edges in x.  The attribute
# is returned as a list, regardless of type.
#
get.edge.value<-function(x, attrname, unlist=TRUE){
  y <- lapply(lapply(x$mel,"[[","atl"),"[[",attrname)
  if(unlist){unlist(y)}else{y}
}


# Retrieve the ID numbers for all edges incident on v, in network x.  
# Outgoing or incoming edges are specified by neighborhood, while na.omit 
# indicates whether or not missing edges should be omitted.  The return value
# is a vector of edge IDs.
#
get.edgeIDs<-function(x, v, alter=NULL, neighborhood=c("out","in","combined"), na.omit=TRUE){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("get.edgeIDs requires an argument of class network.")
  #Extract the edge IDs associated with the given neighborhood
  eid<-switch(match.arg(neighborhood),
    out=x$oel[[v]],
    "in"=x$iel[[v]],
    combined=c(x$oel[[v]],x$iel[[v]])
  )
  eid<-unique(eid)  #Remove duplicates
  el<-x$mel[eid]  #Get the edges themselves
  #If a specified alter is given, remove all edges not containing alter
  if(!is.null(alter))
    elalt<-switch(match.arg(neighborhood),
      out=sapply(el,function(a){alter%in%a$inl}),
      "in"=sapply(el,function(a){alter%in%a$outl}),
      combined=sapply(el,function(a){(alter%in%a$inl)||(alter%in%a$outl)})
    )
  else
    elalt<-rep(TRUE,length(el)) 
  eid<-eid[elalt]
  el<-el[elalt]
  #If needed, remove missing edges
  if((length(el)>0)&&(na.omit)){
    eid<-eid[!get.edge.attribute(el,"na",unlist=TRUE)]
  }
  #Return the result
  eid
}


# Retrieve all edges incident on v, in network x.  Outgoing or incoming
# edges are specified by neighborhood, while na.omit indicates whether
# or not missing edges should be omitted.  The return value is a list of
# edges.
#
get.edges<-function(x, v, alter=NULL, neighborhood=c("out","in","combined"), na.omit=TRUE){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("get.edges requires an argument of class network.")
  #Extract the edges associated with the given neighborhood
  el<-switch(match.arg(neighborhood),
    out=x$mel[x$oel[[v]]],
    "in"=x$mel[x$iel[[v]]],
    combined=x$mel[[unique(c(x$oel[[v]],x$iel[[v]]))]]
  )
  #If a specified alter is given, remove all edges not containing alter
  if(!is.null(alter))
    el<-switch(match.arg(neighborhood),
      out=el[sapply(el,function(a){alter%in%a$inl})],
      "in"=el[sapply(el,function(a){alter%in%a$outl})],
      combined=el[sapply(el,function(a){(alter%in%a$inl)||(alter%in%a$outl)})]
    )
  #If needed, remove missing edges
  if((length(el)>0)&&(na.omit))
    el<-el[!get.edge.attribute(el,"na",unlist=TRUE)]
  #Return the result
  el
}


# Retrieve a specified network-level attribute from network x.  The attribute
# type depends on the underlying storage mode, and cannot be guaranteed.
#
get.network.attribute<-function(x,attrname,unlist=FALSE){
  x <- x$gal[[attrname]]
  if(unlist){unlist(x)}else{x}
}


# Retrieve the neighborhood of v in network x.  Depending on the value of 
# type, the neighborhood in question may be in, out, or the union of the two.
# The return value for the function is a vector containing vertex IDs.
#
get.neighborhood<-function(x,v,type=c("out","in","combined")){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("get.neighborhood requires an argument of class network.")
  #Get requested edges
  if((!is.directed(x))||(match.arg(type)=="combined"))
    el<-c(get.edges(x,v,neighborhood="out",na.omit=TRUE),get.edges(x,v,neighborhood="in",na.omit=TRUE))
  else
    el<-switch(match.arg(type),
      out=get.edges(x,v,neighborhood="out",na.omit=TRUE),
      "in"=get.edges(x,v,neighborhood="in",na.omit=TRUE)
    )
  #Extract the neighborhoods
  raw<-switch(match.arg(type),
    out=unlist(lapply(el,"[[","inl")),
    "in"=unlist(lapply(el,"[[","outl")),
    combined=c(unlist(lapply(el,"[[","inl")),unlist(lapply(el,"[[","outl")))
  )
  #Return the unique members
  neigh<-unique(raw)
  neigh[neigh!=v]
}


# Retrieve a specified vertex attribute (indicated by attrname) from network x.
# Where na.omit==TRUE, values for missing vertices are removed; where
# null.na==TRUE, NULL values are converted to NAs.  The return value of this
# function is a list.
# 
get.vertex.attribute<-function(x,attrname,na.omit=FALSE,null.na=TRUE,
                               unlist=TRUE){
  #Get the list of attribute values
  va<-lapply(x$val,"[[",attrname)
  #If needed, figure out who's missing
  if(na.omit)
    vna<-lapply(x$val,"[[","na")
  else
    vna<-rep(FALSE,length(va))
  #Replace NULL values with NAs, if requested
  if(null.na)
    va[sapply(va,is.null)]<-NA
  #Return the result
  x <- va[!vna]
  if(unlist){unlist(x)}else{x}
}


# Return TRUE iff network x has loops.
#
has.loops<-function(x){
  if(!is.network(x))
    stop("has.loops requires an argument of class network.")
  else
    get.network.attribute(x,"loops")
}


# Return TRUE iff (vi,vj) in network x.  Where na.omit==TRUE, edges flagged
# as missing are ignored.
#
is.adjacent<-function(x,vi,vj,na.omit=TRUE){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("is.adjacent requires an argument of class network.")
  #Get the outedges associated with vi
  ol<-get.edges(x,vi,neighborhood="out",na.omit=FALSE)
  #Determine which edges are missing
  olna<-unlist(get.edge.attribute(ol,"na"))
  #See if any edges match
  olmatch<-sapply(lapply(ol,"[[","inl"),function(y){vj%in%y})
  #Return now, if we can
  if(length(ol)>0){
    if(any(olmatch*(!olna)))                #Actual match
      return(TRUE)
    else if(is.directed(x)){
      if(any(olmatch*olna)&&(!na.omit))     #Matches only on NA
        return(NA)
      else                                  #No valid matches
        return(FALSE)
    }
  }else
    if(is.directed(x))                      #Degenerate case
      return(FALSE)  
  #If we're still here, this must be an undirected graph.  Let's
  #look at vi's inedges
  il<-get.edges(x,vi,neighborhood="in",na.omit=FALSE)
  #Determine which edges are missing
  ilna<-unlist(get.edge.attribute(il,"na"))
  #See if any edges match
  ilmatch<-sapply(lapply(il,"[[","outl"),function(y){vj%in%y})
  #Return as needed
  if(length(il)==0)                       #Degenerate case
    return(FALSE)
  if(any(ilmatch*(!ilna)))                #Actual match
    return(TRUE)
  else if((!na.omit)&&(any(olmatch*olna,ilmatch*ilna)))  #Only NA
    return(NA)
  else
    return(FALSE)
}


# Return TRUE iff network x is directed.
#
is.directed<-function(x){
  if(!is.network(x))
    stop("is.directed requires an argument of class network.")
  else
    get.network.attribute(x,"directed")
}


# Return TRUE iff network x is hypergraphic.
#
is.hyper<-function(x){
  if(!is.network(x))
    stop("is.hyper requires an argument of class network.")
  else
    get.network.attribute(x,"hyper")
}


# Return TRUE iff network x is multiplex.
#
is.multiplex<-function(x){
  if(!is.network(x))
    stop("is.multiplex requires an argument of class network.")
  else
    get.network.attribute(x,"multiple")
}


# Return TRUE iff x is a network.
#
is.network<-function(x){
  inherits(x, "network")
}


# List attributes present on any edge
#
list.edge.attributes<-function(x){
  #First, check to see that this is a graph object
  if(!is.network(x))
    stop("list.network.attributes requires an argument of class network.")
  #Accumulate names
  allnam<-sapply(lapply(x$mel,"[[","atl"),names)
  #Return the sorted, unique attribute names
  sort(unique(as.vector(allnam)))
}


# List network-level attributes
#
list.network.attributes<-function(x){
  #First, check to see that this is a graph object
  if(!is.network(x))
    stop("list.network.attributes requires an argument of class network.")
  #Return the attribute names
  sort(names(x$gal))
}


# List attributes present on any vertex
#
list.vertex.attributes<-function(x){
  #First, check to see that this is a graph object
  if(!is.network(x))
    stop("list.network.attributes requires an argument of class network.")
  #Accumulate names
  allnam<-sapply(x$val,names)
  #Return the sorted, unique attribute names
  sort(unique(as.vector(allnam)))
}


#Retrieve the number of edges in network x.
#
network.edgecount<-function(x,na.omit=TRUE){
  #First, check to see if we can get out immediately
  if(get.network.attribute(x,"mnext")==1)
    return(0)
  #Determine null edges, and possibly missing edges as well
  enull<-sapply(x$mel,is.null)
  if(na.omit){
    ena<-get.edge.attribute(x$mel,"na",unlist=FALSE)
    ena[sapply(ena,is.null)]<-TRUE
    ena<-as.logical(unlist(ena))
  }else
    ena<-rep(FALSE,length(x$mel))
  #Return the edgecount
  sum((!enull)*(!ena))
}


# Retrieve the size (i.e., number of vertices) of network x.
#
network.size<-function(x){
  if(!is.network(x))
    stop("network.size requires an argument of class network.")
  else
    get.network.attribute(x,"n")
}


# Retrieve the vertex names of network x (if present).
#
network.vertex.names<-function(x){
  if(!is.network(x)){
    stop("network.vertex.names requires an argument of class network.")
  }else{
    vnames <- get.vertex.attribute(x,"vertex.names")
    if(is.null(vnames)  | all(is.na(vnames)) ){
      paste(1:network.size(x))
    }else{
      vnames
    }
  }
}


# Permute the internal IDs (ordering) of the vertex set
permute.vertexIDs<-function(x,vids){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("permute.vertexIDs requires an argument of class network.")
  #Sanity check: is this a permutation vector?
  n<-network.size(x)
  if(length(unique(vids))!=n)
    stop("Invalid permutation vector in permute.vertexIDs.")
  #First, determine which vertices need to change
  oord<-1:n
  cvids<-vids[vids!=oord]
  cnpos<-oord[vids!=oord]
  #Now, determine which edges need to change  
  eid<-vector()
  for(i in cvids)
    eid<-c(eid,get.edgeIDs(x,i,neighborhood="combined"))
  eid<-unique(eid)
  #For each such edge, change vertex IDs as needed
  for(i in eid){
    tofix<-x$mel[[i]]$inl%in%cvids
    x$mel[[i]]$inl[tofix]<-cnpos[match(x$mel[[i]]$inl[tofix],cvids)]
    tofix<-x$mel[[i]]$outl%in%cvids
    x$mel[[i]]$outl[tofix]<-cnpos[match(x$mel[[i]]$outl[tofix],cvids)]
  }
  #Now, reorder the vertex properties
  x$val<-x$val[vids]
  x$iel<-x$iel[vids]
  x$oel<-x$oel[vids]
  #Return the result
  x
}


# Set an edge attribute for network x.
#
set.edge.attribute<-function(x,attrname,value,e=1:length(x$mel)){
  for(i in 1:length(e))
    x$mel[[e[i]]]$atl[[attrname]]<-value[[i]]
  x
}


# Set an edge value for network x.
#
set.edge.value<-function(x,attrname,value,e=1:length(x$mel)){
  xmat <- as.matrix.network(x,matrix.type="adjacency")==1
  if(is.directed(x)){
   value <- value[xmat]
  }else{
   value[row(xmat)<col(xmat)] <- NA
   value <- value[xmat]
   value <- value[!is.na(value)]  
  }
  names(value) <- NULL
  for(i in 1:length(e)){
    x$mel[[e[i]]]$atl[[attrname]]<-value[i]
  }
  x
}


# Set a network-level attribute for network x.
#
set.network.attribute<-function(x,attrname,value){
  x$gal[[attrname]]<-value
  x
}


# Set a vertex attribute for network x.
#
set.vertex.attribute<-function(x,attrname,value,v=1:network.size(x)){
  for(i in 1:length(v))
    x$val[[v[i]]][[attrname]]<-value[[i]]
  x
}

