######################################################################
#
# operators.R
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
# This file contains various operators which take networks as inputs.
#
# Contents:
#
# None, as of 8/18/05; %c% was removed due to naming conflict with sna
#
######################################################################


#Graph composition operator; at present, this works crudely, via
#matrix representations
#NOTE: use sna's version instead.
#"%c%"<-function(x,y){
#  #Convert to adjacency form
#  x<-as.matrix.network(x,matrix.type="adjacency")
#  y<-as.matrix.network(y,matrix.type="adjacency")
#  #Check for conformability
#  if(dim(x)[2]!=dim(y)[1])
#    stop("Non-conformable relations in %c%.  Cannot compose.")
#  #Obtain the composed graph
#  network(round((x%*%y)>0),loops=TRUE)
#}

