#required packages
#=================
library(splancs)
library(spatial)
library(network)
library(clue)
library(calibrate)


 
###GENERIC FUNCTIONS for plotting, comparing and matching spatial graphs arising in biometrics###
 
# plotspatialgraph(g) -- plots g as a spatial graph
# plotreggraphs(g1,g2) -- plots 2 registered graphs showing the aligned edge
# rapidscore(g1, g2, tol) -- rapid (and thoughtless) count of matching nodes of two spatial graphs
# shiftrotate <- function(g, edge) -- shifts and rotates a spatial graph such that edge is the origin and x-axis


#################################################################################
############# START plotspatialgraph - plots g as a spatial graph ###############
plotspatialgraph <- function(g){

true_edges <- as.matrix(g, matrix.type="edgelist")
xvals <- get.vertex.attribute(g, "xcoord")
yvals <- get.vertex.attribute(g, "ycoord")
true_vertices <- cbind(xvals, -yvals)
#true_vertices <- cbind(xvals, yvals)
#true_vertices <- cbind(yvals, xvals)

pointmap(true_vertices, add=FALSE, pch=20, col="black") #), xaxt="n",yaxt="n")
textxy(xvals,-yvals,seq(1:dim(true_vertices)[1]) )
#textxy(xvals,yvals,seq(1:dim(true_vertices)[1]) )
 
num_edges <- length(true_edges[,1])
if(num_edges>0){
for (i in 1:num_edges){
   start <- true_edges[i,1]
   end <- true_edges[i,2]
   added_edge <- matrix(c(true_vertices[start,1],true_vertices[end,1],true_vertices[start,2],true_vertices[end,2]), 2, 2)
   lines(added_edge, lwd=0.5, col="black")
}
}#end of if

}
############# END plotspatialgraph - plots g as a spatial graph ###############

##################################################################################################
## START plotreggraphs
###################################################################################################
plotreggraphs<-function(g1,g2){

g1_vertices<-cbind(get.vertex.attribute(g1,"xcoord"), get.vertex.attribute(g1,"ycoord"))
g2_vertices<-cbind(get.vertex.attribute(g2,"xcoord"), get.vertex.attribute(g2,"ycoord"))

g1_edgelist<-as.matrix.network(g1,matrix.type="edgelist")
g2_edgelist<-as.matrix.network(g2,matrix.type="edgelist")

#First "trick" R into making the plot big enough to see both graphs
#-----------------------------------------------------------------------
all_vertices<-rbind(g1_vertices,g2_vertices)
pointmap(all_vertices, add=FALSE, pch=20, col="white")

#plot vertices of g1
#====================

pointmap(g1_vertices, add=TRUE, pch=20, col="blue")
#textxy(g1_vertices[,1],g1_vertices[,2],seq(1:dim(g1_vertices)[1]),dcol="blue" )
 
#plot edges of g1
#===================

for (i in 1:length(g1_edgelist[,1])){
 
   start <- g1_edgelist[i,1]
   end <- g1_edgelist[i,2]

   added_edge <- matrix(c(g1_vertices[start,1],g1_vertices[end,1],g1_vertices[start,2],g1_vertices[end,2]), 2, 2)
 
   lines(added_edge, lwd=2, col="blue", lty=2)
}


#plot vertices of g2
#====================

pointmap(g2_vertices, add=TRUE, pch=10, col="red")
#textxy(g2_vertices[,1],g2_vertices[,2],seq(1:dim(g2_vertices)[1]),dcol="red")
 
#plot edges of g2
#===================

for (i in 1:length(g2_edgelist[,1])){
 
   start <- g2_edgelist[i,1]
   end <- g2_edgelist[i,2]

   added_edge <- matrix(c(g2_vertices[start,1],g2_vertices[end,1],g2_vertices[start,2],g2_vertices[end,2]), 2, 2)
 
   lines(added_edge, lwd=0.1, col="red")
}



}#end of function




##################################################################################################


###########################################################################################################
############# START rapidscore - "rapid" **count** of "common" nodes for two graphs g1 & g2 ###############
rapidscore <- function(g1, g2, tol){

#tol = size of bounding box to determine a vertex match

xvals <- get.vertex.attribute(g1, "xcoord")
yvals <- get.vertex.attribute(g1, "ycoord")
vertices_enr <- cbind(xvals, yvals)
num_vertices_enr <- length(vertices_enr[,1])

xvals <- get.vertex.attribute(g2, "xcoord")
yvals <- get.vertex.attribute(g2, "ycoord")
vertices_que <- cbind(xvals, yvals)
num_vertices_que <- length(vertices_que[,1])

match_score <- 0
vertices_enr_taken <- rep(FALSE, times=num_vertices_enr)
vertices_que_taken <- rep(FALSE, times=num_vertices_que)

for (i in 1:num_vertices_enr){   
   for (j in 1:num_vertices_que){

   if (vertices_enr_taken[i]) next
   if (vertices_que_taken[j]) next

   diff <- sqrt((vertices_enr[i,1]-vertices_que[j,1])**2+(vertices_enr[i,2]-vertices_que[j,2])**2)

   if (diff<=tol) {
                   match_score <- match_score + 1
                   vertices_enr_taken[i]<-TRUE
                   vertices_que_taken[j]<-TRUE
                  }
}
}
match_score <- 1 - match_score/(sqrt(num_vertices_enr*num_vertices_que))
match_score
}
############# END rapidscore - "rapid" **count** of "common" nodes for two graphs g1 & g2 ###############
###########################################################################################################

#########START rapidscore_matchlist - This function also returns the list of g1 nodes that matched g2 nodes############
###
rapidscore_matchlist <- function(g1, g2, tol){
  
  #tol = size of bounding box to determine a vertex match
  
  xvals <- get.vertex.attribute(g1, "xcoord")
  yvals <- get.vertex.attribute(g1, "ycoord")
  vertices_enr <- cbind(xvals, yvals)
  num_vertices_enr <- length(vertices_enr[,1])
  
  xvals <- get.vertex.attribute(g2, "xcoord")
  yvals <- get.vertex.attribute(g2, "ycoord")
  vertices_que <- cbind(xvals, yvals)
  num_vertices_que <- length(vertices_que[,1])
  
  match_score <- 0
  matchflag<-rep(10000, times=num_vertices_enr) # indicates the index of the matching node from query graph or 10000 if no match
  vertices_enr_taken <- rep(FALSE, times=num_vertices_enr)
  vertices_que_taken <- rep(FALSE, times=num_vertices_que)
  
  for (i in 1:num_vertices_enr){   
    for (j in 1:num_vertices_que){
      
      if (vertices_enr_taken[i]) next
      if (vertices_que_taken[j]) next
      
      diff <- sqrt((vertices_enr[i,1]-vertices_que[j,1])**2+(vertices_enr[i,2]-vertices_que[j,2])**2)
      
      if (diff<=tol) {
        match_score <- match_score + 1
        matchflag[i]<-j
        vertices_enr_taken[i]<-TRUE
        vertices_que_taken[j]<-TRUE
      }
    }
  }
  match_score <- 1 - match_score/(sqrt(num_vertices_enr*num_vertices_que))
  
 matchlist<-cbind(which( matchflag != 10000, arr.ind=TRUE), matchflag[ which(matchflag!=10000)])   
  
list(match_score=match_score,matchlist=matchlist)
}



##############################################################################################################################################
############# START shiftrotate - returns a rotated & shifted version of a spatial graph, with respect to a particular edge where ############
#############                     the start of the edge is the origin and the +ve x-axis is in the direction of the edge.         ############
shiftrotate <- function(g, edge){

edges <- as.matrix.network(g, matrix.type="edgelist")

xvals <- get.vertex.attribute(g, "xcoord")
yvals <- get.vertex.attribute(g, "ycoord")
vertices <- cbind(xvals, yvals)
num_vertices <- length(vertices[,1])

#APPLY shift and rotation to g
#=============================

start <- edges[edge,1]
end <- edges[edge,2]

if (vertices[start,1]>vertices[end,1]) {
#reverse labels of "start" and "end", so that start always refers to the node having the smallest x-coordinate in the original co-ordinate system
start <- edges[edge,2]
end <- edges[edge,1]
}

if (vertices[start,1]==vertices[end,1]) theta <- 90
if (vertices[start,1]<vertices[end,1]) theta <- (180/pi)*atan((vertices[end,2]-vertices[start,2])/(vertices[end,1]-vertices[start,1]))
if (vertices[start,1]>vertices[end,1]) theta <- (180/pi)*atan((vertices[start,2]-vertices[end,2])/(vertices[start,1]-vertices[end,1]))

newvertices <- vertices
newervertices <- vertices

#shift origin to co-ordinates of start
newvertices[,1] <- newvertices[,1]-vertices[start,1]
newvertices[,2] <- newvertices[,2]-vertices[start,2]

#rotate all points such that the edge (start, end) becomes a secant coinciding with the positive x-axis.
theta <- -(pi/180)*theta  #switches to radians
for (i in 1:num_vertices){
newervertices[i,1] <- newvertices[i,1]*cos(theta)-newvertices[i,2]*sin(theta)
newervertices[i,2] <- newvertices[i,2]*cos(theta)+newvertices[i,1]*sin(theta)
}


# g_aligned <- as.network.matrix(edges, matrix.type="edgelist", directed=FALSE)
# gsize <- network.size(g_aligned)
# if (num_vertices-gsize>0) add.vertices(g_aligned,num_vertices-gsize)  #add isolated vertices if necessary
g_aligned<-g

set.vertex.attribute(g_aligned, "xcoord", newervertices[,1])
set.vertex.attribute(g_aligned, "ycoord", newervertices[,2])

#update length and slope of g_aligned
#------------------------------------

lengths<-rep(0,times=length(edges[,1]) )
slopes<-rep(0,times=length(edges[,1]) )

for (i in 1:length(edges[,1])){

   start <- edges[i,1]
   end <- edges[i,2]

   lengths[i] <- sqrt((newervertices[start,1]-newervertices[end,1])**2+(newervertices[start,2]-newervertices[end,2])**2)
   
   if(newervertices[start,1]==newervertices[end,1]) slopes[i]<- 90
   if (newervertices[start,1]<newervertices[end,1]) slopes[i] <- (180/pi)*atan((newervertices[end,2]-newervertices[start,2])/(newervertices[end,1]-newervertices[start,1]))
   if (newervertices[start,1]>newervertices[end,1]) slopes[i] <- (180/pi)*atan((newervertices[start,2]-newervertices[end,2])/(newervertices[start,1]-newervertices[end,1]))

   lengths[i]<-round(lengths[i],2)
   slopes[i]<-round(slopes[i],2)

}

set.edge.attribute(g_aligned,"edgelength",lengths)
set.edge.attribute(g_aligned,"edgeslope",slopes)


list(g_aligned=g_aligned, theta=theta)
}
############# END shiftrotate - returns a rotated & shifted version of a spatial graph, with respect to a particular edge where ############
#############                   the start of the edge is the origin and the +ve x-axis is in the direction of the edge.         ############
############################################################################################################################################

##################Function plots the spatial graph with thickness effect in edges###############
plotspatialgraph_t <- function(g){
  
  true_edges <- as.matrix(g, matrix.type="edgelist")
  xvals <- get.vertex.attribute(g, "xcoord")
  yvals <- get.vertex.attribute(g, "ycoord")
  #true_vertices <- cbind(xvals, yvals)
  true_vertices <- cbind(xvals, -yvals)
  th<-get.edge.value(g, "thickness")
  
  
  pointmap(true_vertices, add=FALSE, pch=19, col="black", xaxt="n",yaxt="n")
  #textxy(xvals,yvals,seq(1:dim(true_vertices)[1]) , cx=0.8)
  
  num_edges <- length(true_edges[,1])
  if(num_edges>0){
    for (i in 1:num_edges){
      start <- true_edges[i,1]
      end <- true_edges[i,2]
      added_edge <- matrix(c(true_vertices[start,1],true_vertices[end,1],true_vertices[start,2],true_vertices[end,2]), 2, 2)
      if(th[i]< 1){lines(added_edge, lwd=0.5, col="blue", lty=1)}
      if(th[i]>= 1 && th[i]< 2 ){lines(added_edge, lwd=0.5, col="blue", lty=1)}
      if(th[i]>=2 && th[i]< 3 ){lines(added_edge, lwd=1.5, col="blue", lty=1)}
      if(th[i]>=3 && th[i]< 4 ){lines(added_edge, lwd=1.5, col="blue", lty=1)}
      if(th[i]>=4 && th[i]< 5 ){lines(added_edge, lwd=2.5, col="blue", lty=1)}
      if(th[i]>=5 && th[i]< 6 ){lines(added_edge, lwd=2.5, col="blue", lty=1)}
      if(th[i]>=6 && th[i]< 7 ){lines(added_edge, lwd=4, col="blue", lty=1)}
      if(th[i]>=7 && th[i]< 8 ){lines(added_edge, lwd=4, col="blue", lty=1)}
      
    }
  }#end of if
  
  
}
############# END plotspatialgraph_t - plots g as a spatial graph ###############

##################################################################################################
## START plotreggraphs_t - plots registered graphs with thickness effect
###################################################################################################
plotreggraphs_t<-function(g1,g2){
  
  g1_vertices<-cbind(get.vertex.attribute(g1,"xcoord"), get.vertex.attribute(g1,"ycoord"))
  g2_vertices<-cbind(get.vertex.attribute(g2,"xcoord"), get.vertex.attribute(g2,"ycoord"))
  
  g1_edgelist<-as.matrix.network(g1,matrix.type="edgelist")
  g2_edgelist<-as.matrix.network(g2,matrix.type="edgelist")
  
  g1_th<-get.edge.value(g1, "thickness")
  g2_th<-get.edge.value(g2, "thickness")
  
  #First "trick" R into making the plot big enough to see both graphs
  #-----------------------------------------------------------------------
  all_vertices<-rbind(g1_vertices,g2_vertices)
  pointmap(all_vertices, add=FALSE, pch=20, col="white")
  
  #plot vertices of g1
  #====================
  
  pointmap(g1_vertices, add=TRUE, pch=20, col="blue")
  #textxy(g1_vertices[,1],g1_vertices[,2],seq(1:dim(g1_vertices)[1]),dcol="blue" )
  
  #plot edges of g1
  #===================
  
  for (i in 1:length(g1_edgelist[,1])){
    
    start <- g1_edgelist[i,1]
    end <- g1_edgelist[i,2]
    
    added_edge <- matrix(c(g1_vertices[start,1],g1_vertices[end,1],g1_vertices[start,2],g1_vertices[end,2]), 2, 2)
    
    if(g1_th[i]< 1){lines(added_edge, lwd=0.5, col="blue", lty=1)}
    if(g1_th[i]>= 1 && g1_th[i]< 2 ){lines(added_edge, lwd=1, col="blue", lty=1)}
    if(g1_th[i]>=2 && g1_th[i]< 3 ){lines(added_edge, lwd=1.5, col="blue", lty=1)}
    if(g1_th[i]>=3 && g1_th[i]< 4 ){lines(added_edge, lwd=2, col="blue", lty=1)}
    if(g1_th[i]>=4 && g1_th[i]< 5 ){lines(added_edge, lwd=2.5, col="blue", lty=1)}
    if(g1_th[i]>=5 && g1_th[i]< 6 ){lines(added_edge, lwd=3, col="blue", lty=1)}
    if(g1_th[i]>=6 && g1_th[i]< 7 ){lines(added_edge, lwd=3.5, col="blue", lty=1)}
    if(g1_th[i]>=7 && g1_th[i]< 8 ){lines(added_edge, lwd=4, col="blue", lty=1)}
    
  }
  
  
  #plot vertices of g2
  #====================
  
  pointmap(g2_vertices, add=TRUE, pch=10, col="red")
  #textxy(g2_vertices[,1],g2_vertices[,2],seq(1:dim(g2_vertices)[1]),dcol="red")
  
  #plot edges of g2
  #===================
  
  for (i in 1:length(g2_edgelist[,1])){
    
    start <- g2_edgelist[i,1]
    end <- g2_edgelist[i,2]
    
    added_edge <- matrix(c(g2_vertices[start,1],g2_vertices[end,1],g2_vertices[start,2],g2_vertices[end,2]), 2, 2)
    
    
    if(g2_th[i]< 1){lines(added_edge, lwd=0.5, col="red", lty=2)}
    if(g2_th[i]>= 1 && g2_th[i]< 2 ){lines(added_edge, lwd=1, col="red", lty=2)}
    if(g2_th[i]>=2 && g2_th[i]< 3 ){lines(added_edge, lwd=1.5, col="red", lty=2)}
    if(g2_th[i]>=3 && g2_th[i]< 4 ){lines(added_edge, lwd=2, col="red", lty=2)}
    if(g2_th[i]>=4 && g2_th[i]< 5 ){lines(added_edge, lwd=2.5, col="red", lty=2)}
    if(g2_th[i]>=5 && g2_th[i]< 6 ){lines(added_edge, lwd=3, col="red", lty=2)}
    if(g2_th[i]>=6 && g2_th[i]< 7 ){lines(added_edge, lwd=3.5, col="red", lty=2)}
    if(g2_th[i]>=7 && g2_th[i]< 8 ){lines(added_edge, lwd=4, col="red", lty=2)}
          
  }
  
  
  
}#end of function

#############################################################################################################################

############# START rapidmatchscore - "rapid" **count** of "common" nodes for two graphs g1 & g2, returns the match count normalised by GM of vertex count of g1 and g2 ###############
rapidmatchscore <- function(g1, g2, tol){
  
  #tol = size of bounding box to determine a vertex match
  
  xvals <- get.vertex.attribute(g1, "xcoord")
  yvals <- get.vertex.attribute(g1, "ycoord")
  vertices_enr <- cbind(xvals, yvals)
  num_vertices_enr <- length(vertices_enr[,1])
  
  xvals <- get.vertex.attribute(g2, "xcoord")
  yvals <- get.vertex.attribute(g2, "ycoord")
  vertices_que <- cbind(xvals, yvals)
  num_vertices_que <- length(vertices_que[,1])
  
  match_score <- 0
  vertices_enr_taken <- rep(FALSE, times=num_vertices_enr)
  vertices_que_taken <- rep(FALSE, times=num_vertices_que)
  
  for (i in 1:num_vertices_enr){   
    for (j in 1:num_vertices_que){
      
      if (vertices_enr_taken[i]) next
      if (vertices_que_taken[j]) next
      
      diff <- sqrt((vertices_enr[i,1]-vertices_que[j,1])**2+(vertices_enr[i,2]-vertices_que[j,2])**2)
      
      if (diff<=tol) {
        match_score <- match_score + 1
        vertices_enr_taken[i]<-TRUE
        vertices_que_taken[j]<-TRUE
      }
    }
  }
  match_score
 
}

########################################################################################################
# Function to detect degree 2 nodes, remove them and and connect their neighbors in palm vein graphs ## 
########################################################################################################
getSmoothGraph_Pv<-function(g){

  N<-network.size(g)
  g_new<-g  
 vertices_todelete<-numeric()
for(i in 1:N){  
  n_v<-get.neighborhood(g,i)
  if(length(n_v)!=2) next
  vertices_todelete<-c( vertices_todelete, i)
   add.edge(g_new, n_v[1], n_v[2])
}#end of i loop
  delete.vertices(g_new, vertices_todelete)
  g_new  
}#end of function
