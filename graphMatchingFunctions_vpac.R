#required packages
#=================
#setwd("I:/Academic/Research/Projects and Grants/Biometric Databases/Graph Databases")

#setwd("D:/Arathi/Graph Databases") #copied graphs to local machine to increase speed
library(splancs)
library(spatial)
library(network)
library(clue)
source("basicGraphFunctions_vpac.R")


#Functions defined in this file
#================================
# getGraph<-function(biometric,sourceid,person,sample)
# registerGraphs<-function(g1,g2,tol,f_plot)
# getGraphEditPath<-function(g1,g2,cid, f_edgeinclude)
# getEdgeEditCost<-function(g1, i, g2, j, cid)
# getMCS<-function(g1,g2, edits)
########################################################
# getDist<-function(g1,g2,mcs) -- function that computes several types of distances between 2 graphs and returns that list
# getGenImpScores<-function(biometric,sourceid,indexfilelist,tol,cid,f_edgeinclude,f_plot)
# getEER(genuines, imposters, genscorematrix, impscorematrix) -- gets an equal error rate from two sets of normalised distances
# getROC(genuines, imposters, genscorematrix, impscorematrix, f_plot) --returns a matrix of FNMR and FMR at a range of thresholds
# registerGraphs_maxtheta<-function(g1,g2,tol,maxtheta, f_plot) -- registers based on a max bound for rotation of frame

#============================================================================


##############################################################################
#This function creates a graph object with edge and vertex attributes.
# source->AutoJJ, AutoSL, Manual
# biometric->FP,BP,RE,HV
##############################################################################

getGraph<-function(biometric,sourceid,person,sample){

inputNodesFile<-""
inputEdgesFile<-""

# Identify the directory where graphs are stored from biometric and source
#-------------------------------------------------------------------------
if (biometric=="FP") {
inputNodesFile<-paste(getwd(),"/FP_",person,"_",sample,"_nodes.txt",sep="")
inputEdgesFile<-paste(getwd(),"/FP_",person,"_",sample,"_edges.txt", sep="")
} # end of if (biometric=="FP")

if(biometric=="RE" && sourceid=="Manual"){
inputNodesFile<-paste(getwd(),"/RE_",person,"_",sample,"_nodes.txt",sep="")
inputEdgesFile<-paste(getwd(),"/RE_",person,"_",sample,"_edges.txt", sep="")
} # end of if (biometric=="RE")

if(biometric=="RE" && sourceid=="AutoSL_TIP"){
  inputNodesFile<-paste(getwd(),"/R_00",person,"_",sample,"_bpoints.txt",sep="")
  inputEdgesFile<-paste(getwd(),"/R_00",person,"_",sample,"_links.txt", sep="")
} # end of if (biometric=="RE")

if(biometric=="RE" && sourceid=="AutoSLJK"){
  
  if(person<10) personid<-paste("0", person, sep="")
  if(person>=10 ) personid<-person
  
  inputNodesFile<-paste(wdname,"/RE_",personid,"_",sample,"_nodes.txt",sep="")
  inputEdgesFile<-paste(wdname,"/RE_",personid,"_",sample,"_edges.txt", sep="")
  
  
} # end of if (biometric=="RE")

if(biometric=="HV"){
  inputNodesFile<-paste(getwd(),"/p_",person,"_",sourceid[2],"_",sample,"_bpoints.txt",sep="")
  inputEdgesFile<-paste(getwd(),"/p_",person,"_",sourceid[2],"_",sample,"_links.txt",sep="")
  
}

#Create graph
#-------------
v <- read.table(inputNodesFile,header=FALSE, sep="\t", fill=TRUE)
N<-length(v[,1])
e<-read.table(inputEdgesFile,header=FALSE, sep="\t", fill=TRUE)
e<-as.matrix(e)
e<-e[,1:2] #for SNIR_T_LR mainly , lets ignore T for now.

#arrange edges with lower vertex index first and higher vertex index second #added on 24/02/2015
toswap<-which(e[,1]>e[,2])
for(k in toswap){
  val<-e[k,]
  newval<-c(val[2], val[1])
  e[k,]<-newval
}
#remove repeated edges
e<-unique(e)

g <- as.network.matrix(e, matrix.type="edgelist",directed=FALSE)
gsize <- network.size(g)
if (N-gsize>0) add.vertices(g,N-gsize)  #add isolated vertices if necessary
set.vertex.attribute(g, "xcoord", v[,2])
set.vertex.attribute(g, "ycoord", v[,3])
if(biometric=="HV"){
  set.vertex.attribute(g,"type", v[,4]) 
  
}


#delete edges where start and end nodes are same.
todelete<-which(e[,1]==e[,2])
delete.edges(g,todelete)


lengths<-rep(0,times=length(e[,1]) )
slopes<-rep(0,times=length(e[,1]) )

for (i in 1:length(e[,1])){

   start <- e[i,1]
   end <- e[i,2]

   lengths[i] <- sqrt((v[start,2]-v[end,2])**2+(v[start,3]-v[end,3])**2)
   
   if((v[end,2]-v[start,2])==0){slopes[i]=90}
   if((v[end,2]-v[start,2])!=0){
     if (v[start,2]<v[end,2]) slopes[i] <- (180/pi)*atan((v[end,3]-v[start,3])/(v[end,2]-v[start,2]))
     if (v[start,2]>=v[end,2]) slopes[i] <- (180/pi)*atan((v[start,3]-v[end,3])/(v[start,2]-v[end,2]))
         
   }
   
   lengths[i]<-round(lengths[i],2)
   slopes[i]<-round(slopes[i],2)
 #print(i)
}

set.edge.attribute(g,"edgelength",lengths)
set.edge.attribute(g,"edgeslope",slopes)

 #add ridgecount for BP




#return graph
#-------------
g

}#end of getGraph function
###########################################################################


################################################################################
## START registerGraphs function.
## This function takes 2 graphs as input, aligns them on the best aligning edge
## Also requires tol, the maximum matching tolerance for graph vertices
## It returns a list of 4 objects: the shifted and rotated graphs corresponding to 
## matching edge pair, the best matching edge pair and the best match score 
## (lower is better as score measures distance)
## It can show the best alignment in a plot, if the flag(f_plot) is set to 1.
####################################################################################

registerGraphs<-function(g1,g2,tol,f_plot){

g1_edgelength<-get.edge.value(g1,"edgelength")
g2_edgelength<-get.edge.value(g2,"edgelength")
g1_edgeslope<-get.edge.value(g1,"edgeslope")
g2_edgeslope<-get.edge.value(g2,"edgeslope")

g1_edgelist<-as.matrix.network(g1, matrix.type="edgelist")
g2_edgelist<-as.matrix.network(g2, matrix.type="edgelist")


#match edges by comparing lengths and slopes
#===========================================

edge_matching <- matrix(0, length(g1_edgelength), length(g2_edgelength) )


for (i in 1:length(edge_matching[,1])){
   for (j in 1:length(edge_matching[1,])){
     
#      d1<-length(get.neighborhood(g1,g1_edgelist[i,1], type="combined"))+length(get.neighborhood(g1,g1_edgelist[i,2], type="combined"))
#      d2<-length(get.neighborhood(g2,g2_edgelist[j,1], type="combined"))+length(get.neighborhood(g2,g2_edgelist[j,2], type="combined"))
#      

     #edge_matching[i,j] <- abs(g1_edgelength[i]-g2_edgelength[j])
     
   edge_matching[i,j] <- sqrt((g1_edgelength[i]-g2_edgelength[j])**2+(g1_edgeslope[i]-g2_edgeslope[j])**2)
   edge_matching[i,j] <- edge_matching[i,j]/(0.5*g1_edgelength[i]+0.5*g2_edgelength[j])
     #edge_matching[i,j] <- edge_matching[i,j]/(d1+d2)
  

}
}

ranking <- order(edge_matching, decreasing=FALSE)



#take first 50 closest edge pairs to realign and get best match score
#================================================================
NEdgePairs<-0.1*length(ranking)
#NEdgePairs<-length(ranking)
if(length(ranking)<NEdgePairs) { NEdgePairs<- length(ranking) }
bestMatchScore<-1
bestAlignEdge<-rep(0,2)
g1_shiftrot<-g1
g2_shiftrot<-g2
scores<-numeric()
alignedEdge<-numeric()
bestRank<-0
                 

for(i in 1:NEdgePairs) {
#for(i in 2:2){
  
rank<-i

column <- ceiling(ranking[rank]/length(edge_matching[,1]))
row <- ranking[rank]%%length(edge_matching[,1])
if (row==0) row <- length(edge_matching[,1])
edgepair <- c(row, column)
alignedEdge<-rbind(alignedEdge,edgepair)

#print(edgepair)

g1_temp<-shiftrotate(g1,row)$g_aligned #new, after shiftrot changed to return list
g2_temp<-shiftrotate(g2,column)$g_aligned
tempscore<-rapidscore(g1_temp,g2_temp,tol)

scores<-c(scores,tempscore)
   
             

if(tempscore<bestMatchScore){
bestMatchScore<-tempscore
bestAlignEdge<-c(row,column)
g1_shiftrot<-g1_temp
g2_shiftrot<-g2_temp
bestRank<-rank
}


} #end of i loop

if(f_plot==1){
windows()
plotreggraphs(g1_shiftrot,g2_shiftrot)
title("registered graphs")
}

output<-list(graph1=g1_shiftrot,graph2=g2_shiftrot,bestEdgePair=bestAlignEdge,bestScore=bestMatchScore,bestRank=bestRank)

}# end of function
########################################################################################################

###################################################################################
########### BEGIN getGraphEditPath ################################################
## This function takes 2 ALIGNED graphs and cid, the insertion/deletion cost as input
## It has a fourth input, a flag that selects if edges must be included in deciding the edit path (1 include, 0 exculde)
## It outputs the cheapest graph edit path in going from g1 to g2
## The cost matrix uses euclidean distance between vertex attributes as the basis 
## for computing cost


getGraphEditPath<-function(g1,g2,cid,f_edgeinclude){
g1_vertices<-cbind(get.vertex.attribute(g1,"xcoord"),get.vertex.attribute(g1,"ycoord"))
g2_vertices<-cbind(get.vertex.attribute(g2,"xcoord"),get.vertex.attribute(g2,"ycoord"))

N1<-network.size(g1)
N2<-network.size(g2)

##Build the cost matrix
############################

block1<-matrix(0,N1, N2)

for (i in 1:N1){
     for(j in 1:N2){
       
	sij<-sqrt( (g1_vertices[i,1]-g2_vertices[j,1])**2+(g1_vertices[i,2]-g2_vertices[j,2])**2 )    
#   if(sij<=5){sij=0}
#   if(sij>5){sij=100000}
  
	if(f_edgeinclude){
		edgeedit_cost<-getEdgeEditCost(g1,i,g2,j,cid)
		}
	else{
	edgeedit_cost<-0
		}

	block1[i,j]<-sij+edgeedit_cost
	}#end of j loop
}#end of i loop

block2<-matrix(100000,N1,N1)
for(i in 1:N1) 
	block2[i,i]<- cid+ (cid*length(get.neighborhood(g1,i)))

block3<-matrix(100000,N2,N2)
for(j in 1:N2)
	block3[j,j]<- cid + (cid*length(get.neighborhood(g2,j)))

block4<-matrix(0,N2,N1)

upper<-cbind(block1, block2)
lower<-cbind(block3, block4)
C<-rbind(upper,lower)

x<-solve_LSAP(C, maximum=FALSE)
edits <- cbind(seq_along(x), x)
GED <- sum(C[edits])

list(editPath=edits, editCost=GED)

}#end of graph edit path
###################################################################################

getGraphEditPath_twocids<-function(g1,g2,cid_n, cid_e,f_edgeinclude){
  g1_vertices<-cbind(get.vertex.attribute(g1,"xcoord"),get.vertex.attribute(g1,"ycoord"))
  g2_vertices<-cbind(get.vertex.attribute(g2,"xcoord"),get.vertex.attribute(g2,"ycoord"))
  
  N1<-network.size(g1)
  N2<-network.size(g2)
  
  ##Build the cost matrix
  ############################
  
  block1<-matrix(0,N1, N2)
  
  for (i in 1:N1){
    for(j in 1:N2){
      
      sij<-sqrt( (g1_vertices[i,1]-g2_vertices[j,1])**2+(g1_vertices[i,2]-g2_vertices[j,2])**2 )    
      #   if(sij<=5){sij=0}
      #   if(sij>5){sij=100000}
      
      if(f_edgeinclude){
        edgeedit_cost<-getEdgeEditCost(g1,i,g2,j,cid_e)
      }
      else{
        edgeedit_cost<-0
      }
      
      block1[i,j]<-sij+edgeedit_cost
    }#end of j loop
  }#end of i loop
  
  block2<-matrix(100000,N1,N1)
  for(i in 1:N1) 
    block2[i,i]<- cid_n+ (cid_e*length(get.neighborhood(g1,i)))
  
  block3<-matrix(100000,N2,N2)
  for(j in 1:N2)
    block3[j,j]<- cid_n + (cid_e*length(get.neighborhood(g2,j)))
  
  block4<-matrix(0,N2,N1)
  
  upper<-cbind(block1, block2)
  lower<-cbind(block3, block4)
  C<-rbind(upper,lower)
  
  x<-solve_LSAP(C, maximum=FALSE)
  edits <- cbind(seq_along(x), x)
  GED <- sum(C[edits])
  
  list(editPath=edits, editCost=GED)
  
}#end of graph edit path_twocids

########start ###################################################################
## This function calculates the edge cost incurred when a node in g1 is substituted 
## for a node in g2, node in g1 is deleted or a node in g2 is inserted.

getEdgeEditCost<-function(g1, i, g2, j, cid){

xvals <- get.vertex.attribute(g1, "xcoord")
yvals <- get.vertex.attribute(g1, "ycoord")
vertices_enr <- cbind(xvals, yvals)
xvals <- get.vertex.attribute(g2, "xcoord")
yvals <- get.vertex.attribute(g2, "ycoord")
vertices_que <- cbind(xvals, yvals)

neighbours1 <- get.neighborhood(g1,i)
neighbours2 <- get.neighborhood(g2,j)

E <- length(get.neighborhood(g1,i))
Q <- length(get.neighborhood(g2,j))

extra_cost<-0 #initialise

if (E==0 && Q==0) {extra_cost<-0}

if (E==0 && Q>0) {extra_cost<-cid*Q} # cost of inserting Q edges

if (E>0 && Q==0) {extra_cost<-cid*E} # cost of deleting these edges

if (E>0 && Q>0){

block1 <- array(0, c(E,Q))
for (i in 1:E){   
   for (j in 1:Q){

   n1 <- neighbours1[i]
   n2 <- neighbours2[j]
   sij <- sqrt((vertices_enr[n1,1]-vertices_que[n2,1])**2+(vertices_enr[n1,2]-vertices_que[n2,2])**2)
#    if(sij<=5){sij<-0}
#    if(sij>5){sij<-100000}
   block1[i,j]<- sij

}
}

block2 <- array(100000, c(E,E))
diag(block2) <- cid

block3 <- array(100000, c(Q,Q))
diag(block3) <- cid

block4 <- array(0, c(Q,E))

upper <- cbind(block1, block2)
lower <- cbind(block3, block4)
C <- rbind(upper, lower)

z <- solve_LSAP(C, maximum=FALSE)

extra_cost <- sum(C[cbind(seq_along(z), z)])
}

extra_cost



}#end of function

#####################################################################################################
########### Edit Path Using bounding box cost function


getGraphEditPath_bb<-function(g1,g2,cid,tol,f_edgeinclude){
  g1_vertices<-cbind(get.vertex.attribute(g1,"xcoord"),get.vertex.attribute(g1,"ycoord"))
  g2_vertices<-cbind(get.vertex.attribute(g2,"xcoord"),get.vertex.attribute(g2,"ycoord"))
  
  N1<-network.size(g1)
  N2<-network.size(g2)
  
  ##Build the cost matrix
  ############################
  
  block1<-matrix(0,N1, N2)
  
  for (i in 1:N1){
    for(j in 1:N2){
      
      sij<-sqrt( (g1_vertices[i,1]-g2_vertices[j,1])**2+(g1_vertices[i,2]-g2_vertices[j,2])**2 )    
        if(sij<=tol){sij=0}
        if(sij>tol){sij=100000}
      
      if(f_edgeinclude){
        edgeedit_cost<-getEdgeEditCost_bb(g1,i,g2,j,cid,tol)
      }
      else{
        edgeedit_cost<-0
      }
      
      block1[i,j]<-sij+edgeedit_cost
    }#end of j loop
  }#end of i loop
  
  block2<-matrix(100000,N1,N1)
  for(i in 1:N1) 
    block2[i,i]<- cid+ (cid*length(get.neighborhood(g1,i)))
  
  block3<-matrix(100000,N2,N2)
  for(j in 1:N2)
    block3[j,j]<- cid + (cid*length(get.neighborhood(g2,j)))
  
  block4<-matrix(0,N2,N1)
  
  upper<-cbind(block1, block2)
  lower<-cbind(block3, block4)
  C<-rbind(upper,lower)
  
  x<-solve_LSAP(C, maximum=FALSE)
  edits <- cbind(seq_along(x), x)
  GED <- sum(C[edits])
  
  list(editPath=edits, editCost=GED)
  
}#end of graph edit path
###################################################################################
########start ###################################################################
## This function calculates the edge cost incurred when a node in g1 is substituted 
## for a node in g2, node in g1 is deleted or a node in g2 is inserted.

getEdgeEditCost_bb<-function(g1, i, g2, j, cid, tol){
  
  xvals <- get.vertex.attribute(g1, "xcoord")
  yvals <- get.vertex.attribute(g1, "ycoord")
  vertices_enr <- cbind(xvals, yvals)
  xvals <- get.vertex.attribute(g2, "xcoord")
  yvals <- get.vertex.attribute(g2, "ycoord")
  vertices_que <- cbind(xvals, yvals)
  
  neighbours1 <- get.neighborhood(g1,i)
  neighbours2 <- get.neighborhood(g2,j)
  
  E <- length(get.neighborhood(g1,i))
  Q <- length(get.neighborhood(g2,j))
  
  extra_cost<-0 #initialise
  
  if (E==0 && Q==0) {extra_cost<-0}
  
  if (E==0 && Q>0) {extra_cost<-cid*Q} # cost of inserting Q edges
  
  if (E>0 && Q==0) {extra_cost<-cid*E} # cost of deleting these edges
  
  if (E>0 && Q>0){
    
    block1 <- array(0, c(E,Q))
    for (i in 1:E){   
      for (j in 1:Q){
        
        n1 <- neighbours1[i]
        n2 <- neighbours2[j]
        sij <- sqrt((vertices_enr[n1,1]-vertices_que[n2,1])**2+(vertices_enr[n1,2]-vertices_que[n2,2])**2)
           if(sij<=tol){sij<-0}
           if(sij>tol){sij<-100000}
        block1[i,j]<- sij
        
      }
    }
    
    block2 <- array(100000, c(E,E))
    diag(block2) <- cid
    
    block3 <- array(100000, c(Q,Q))
    diag(block3) <- cid
    
    block4 <- array(0, c(Q,E))
    
    upper <- cbind(block1, block2)
    lower <- cbind(block3, block4)
    C <- rbind(upper, lower)
    
    z <- solve_LSAP(C, maximum=FALSE)
    
    extra_cost <- sum(C[cbind(seq_along(z), z)])
  }
  
  extra_cost
  
  
  
}#end of function







######################################################################################################
############# START getmcs - takes an edit path and returns an mcs from two graphs ###################
## MCS is defined as a node induced subgraph of the query graph g2. 
## The MCS will have edges between the 2 nodes if corresponding nodes in g1 had an edge between them.
#######################################################################################################

getMCS <- function(g1, g2, edits){

E <- network.size(g1)
Q <- network.size(g2)


#create a new network for the mcs, with no edges, and same size as query
#========================================================================

mcs <- network.initialize(Q, directed=FALSE)


# now consider all edges in g1 and ADD an edge to the mcs if the corresponding nodes 
# in g2 are substituted as per the optimal assignment
#============================================================

edges <- as.matrix(g1, matrix.type="edgelist")
num_edges <- dim(edges)[1]
starts <- numeric()
ends <- numeric()

for (i in 1:num_edges) {
	edge <- edges[i,]
	base <- as.integer(edits[edge[1],2])
	tip <- as.integer(edits[edge[2],2])
           
	if (base<Q & tip<Q){

		possibles <- get.neighborhood(g2,base)
		degree <- length(possibles)

		if (degree>0) {

			flag <- 0
			for (j in 1:degree){
                    	if (possibles[j]==tip) {add.edge(mcs, base, tip)}
                   }

		}
	}
}


#paste in vertex attributes from g2 into mcs
#================================================


xvals <- get.vertex.attribute(g2, "xcoord")
yvals <- get.vertex.attribute(g2, "ycoord")
vertices <- cbind(xvals, yvals)

set.vertex.attribute(mcs, "xcoord", vertices[,1])
set.vertex.attribute(mcs, "ycoord", vertices[,2])



# delete from the mcs, nodes of g2 that have to be 
# inserted as per the optimal assignment
#====================================================

vertices_to_delete <- numeric()

for (i in 1:Q){

r <- E+i
fate <-edits[r,2]
if (fate<=Q) {vertices_to_delete <- c(vertices_to_delete, i)}
}

delete.vertices(mcs, vertices_to_delete)

if(network.size(mcs) > 0 ){
  #add edge and node attributes - this part newly added on 10/07/2013
  mcs_edges<-as.matrix.network(mcs, matrix.type="edgelist")
  if(length(mcs_edges[,1])>0){
    
    lengths<-rep(0,times=length(mcs_edges[,1]) )
    slopes<-rep(0,times=length(mcs_edges[,1]) )
    v<-cbind(get.vertex.attribute(mcs, "xcoord") , get.vertex.attribute(mcs, "ycoord") )
    
    for (i in 1:length(mcs_edges[,1])){
      
      start <- mcs_edges[i,1]
      end <- mcs_edges[i,2]
      
      lengths[i] <- sqrt((v[start,1]-v[end,1])**2+(v[start,2]-v[end,2])**2)
      
      if((v[end,1]-v[start,1])==0){slopes[i]=90}
      if((v[end,1]-v[start,1])!=0){
        if (v[start,1]<v[end,1]) slopes[i] <- (180/pi)*atan((v[end,2]-v[start,2])/(v[end,1]-v[start,1]))
        if (v[start,1]>=v[end,1]) slopes[i] <- (180/pi)*atan((v[start,2]-v[end,2])/(v[start,1]-v[end,1]))
        
      }
      
      lengths[i]<-round(lengths[i],2)
      slopes[i]<-round(slopes[i],2)
      #print(i)
    }
    
    set.edge.attribute(mcs,"edgelength",lengths)
    set.edge.attribute(mcs,"edgeslope",slopes)
    
    
    
  }
}


mcs

}#end of function
##########################################################################################################

############################################################################################################
# START getDist - function that computes several types of distances between 2 graphs and returns that list

getDist<-function(g1,g2,mcs,cid){

g1_nodes<-network.size(g1)
g1_edges<-network.edgecount(g1)
g2_nodes<-network.size(g2)
g2_edges<-network.edgecount(g2)
mcs_nodes<-network.size(mcs)
mcs_edges<-network.edgecount(mcs)

NodeInsDelCost<-cid



#Common subraph distance
#-------------------------
d_1_infty_nodes<-max(g1_nodes,g2_nodes)-mcs_nodes
d_1_infty_edges<-max(g1_edges,g2_edges)-mcs_edges
dbar_1_infty_nodes<-d_1_infty_nodes/max(g1_nodes,g2_nodes)
dbar_1_infty_edges<-d_1_infty_edges/max(g1_edges,g2_edges)

#Bunke Miro combo and weighted combo= d(1,1)nodes+d(1,1)edges
#---------------------------------------------------------------
BM<-g1_nodes+g2_nodes-(2*mcs_nodes)+g1_edges+g2_edges-(2*mcs_edges)
# in weighted assime c(x)=x, edge subs cost=0. Edge deletion cost, edge insertion cost determined by NodeInsDelCost
#---------------------------------------------------------------------------------------------------------------------
BMweighted<-g1_nodes+g2_nodes-(2*mcs_nodes)+0+((g1_edges-mcs_edges)*NodeInsDelCost)+((g2_edges-mcs_edges)*NodeInsDelCost)


#Other distance measures from the max min family
#-------------------------------------------------
d_1_2_nodes<-d_1_infty_nodes+(0.5*(min(g1_nodes,g2_nodes)-mcs_nodes ))
d_1_2_edges<-d_1_infty_edges+(0.5*(min(g1_edges,g2_edges)-mcs_edges))
dbar_1_2_nodes<-d_1_2_nodes/(mcs_nodes+d_1_2_nodes)
dbar_1_2_edges<-d_1_2_edges/(mcs_edges+d_1_2_edges)

d_1_1_nodes<-d_1_infty_nodes+(1*(min(g1_nodes,g2_nodes)-mcs_nodes ))
d_1_1_edges<-d_1_infty_edges+(1*(min(g1_edges,g2_edges)-mcs_edges))
dbar_1_1_nodes<-d_1_1_nodes/(mcs_nodes+d_1_1_nodes)
dbar_1_1_edges<-d_1_1_edges/(mcs_edges+d_1_1_edges)


#Other distance measures from the Minkowski family
#-------------------------------------------------
d_2_2_nodes<-sqrt((g1_nodes-mcs_nodes)^2+(g2_nodes-mcs_nodes)^2)
d_2_2_edges<-sqrt((g1_edges-mcs_edges)^2+(g2_edges-mcs_edges)^2)
dbar_2_2_nodes<-d_2_2_nodes/(mcs_nodes+d_2_2_nodes)
dbar_2_2_edges<-d_2_2_edges/(mcs_edges+d_2_2_edges)

#Non metric distance measure (geometric mean based)
#---------------------------------------------------
d_sqrt_nodes<-sqrt(g1_nodes*g2_nodes)- mcs_nodes
d_sqrt_edges<-sqrt(g1_edges*g2_edges)- mcs_edges
dbar_sqrt_nodes<-d_sqrt_nodes/sqrt(g1_nodes*g2_nodes)
dbar_sqrt_edges<-d_sqrt_edges/sqrt(g1_edges*g2_edges)


list(d_1_infty_n=d_1_infty_nodes, d_1_infty_e=d_1_infty_edges, dbar_1_infty_n=dbar_1_infty_nodes, dbar_1_infty_e=dbar_1_infty_edges, bm=BM, bmweighted=BMweighted, d_1_2_n=d_1_2_nodes, d_1_2_e=d_1_2_edges, dbar_1_2_n=dbar_1_2_nodes, dbar_1_2_e=dbar_1_2_edges, d_1_1_n=d_1_1_nodes, d_1_1_e=d_1_1_edges, dbar_1_1_n=dbar_1_1_nodes, dbar_1_1_e=dbar_1_1_edges, d_2_2_n=d_2_2_nodes, d_2_2_e=d_2_2_edges, dbar_2_2_n=dbar_2_2_nodes, dbar_2_2_e=dbar_2_2_edges, g1_n=g1_nodes, g1_e=g1_edges, g2_n=g2_nodes, g2_e=g2_edges, mcs_n=mcs_nodes, mcs_e=mcs_edges, d_sqrt_n=d_sqrt_nodes, d_sqrt_e=d_sqrt_edges, dbar_sqrt_n=dbar_sqrt_nodes, dbar_sqrt_e=dbar_sqrt_edges)
}# end of function




##########################################################################################################
# START getGenImpScores -
# This function takes the biometric (biometric->FP,BP,RE,HV), sourceid (source->AutoJJ, AutoSL, Manual),
# indexfilelist - filename of a list of person ids , sample ids and filenames in a 3 column tab limited file. 
# Other inputs are :
# tol- tolerance for rough alignment
# cid - insertion deletion cost in Munkre's algorithm
# f_edgeinclude - flag to include (1) or exclude(0) edge costs in Munkre's algorithm
# It outputs two column vectors - gen, listing genuine scores and imp, listing imposter scores
############################################################################################################


getGenImpScores<-function(biometric,sourceid,indexfilelist,tol,cid,f_edgeinclude,f_plot){
genuine<-numeric()
imposter<-numeric()

filelist<-read.table(indexfilelist,header=FALSE, sep="\t", fill=TRUE)
filelist<-as.matrix(filelist)

DBSize<-nrow(filelist)

genscorematrix<-numeric()
impscorematrix<-numeric()

for(i in 1:(DBSize-1)){
	person1<-filelist[i,1]
	sample1<-filelist[i,2]
	G1<-getGraph(biometric,sourceid,person1,sample1)

	for(j in (i+1):DBSize){
		person2<-filelist[j,1]
		sample2<-filelist[j,2]
		G2<-getGraph(biometric,sourceid,person2,sample2)

		RG<-registerGraphs(G1,G2,tol,f_plot)
		Path<-getGraphEditPath(RG$graph1,RG$graph2,cid,f_edgeinclude)
		MCS<-getMCS(RG$graph1,RG$graph2,Path$editPath)
		D<-getDist(RG$graph1,RG$graph2,MCS,cid) 
		if(person1==person2){
		genuine<-c(genuine,D$dbar_sqrt_n)
		genscorematrix<-rbind(genscorematrix,cbind(person1,sample1,person2,sample2,D$dbar_sqrt_n))
		}
		if(person1!=person2 && sample1==1 && sample2==1){
		imposter<-c(imposter,D$dbar_sqrt_n)
		impscorematrix<-rbind(impscorematrix,cbind(person1,sample1,person2,sample2,D$dbar_sqrt_n))
		}
print(paste(person1,sample1,"vs",person2,sample2,D$dbar_sqrt_n,sep=" "))
op_row<-cbind(person1,sample1,person2,sample2,D$dbar_sqrt_n)
write.table(op_row, file="H:/Uni Stuff/Uni Stuff/2012/PostDoc/Research Projects/ARC 2012-2015/R Code/Outputs/WithBGM.txt",append=TRUE, sep=" ",row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
		
#print(i)
#print("vs")
#print(j)
	}# end of j loop
}#end of i loop

list(gen=genuine,imp=imposter,gsm=genscorematrix,ism=impscorematrix)
}# end of function
###########################################################################################################
# featureNumber : 1 to 16 based on what feature of mcs you want to compare
##########################################################################


###########################################################################################################
############# START getEER - gets an equal error rate from two sets of normalised distances ###############
getEER <- function(genuines, imposters, genscorematrix,impscorematrix){

x <- rep(0, times=999)
y <- rep(0, times=999)
w <- rep(0, times=999)#stores the operating threshold

for (i in 1:999) {
Trues <- genuines
Fakes <- imposters
p <- i/1000
w[i]<-p


x[i] <- length(which(Fakes<=p, arr.ind=TRUE))/length(Fakes)
#x[i] <- 1-sum(Fakes)/dim(Fakes)[1]

#FRR is the proportion of T that is below p

y[i] <- length(which(Trues>p, arr.ind=TRUE))/length(Trues)
}
#plot(x,y, xlab="FAR", ylab="FRR")

z <- abs(x-y)
intersection <- order(z)[1]

EER<-0.5*(x[intersection]+y[intersection])
threshold<-w[intersection]

#now get list of genuine comparisons that had score above the threshold and imposters who had score below threshold.
generrors<-numeric()
imperrors<-numeric()

for(i in 1:length(genuines)){
  if(genuines[i]>threshold) generrors<-rbind(generrors,genscorematrix[i,])
  
}#end of i loop

for(i in 1:length(imposters)){
  if(imposters[i]<=threshold) imperrors<-rbind(imperrors,impscorematrix[i,])
}


list(eer=EER,th=threshold,genrejects=generrors,impaccepts=imperrors)
}
############# END getEER - gets an equal error rate from two sets of normalised distances #################
###########################################################################################################

############# START getEER - gets an equal error rate from two sets of normalised distances ###############
getEER_nomx <- function(genuines, imposters){
  
  x <- rep(0, times=999)
  y <- rep(0, times=999)
  w <- rep(0, times=999)#stores the operating threshold
  
  for (i in 1:999) {
    Trues <- genuines
    Fakes <- imposters
    p <- i/1000
    w[i]<-p
    
    
    x[i] <- length(which(Fakes<=p, arr.ind=TRUE))/length(Fakes)
    #x[i] <- 1-sum(Fakes)/dim(Fakes)[1]
    
    #FRR is the proportion of T that is below p
    
    y[i] <- length(which(Trues>p, arr.ind=TRUE))/length(Trues)
  }
  #plot(x,y, xlab="FAR", ylab="FRR")
  
  z <- abs(x-y)
  intersection <- order(z)[1]
  
  EER<-0.5*(x[intersection]+y[intersection])
  threshold<-w[intersection]
  
  
  
  list(eer=EER,th=threshold)
}
############# END getEER - gets an equal error rate from two sets of normalised distances #################
###########################################################################################################


################START getROC
getROC<-function(genuines, imposters, genscorematrix, impscorematrix, f_plot){

  x <- rep(0, times=999)
  y <- rep(0, times=999)
  w <- rep(0, times=999)#stores the operating threshold
  
  for (i in 1:999) {
    Trues <- genuines
    Fakes <- imposters
    p <- i/1000
    w[i]<-p
    
#     #FAR is the proportion of F that is below p
#     F[F>p] <- 1
#     F[F<=p] <- 0
#     x[i] <- 1-sum(F)/length(F)
    x[i] <- length(which(Fakes<=p, arr.ind=TRUE))/length(Fakes)
    
#     #FRR is the proportion of T that is above p
#     T[T>p] <- 1
#     T[T<=p] <- 0
#     y[i] <- sum(T)/length(T)
    y[i] <- length(which(Trues>p, arr.ind=TRUE))/length(Trues)
  }
  
  ROCmx<-cbind(w,100*y,100*x)
  
  if(f_plot==1){
    plot(ROCmx[,3],ROCmx[,2], type="l", xlab="FAR", ylab="FNMR")
   title("ROC")
    abline(a=0, b=1, untf=FALSE)
    abline(v=0, col="red")
  
  }
  
  ROCmx
  
}# end of function




##################END getROC


############# START getEER_unnorm - gets an equal error rate from two sets of unnormalised distances ###############
getEER_unnorm<- function(genuines, imposters, genscorematrix,impscorematrix){
  
  range(genuines)
  range(imposters)
  thMin<-min( c( range(genuines), range(imposters) ) )
  thMax<-max( c( range(genuines), range(imposters) ) )
  
  x <- rep(0, times=999)
  y <- rep(0, times=999)
  w <- rep(0, times=999)#stores the operating threshold
  
  for (i in 1:999) {
    Trues <- genuines
    Fakes <- imposters
    p <- (i/1000)*(thMax-thMin)
    p<-thMin+p
    w[i]<-p
    
    #set all values below p will be accepted
    
    #FAR (unnorm) is the proportion of F that is below p
   
    x[i] <- length(which(Fakes<=p, arr.ind=TRUE))/length(Fakes)
    #x[i] <- 1-sum(Fakes)/dim(Fakes)[1]
    
    #FRR is the proportion of T that is below p
    
    y[i] <- length(which(Trues>p, arr.ind=TRUE))/length(Trues)
  }
  #plot(x,y, xlab="FAR", ylab="FRR")
  
  z <- abs(x-y)
  intersection <- order(z)[1]
  
  EER<-0.5*(x[intersection]+y[intersection])
  threshold<-w[intersection]
  
  #now get list of genuine comparisons that had score above the threshold and imposters who had score below threshold.
  generrors<-numeric()
  imperrors<-numeric()
  
  for(i in 1:length(genuines)){
    if(genuines[i]>threshold) generrors<-rbind(generrors,genscorematrix[i,])
   
  }#end of i loop
  
  for(i in 1:length(imposters)){
    if(imposters[i]<=threshold) imperrors<-rbind(imperrors,impscorematrix[i,])
  }
  
  list(eer=EER,th=threshold,genrejects=generrors,impaccepts=imperrors)
}
############# END getEER - gets an equal error rate from two sets of normalised distances #################
###########################################################################################################

############# START getEER_unnorm - gets an equal error rate from two sets of unnormalised similarities ###############
getEER_sim<- function(genuines, imposters, genscorematrix,impscorematrix){
  
#   range(genuines)
#   range(imposters)
  thMin<-min( c( range(genuines), range(imposters) ) )
  thMax<-max( c( range(genuines), range(imposters) ) )
  
  x <- rep(0, times=999)
  y <- rep(0, times=999)
  w <- rep(0, times=999)#stores the operating threshold
  
  for (i in 1:999) {
    Trues <- genuines
    Fakes <- imposters
    p <- (i/1000)*(thMax-thMin)
    p<-thMin+p
    w[i]<-p
    
    #set all values below p will be accepted
    
    #FAR (unnorm) is the proportion of F that is below p
    
    x[i] <- length(which(Fakes>p, arr.ind=TRUE))/length(Fakes)
    #x[i] <- 1-sum(Fakes)/dim(Fakes)[1]
    
    #FRR is the proportion of T that is below p
    
    y[i] <- length(which(Trues<=p, arr.ind=TRUE))/length(Trues)
  }
  #plot(x,y, xlab="FAR", ylab="FRR")
  
  z <- abs(x-y)
  intersection <- order(z)[1]
  
  EER<-0.5*(x[intersection]+y[intersection])
  threshold<-w[intersection]
  
  #now get list of genuine comparisons that had score below the threshold and imposters who had score above threshold.
  generrors<-numeric()
  imperrors<-numeric()
  
  for(i in 1:length(genuines)){
    if(genuines[i]<=threshold) generrors<-rbind(generrors,genscorematrix[i,])   
  }#end of i loop
  for(i in 1:length(imposters)){
    if(imposters[i]>threshold) imperrors<-rbind(imperrors,impscorematrix[i,])
  }
  
  list(eer=EER,th=threshold,genrejects=generrors,impaccepts=imperrors)
}
############# END getEER - gets an equal error rate from two sets of normalised distances #################
###########################################################################################################

getError<-function(genuines, imposters, genscorematrix,impscorematrix, Th){
  
  #now get list of genuine comparisons that had score above the threshold and imposters who had score below threshold.
  generrors<-numeric()
  imperrors<-numeric()
  FNM<-0
  FM<-0
  FNMR<-0
  FMR<-0
  TE<-0
  
  gen_ctr<-0
  imp_ctr<-0
  for(i in 1:length(genuines)){
    if(genuines[i]>Th) {generrors<-rbind(generrors,genscorematrix[i,])
                        gen_ctr<-gen_ctr+1 }
    
  }#end of i loop
  
  for(i in 1:length(imposters)){
    if(imposters[i]<=Th) { imperrors<-rbind(imperrors,impscorematrix[i,])
                           imp_ctr<-imp_ctr+1}
    
  }
  
  if(gen_ctr>0){
    FNM<-dim(generrors)[1]
    FNMR<-FNM/dim(genscorematrix)[1] *100
    
  } 
    
  if(imp_ctr>0){
    FM<-dim(imperrors)[1]
    FMR<-FM/dim(impscorematrix)[1] *100
  } 
  
  
  TE<- (FNM+FM) / ( dim(genscorematrix)[1]+dim(impscorematrix)[1] ) * 100
  list(FNMR=FNMR, FMR=FMR, TE=TE, Th=Th, gsm=generrors, ism=imperrors)
}

###################################################################################
########### BEGIN getGraphEditPath_slack ################################################
## This function takes 2 ALIGNED graphs and cid, the insertion/deletion cost as input
## It has a fourth input, a flag that selects if edges must be included in deciding the edit path (1 include, 0 exculde)
## It outputs the cheapest graph edit path in going from g1 to g2
## The cost matrix uses euclidean distance between vertex attributes as the basis 
## for computing cost
## It is also slack in that it encourages imposters to match.


getGraphEditPath_slack<-function(g1,g2,cid,f_edgeinclude){
  g1_vertices<-cbind(get.vertex.attribute(g1,"xcoord"),get.vertex.attribute(g1,"ycoord"))
  g2_vertices<-cbind(get.vertex.attribute(g2,"xcoord"),get.vertex.attribute(g2,"ycoord"))
  
  N1<-network.size(g1)
  N2<-network.size(g2)
  
  ##Build the cost matrix
  ############################
  
  block1<-matrix(0,N1, N2)
  
  for (i in 1:N1){
    for(j in 1:N2){
      
      sij<-sqrt( (g1_vertices[i,1]-g2_vertices[j,1])**2+(g1_vertices[i,2]-g2_vertices[j,2])**2 )    
        if(sij<=(cid*3)){sij=0}
        #if(sij>(cid*3)){sij=100000}
      
      if(f_edgeinclude){
        edgeedit_cost<-getEdgeEditCost_slack(g1,i,g2,j,cid)
      }
      else{
        edgeedit_cost<-0
      }
      
      block1[i,j]<-sij+edgeedit_cost
    }#end of j loop
  }#end of i loop
  
  block2<-matrix(100000,N1,N1)
  for(i in 1:N1) 
    block2[i,i]<- cid+ (cid*length(get.neighborhood(g1,i)))
  
  block3<-matrix(100000,N2,N2)
  for(j in 1:N2)
    block3[j,j]<- cid + (cid*length(get.neighborhood(g2,j)))
  
  block4<-matrix(0,N2,N1)
  
  upper<-cbind(block1, block2)
  lower<-cbind(block3, block4)
  C<-rbind(upper,lower)
  
  x<-solve_LSAP(C, maximum=FALSE)
  edits <- cbind(seq_along(x), x)
  GED <- sum(C[edits])
  
  list(editPath=edits, editCost=GED)
  
}#end of graph edit path
###################################################################################
########start ###################################################################
## This function calculates the edge cost incurred when a node in g1 is substituted 
## for a node in g2, node in g1 is deleted or a node in g2 is inserted.

getEdgeEditCost_slack<-function(g1, i, g2, j, cid){
  
  xvals <- get.vertex.attribute(g1, "xcoord")
  yvals <- get.vertex.attribute(g1, "ycoord")
  vertices_enr <- cbind(xvals, yvals)
  xvals <- get.vertex.attribute(g2, "xcoord")
  yvals <- get.vertex.attribute(g2, "ycoord")
  vertices_que <- cbind(xvals, yvals)
  
  neighbours1 <- get.neighborhood(g1,i)
  neighbours2 <- get.neighborhood(g2,j)
  
  E <- length(get.neighborhood(g1,i))
  Q <- length(get.neighborhood(g2,j))
  
  extra_cost<-0 #initialise
  
  if (E==0 && Q==0) {extra_cost<-0}
  
  if (E==0 && Q>0) {extra_cost<-cid*Q} # cost of inserting Q edges
  
  if (E>0 && Q==0) {extra_cost<-cid*E} # cost of deleting these edges
  
  if (E>0 && Q>0){
    
    block1 <- array(0, c(E,Q))
    for (i in 1:E){   
      for (j in 1:Q){
        
        n1 <- neighbours1[i]
        n2 <- neighbours2[j]
        sij <- sqrt((vertices_enr[n1,1]-vertices_que[n2,1])**2+(vertices_enr[n1,2]-vertices_que[n2,2])**2)
            if(sij<=(cid*3)){sij<-0}
        #    if(sij>5){sij<-100000}
        block1[i,j]<- sij
        
      }
    }
    
    block2 <- array(100000, c(E,E))
    diag(block2) <- cid
    
    block3 <- array(100000, c(Q,Q))
    diag(block3) <- cid
    
    block4 <- array(0, c(Q,E))
    
    upper <- cbind(block1, block2)
    lower <- cbind(block3, block4)
    C <- rbind(upper, lower)
    
    z <- solve_LSAP(C, maximum=FALSE)
    
    extra_cost <- sum(C[cbind(seq_along(z), z)])
  }
  
  extra_cost
  
  
  
}#end of function

#####################################################################################################
########### Edit Path Using bounding box cost function

################START getROC
getROC_nomx<-function(genuines, imposters, f_plot){
  
  x <- rep(0, times=999)
  y <- rep(0, times=999)
  w <- rep(0, times=999)#stores the operating threshold
  
  for (i in 1:9999) {
    Trues <- genuines
    Fakes <- imposters
    p <- i/10000
    w[i]<-p
    
    #     #FAR is the proportion of F that is below p
    #     F[F>p] <- 1
    #     F[F<=p] <- 0
    #     x[i] <- 1-sum(F)/length(F)
    x[i] <- length(which(Fakes<=p, arr.ind=TRUE))/length(Fakes)
    
    #     #FRR is the proportion of T that is above p
    #     T[T>p] <- 1
    #     T[T<=p] <- 0
    #     y[i] <- sum(T)/length(T)
    y[i] <- length(which(Trues>p, arr.ind=TRUE))/length(Trues)
  }
  
  ROCmx<-cbind(w,100*y,100*x)
  
  if(f_plot==1){
    windows()
    pts<-cbind(ROCmx[,3],ROCmx[,2])
    plot(lowess(pts), type="b", xlab="FMR", ylab="FNMR")
    lines(lowess(pts), col=3)
    abline(a=0, b=1, untf=FALSE)
    abline(v=0, col="red")
    
  }
  
  ROCmx
  
}# end of function




##################END getROC

############# START getEER_unnorm - gets an equal error rate from two sets of unnormalised distances ###############
getEER_unnorm_nomx<- function(genuines, imposters){
  
  range(genuines)
  range(imposters)
  thMin<-min( c( range(genuines), range(imposters) ) )
  thMax<-max( c( range(genuines), range(imposters) ) )
  
  x <- rep(0, times=999)
  y <- rep(0, times=999)
  w <- rep(0, times=999)#stores the operating threshold
  
  for (i in 1:999) {
    Trues <- genuines
    Fakes <- imposters
    p <- (i/1000)*(thMax-thMin)
    p<-thMin+p
    w[i]<-p
    
    #set all values below p will be accepted
    
    #FAR (unnorm) is the proportion of F that is below p
    
    x[i] <- length(which(Fakes<=p, arr.ind=TRUE))/length(Fakes)
    #x[i] <- 1-sum(Fakes)/dim(Fakes)[1]
    
    #FRR is the proportion of T that is below p
    
    y[i] <- length(which(Trues>p, arr.ind=TRUE))/length(Trues)
  }
  #plot(x,y, xlab="FAR", ylab="FRR")
  
  z <- abs(x-y)
  intersection <- order(z)[1]
  
  EER<-0.5*(x[intersection]+y[intersection])
  threshold<-w[intersection]
  
  
  
  list(eer=EER,th=threshold)
}
############# END getEER - gets an equal error rate from two sets of normalised distances #################
###########################################################################################################

############# START getEER_unnorm - gets an equal error rate from two sets of unnormalised similarities ###############
getEER_sim_nomx<- function(genuines, imposters){
  
  #   range(genuines)
  #   range(imposters)
  thMin<-min( c( range(genuines), range(imposters) ) )
  thMax<-max( c( range(genuines), range(imposters) ) )
  
  x <- rep(0, times=999)
  y <- rep(0, times=999)
  w <- rep(0, times=999)#stores the operating threshold
  
  for (i in 1:999) {
    Trues <- genuines
    Fakes <- imposters
    p <- (i/1000)*(thMax-thMin)
    p<-thMin+p
    w[i]<-p
    
    #set all values below p will be accepted
    
    #FAR (unnorm) is the proportion of F that is below p
    
    x[i] <- length(which(Fakes>p, arr.ind=TRUE))/length(Fakes)
    #x[i] <- 1-sum(Fakes)/dim(Fakes)[1]
    
    #FRR is the proportion of T that is below p
    
    y[i] <- length(which(Trues<=p, arr.ind=TRUE))/length(Trues)
  }
  #plot(x,y, xlab="FAR", ylab="FRR")
  
  z <- abs(x-y)
  intersection <- order(z)[1]
  
  EER<-0.5*(x[intersection]+y[intersection])
  threshold<-w[intersection]
  
  list(eer=EER,th=threshold)
}
############# END getEER - gets an equal error rate from two sets of unnormalised similarities #################

############# START getROC_unnorm - gets an equal error rate from two sets of unnormalised distances ###############
getROC_unnorm_nomx<- function(genuines, imposters, fplot){
  
  range(genuines)
  range(imposters)
  thMin<-min( c( range(genuines), range(imposters) ) )
  thMax<-max( c( range(genuines), range(imposters) ) )
  
  x <- rep(0, times=999)
  y <- rep(0, times=999)
  w <- rep(0, times=999)#stores the operating threshold
  
  for (i in 1:999) {
    Trues <- genuines
    Fakes <- imposters
    p <- (i/1000)*(thMax-thMin)
    p<-thMin+p
    w[i]<-p
    
    #set all values below p will be accepted
    
    #FAR (unnorm) is the proportion of F that is below p
    
    x[i] <- length(which(Fakes<=p, arr.ind=TRUE))/length(Fakes)
    #x[i] <- 1-sum(Fakes)/dim(Fakes)[1]
    
    #FRR is the proportion of T that is below p
    
    y[i] <- length(which(Trues>p, arr.ind=TRUE))/length(Trues)
  }
 
  
  z <- abs(x-y)
  intersection <- order(z)[1]
  
  EER<-0.5*(x[intersection]+y[intersection])
  threshold<-w[intersection]
  
  if(fplot==1){
    windows()
    #plot(x,y, type="l", lty=1,xlab="FAR", ylab="FNMR", lwd=1)
    plot(log(x),log(y), type="l", lty=1,xlab="FAR (log scale)", ylab="FNMR (log scale)", lwd=1)
    abline(a=0, b=1, col="green", lwd=2)
    abline(h=EER, col="red", lty=2)
    abline(v=EER, col="red", lty=2)
    textxy(X=EER, Y=EER, labs=round(EER, digits=4), offset=0.8, cex=1)
  }
  
  
  list(w=w, x=x, y=y, eer=EER,th=threshold)
}
############# END getEER - gets an equal error rate from two sets of normalised distances #################
###########################################################################################################

############# START getEER_unnorm - gets an equal error rate from two sets of unnormalised similarities ###############
getROC_sim_nomx<- function(genuines, imposters, fplot){
  
  #   range(genuines)
  #   range(imposters)
  thMin<-min( c( range(genuines), range(imposters) ) )
  thMax<-max( c( range(genuines), range(imposters) ) )
  
  x <- rep(0, times=999)
  y <- rep(0, times=999)
  w <- rep(0, times=999)#stores the operating threshold
  
  for (i in 1:999) {
    Trues <- genuines
    Fakes <- imposters
    p <- (i/1000)*(thMax-thMin)
    p<-thMin+p
    w[i]<-p
    
    #set all values below p will be accepted
    
    #FAR (unnorm) is the proportion of F that is below p
    
    x[i] <- length(which(Fakes>p, arr.ind=TRUE))/length(Fakes)
    #x[i] <- 1-sum(Fakes)/dim(Fakes)[1]
    
    #FRR is the proportion of T that is below p
    
    y[i] <- length(which(Trues<=p, arr.ind=TRUE))/length(Trues)
  }
  #plot(x,y, xlab="FAR", ylab="FRR")
  
  z <- abs(x-y)
  intersection <- order(z)[1]
  
  EER<-0.5*(x[intersection]+y[intersection])
  threshold<-w[intersection]
  
  if(fplot==1){
    windows()
    plot(x,y, type="l", lty=1,xlab="FAR", ylab="FRR", lwd=1)
    abline(a=0, b=1, col="green", lwd=2)
    abline(h=EER, col="red", lty=2)
    abline(v=EER, col="red", lty=2)
    textxy(X=EER, Y=EER, labs=round(EER, digits=4), offset=0.8, cex=1)
    
  }
  
  list(w=w, x=x, y=y, eer=EER,th=threshold)
}
############# END getEER - gets an equal error rate from two sets of unnormalised similarities #################

###########plothist#########################################

plothist<-function(genuines, imposters){
  
  gen_hist<-hist(genuines, breaks=10, plot=FALSE )
  imp_hist<-hist(imposters, breaks=10, plot=FALSE)
  eer<-getEER_nomx(genuines,imposters)
  #windows()
  plot(gen_hist$mids, gen_hist$counts/sum(gen_hist$counts)*100, xlim=c(0,max( max(genuines), max(imposters) )), ylim=c(0,100), col="green", type='b', pch=2, xlab="scores", ylab="distribution")
  lines(imp_hist$mids, imp_hist$counts/sum(imp_hist$counts)*100, col="red", type='b', pch=2)
  title(main=paste("EER =", eer$eer, sep=" "))
  
  
}


##########end plothist function #############################
