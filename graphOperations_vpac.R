#This file defines operations between two noisy spatial graphs
# Functions defined are
# Graph Union

#This function takes two graphs, aligns them identifies vertices of g2 that are within tol of g1, and outputs the union.
noisySpatialUnion<-function(g1,g2,tol, maxtheta, f_plot){
  RG<-registerGraphs_maxtheta(g1,g2,tol,maxtheta,f_plot)
  #take only alignment on the first best pair of edges
 
  g1<-RG$graph1.g_aligned[[1]]
  g2<-RG$graph2.g_aligned[[1]]
#   windows()
#   plotspatialgraph(g1)
#   windows()
#   plotspatialgraph(g2)
  set.vertex.attribute(g1, "vertex.names", 1:network.size(g1))
  set.vertex.attribute(g2, "vertex.names", 1:network.size(g2))
  
  RS<-rapidscore_matchlist(g1, g2, tol)
  
  v1<-cbind( get.vertex.attribute(g1, "xcoord"), get.vertex.attribute(g1, "ycoord") )
  v2<-cbind( get.vertex.attribute(g2, "xcoord"), get.vertex.attribute(g2, "ycoord") )
  v_new<-rbind(v1,v2)
  a1<-as.matrix.network(g1, matrix.type="adjacency")
  a2<-as.matrix.network(g2, matrix.type="adjacency")
  a12<-matrix(0, nrow=dim(v1)[1], ncol=dim(v2)[1])
  a21<-matrix(0, nrow=dim(v2)[1], ncol=dim(v1)[1])
  
  a_new<-rbind( as.matrix( cbind(a1,a12)) , cbind(a21,a2) )
   
  #connect edges of g2 to matching vertices of g1
  listmatchedvertices<-cbind(RS$matchlist[,1], RS$matchlist[,2]+network.size(g1) )
  for(i in 1:dim(listmatchedvertices)[1] ){
    g2row<-listmatchedvertices[i,2]
    g1row<-listmatchedvertices[i,1]
    g2col<- which(a_new[g2row,] !=0 )
    a_new[g1row,g2col]<-1
    
    
    
  }
  g_new<-as.network.matrix(a_new, matrix.type="adjacency")
  set.vertex.attribute(g_new, "xcoord", v_new[,1])
  set.vertex.attribute(g_new, "ycoord", v_new[,2])
  #delete the matched vertices of g2, to avoid vertices very close ti each other in the new graph.
  delete.vertices(g_new, listmatchedvertices[,2])
  
  #set the slope and length attributes
  edges<-as.matrix.network(g_new, matrix.type="edgelist")
  newervertices<-cbind( get.vertex.attribute(g_new, "xcoord") , get.vertex.attribute(g_new, "ycoord") )
  
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
  
  set.edge.attribute(g_new,"edgelength",lengths)
  set.edge.attribute(g_new,"edgeslope",slopes)
  
  
  if(f_plot==1){
    windows()
    plotspatialgraph(g_new)
    title("Noisy Union Graph")
    
  }
 
  g_new
#   
}

###########
# distance matrix between graphs and prototypes
compute_tpmatrix<-function(d_matrix_t, TrainingIndices, PrototypeIndices, lessthan3ind_tokeep){
  
  tp_matrix<-matrix(0, nrow=length(TrainingIndices), ncol=length(PrototypeIndices))
  colInd<-numeric() #column indices for the prototypes
  for(j in 1:dim(tp_matrix)[2]){
    pInd<-PrototypeIndices[j]
    col<-which(lessthan3ind_tokeep==pInd)   
    colInd<-c(colInd,col)
  }#end of j loop
    temp<-as.matrix(d_matrix_t[,colInd]) #select the columns
    #remove the rows corresponding to all the prototypes
    temp<-temp[-colInd,] 
    tp_matrix<-as.matrix(temp)
    
 
  
  tp_matrix
  
  
}#end of function
