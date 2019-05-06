#required packages
#=================
#setwd("I:/Academic/Research/Projects and Grants/Biometric Databases/Graph Databases")

#setwd("D:/Arathi/Graph Databases") #copied graphs to local machine to increase speed
library(splancs)
library(spatial)
library(network)
library(clue)
library(iterators)
library(foreach)
library(doSNOW)

source("basicGraphFunctions_vpac.R")
source("graphMatchingFunctions_vpac.R")



#Functions defined in this file
#================================
# registerGraphs_maxtheta<-function(g1,g2,tol,maxtheta, f_plot) -- registers based on a max bound for rotation of frame
# registerGraphs_hungarian<-function(g1,g2,cid,f_plot,f_edgeinclude)

#============================================================================

#############START registerGraphs_maxtheta  , maxtheta is the maximum rotation allowed for a pair of matching edges
#######################################################
registerGraphs_maxtheta<-function(g1,g2,tol,maxtheta, f_plot){
  
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
      #edge_matching[i,j] <- edge_matching[i,j]/(0.5*g1_edgelength[i]+0.5*g2_edgelength[j])
      #edge_matching[i,j] <- edge_matching[i,j]/(d1+d2)
      
      
    }
  }
  
  ranking <- order(edge_matching, decreasing=FALSE)
  
  
  #take first 2500 closest edge pairs to realign and get best match score
  #================================================================
  #NEdgePairs<-1000 #for SFIR, NEdgePairs=1000 for SNIR
  #NEdgePairs<-length(ranking)
  NEdgePairs<-max(10, round(0.02*length(ranking)) )
  if(length(ranking)<NEdgePairs) { NEdgePairs<- length(ranking) }
  bestMatchScore<-numeric()
  bestAlignEdge<-numeric()
#   g1_shiftrot<-g1
#   g2_shiftrot<-g2
  g1_shiftrot<-numeric()
  g2_shiftrot<-numeric()
  scores<-numeric()
  alignedEdge<-numeric()
  bestRank<-numeric()
  bestIndex<-numeric()
 
 
  
#   outputlist<-numeric()
#   cl<-makeCluster(6, type = "SOCK")
#   #clusterExport(cl, c("main.fun", "helper1", "helper2"))
#   registerDoSNOW(cl)
  #registerDoMC(8)
  outputlist<-foreach(i=1:NEdgePairs, .combine="cbind", .packages=c("network","sna", "clue")  )%dopar%{
   # outputlist<-foreach(i=1:NEdgePairs, .combine="cbind", .packages=c("network","sna", "clue")  )%do%{
    source("basicGraphFunctions_vpac.R")
 #for(i in 1:NEdgePairs){
    
    rank<-i
    
    column <- ceiling(ranking[rank]/length(edge_matching[,1]))
    row <- ranking[rank]%%length(edge_matching[,1])
    if (row==0) row <- length(edge_matching[,1])
    edgepair <- c(row, column)
        
    g1_temp<-shiftrotate(g1,row)
    g2_temp<-shiftrotate(g2,column)
    tempscore<-rapidscore(g1_temp$g_aligned,g2_temp$g_aligned,tol)
     
    op<-list(score=tempscore, alignedEdgePair=edgepair, rank=rank, g1_shiftrot=g1_temp, g2_shiftrot=g2_temp)
    op
      
#     outputlist<-rbind(outputlist,op)
#     print(cbind(op$rank,op$score))
  
  #}
 
} # end of foreach loop
  
#   stopCluster(cl)
  
  
  #find lowest score here by converting output to dataframe, geting lowest score and corresponding features
  
  outputlist<-as.data.frame(outputlist)
  
  
  scorelist<-numeric()
  
  for(i in 1:dim(outputlist)[2]){
    
    scorelist<-c(scorelist,outputlist[[i]]$score)   
    
  }
  
  
  theta1_list<-numeric()
  for(i in 1:dim(outputlist)[2]){
    op<-outputlist[[i]]$g1_shiftrot$theta * (180)/pi 
    theta1_list<-c(theta1_list,op)
    }

  
  theta2_list<-numeric()
  for(i in 1:dim(outputlist)[2]){
    op<-outputlist[[i]]$g2_shiftrot$theta * (180)/pi
    theta2_list<-c(theta2_list,op)
    }

  orderIndex<-order(scorelist, decreasing=FALSE)
  thetaDiff<-theta1_list-theta2_list
  thetaSatisfy<-which(abs(thetaDiff) < maxtheta )
  
  orderIndexThetaSatisfy<-which(orderIndex %in% thetaSatisfy) #the ordered indices of orderIndex, that satisfy theta contraint
  filteredOrderIndex<-orderIndex[orderIndexThetaSatisfy] #index of edgepairs that satisfy theta contsraint, ordered in increasing order of score.
  lowestScore<-outputlist[[filteredOrderIndex[1]]]$score
  
  minIndices<-which(scorelist[filteredOrderIndex]==lowestScore) #if there are more than 1 best match score, get all of the associated edges, graphs and ranks
  
  for(i in minIndices){
    bestIndex<- rbind(bestIndex,filteredOrderIndex[i])
    bestMatchScore<-rbind(bestMatchScore,outputlist[[filteredOrderIndex[i]]]$score)
    bestAlignEdge<-rbind(bestAlignEdge,outputlist[[filteredOrderIndex[i]]]$alignedEdgePair)
    g1_shiftrot<-rbind(g1_shiftrot,outputlist[[filteredOrderIndex[i]]]$g1_shiftrot)
    g2_shiftrot<-rbind(g2_shiftrot,outputlist[[filteredOrderIndex[i]]]$g2_shiftrot)
    bestRank<-rbind(bestRank,outputlist[[filteredOrderIndex[i]]]$rank)    
    
  }
  
  
  output<-data.frame(graph1=g1_shiftrot, graph2=g2_shiftrot, bestEdgePair=bestAlignEdge, bestScore=bestMatchScore, bestRank=bestRank) 
  
  
    
  if(f_plot==1){
    for(i in 1:length(minIndices)){
      windows()
      plotreggraphs(output$graph1.g_aligned[[i]],output$graph2.g_aligned[[i]])
      #plotreggraphs_t(output$graph1.g_aligned[[i]],output$graph2.g_aligned[[i]]) #for SNIR
      titlestring<-paste("Alignment on rank", i," best edge pair (# ",bestIndex[i], ") d=", round(bestMatchScore,2) , sep="")
      title(titlestring)
      
    }#end of i loop
    
  }
  
  output
  
}# end of function
########################################################################################################

#testhere
#RG<-registerGraphs_maxtheta(g1,g2,tol, maxtheta, f_plot)



################################################################################################
registerGraphs_maxtheta_t<-function(g1,g2,tol,maxtheta, f_plot){
  
  g1_edgelength<-get.edge.value(g1,"edgelength")
  g2_edgelength<-get.edge.value(g2,"edgelength")
  g1_edgeslope<-get.edge.value(g1,"edgeslope")
  g2_edgeslope<-get.edge.value(g2,"edgeslope")
  g1_thickness<-get.edge.value(g1,"thickness")
  g2_thickness<-get.edge.value(g2,"thickness")
  
  g1_edgelist<-as.matrix.network(g1, matrix.type="edgelist")
  g2_edgelist<-as.matrix.network(g2, matrix.type="edgelist")
  
  
  #match edges by comparing lengths and slopes
  #===========================================
  
  edge_matching <- matrix(0, length(g1_edgelength), length(g2_edgelength) )
  
  for (i in 1:length(edge_matching[,1])){
    for (j in 1:length(edge_matching[1,])){
      #1
      edge_matching[i,j] <- sqrt((g1_edgelength[i]-g2_edgelength[j])**2)
      #2
      #edge_matching[i,j] <- sqrt((g1_edgelength[i]-g2_edgelength[j])**2+(10*(g1_thickness[i]-g2_thickness[j])**2))
      #3
      #edge_matching[i,j] <- sqrt((g1_thickness[i]-g2_thickness[j])**2)
      
      #edge_matching[i,j] <- edge_matching[i,j]/(0.5*g1_edgelength[i]+0.5*g2_edgelength[j])
      #edge_matching[i,j] <- edge_matching[i,j]/(d1+d2)
      
      
    }
  }
  
  ranking <- order(edge_matching, decreasing=FALSE)
  
  
  #take first 2500 closest edge pairs to realign and get best match score
  #================================================================
  NEdgePairs<-1000 #for SFIR, NEdgePairs=1000 for SNIR
  #NEdgePairs<-length(ranking)
  if(length(ranking)<NEdgePairs) { NEdgePairs<- length(ranking) }
  bestMatchScore<-numeric()
  bestAlignEdge<-numeric()
  #   g1_shiftrot<-g1
  #   g2_shiftrot<-g2
  g1_shiftrot<-numeric()
  g2_shiftrot<-numeric()
  scores<-numeric()
  alignedEdge<-numeric()
  bestRank<-numeric()
  bestIndex<-numeric()
  
  
  
  #   outputlist<-numeric()
  #   cl<-makeCluster(6, type = "SOCK")
  #   #clusterExport(cl, c("main.fun", "helper1", "helper2"))
  #   registerDoSNOW(cl)
  
  outputlist<-foreach(i=1:NEdgePairs, .combine="cbind", .packages=c("network","sna", "clue")  )%dopar%{
    source("H:/Uni Stuff/Uni Stuff/2012/PostDoc/Research Projects/ARC 2012-2015/R Code/Generic code AA/basicGraphFunctions.R")
    #for(i in 1:NEdgePairs){
    
    rank<-i
    
    column <- ceiling(ranking[rank]/length(edge_matching[,1]))
    row <- ranking[rank]%%length(edge_matching[,1])
    if (row==0) row <- length(edge_matching[,1])
    edgepair <- c(row, column)
    
    g1_temp<-shiftrotate(g1,row)
    g2_temp<-shiftrotate(g2,column)
    tempscore<-rapidscore(g1_temp$g_aligned,g2_temp$g_aligned,tol)
    
    op<-list(score=tempscore, alignedEdgePair=edgepair, rank=rank, g1_shiftrot=g1_temp, g2_shiftrot=g2_temp)
    op
    
    #     outputlist<-rbind(outputlist,op)
    #     print(cbind(op$rank,op$score))
    
    #}
    
  } # end of foreach loop
  
  #   stopCluster(cl)
  
  
  #find lowest score here by converting output to dataframe, geting lowest score and corresponding features
  
  outputlist<-as.data.frame(outputlist)
  
  
  scorelist<-numeric()
  
  for(i in 1:dim(outputlist)[2]){
    
    scorelist<-c(scorelist,outputlist[[i]]$score)   
    
  }
  
  
  theta1_list<-numeric()
  for(i in 1:dim(outputlist)[2]){
    op<-outputlist[[i]]$g1_shiftrot$theta * (180)/pi 
    theta1_list<-c(theta1_list,op)
  }
  
  
  theta2_list<-numeric()
  for(i in 1:dim(outputlist)[2]){
    op<-outputlist[[i]]$g2_shiftrot$theta * (180)/pi
    theta2_list<-c(theta2_list,op)
  }
  
  orderIndex<-order(scorelist, decreasing=FALSE)
  thetaDiff<-theta1_list-theta2_list
  thetaSatisfy<-which(abs(thetaDiff) < maxtheta )
  
  orderIndexThetaSatisfy<-which(orderIndex %in% thetaSatisfy) #the ordered indices of orderIndex, that satisfy theta contraint
  filteredOrderIndex<-orderIndex[orderIndexThetaSatisfy] #index of edgepairs that satisfy theta contsraint, ordered in increasing order of score.
  lowestScore<-outputlist[[filteredOrderIndex[1]]]$score
  
  minIndices<-which(scorelist[filteredOrderIndex]==lowestScore) #if there are more than 1 best match score, get all of the associated edges, graphs and ranks
  
  for(i in minIndices){
    bestIndex<- rbind(bestIndex,filteredOrderIndex[i])
    bestMatchScore<-rbind(bestMatchScore,outputlist[[filteredOrderIndex[i]]]$score)
    bestAlignEdge<-rbind(bestAlignEdge,outputlist[[filteredOrderIndex[i]]]$alignedEdgePair)
    g1_shiftrot<-rbind(g1_shiftrot,outputlist[[filteredOrderIndex[i]]]$g1_shiftrot)
    g2_shiftrot<-rbind(g2_shiftrot,outputlist[[filteredOrderIndex[i]]]$g2_shiftrot)
    bestRank<-rbind(bestRank,outputlist[[filteredOrderIndex[i]]]$rank)    
    
  }
  
  
  output<-data.frame(graph1=g1_shiftrot, graph2=g2_shiftrot, bestEdgePair=bestAlignEdge, bestScore=bestMatchScore, bestRank=bestRank) 
  
  
  
  if(f_plot==1){
    for(i in 1:length(minIndices)){
      windows()
      #plotreggraphs(output$graph1.g_aligned[[i]],output$graph2.g_aligned[[i]])
      plotreggraphs_t(output$graph1.g_aligned[[i]],output$graph2.g_aligned[[i]]) #for SNIR
      titlestring<-paste("Alignment with thickness on rank", i," best edge pair (# ",bestIndex[i], ")", sep="")
      title(titlestring)
      
    }#end of i loop
    
  }
  
  output
  
}# end of function
########################################################################################################


##############################################################################################

#########

registerGraphs_hungarian<-function(g1,g2,cid,f_plot,f_edgeinclude, maxtheta){
  
  g1_edgelength<-get.edge.value(g1,"edgelength")
  g2_edgelength<-get.edge.value(g2,"edgelength")
  g1_edgeslope<-get.edge.value(g1,"edgeslope")
  g2_edgeslope<-get.edge.value(g2,"edgeslope")
  
  
  #match edges by comparing lengths and slopes
  #===========================================
  
  edge_matching <- matrix(0, length(g1_edgelength), length(g2_edgelength) )
  
  for (i in 1:length(edge_matching[,1])){
    for (j in 1:length(edge_matching[1,])){
      
      edge_matching[i,j] <- sqrt((g1_edgelength[i]-g2_edgelength[j])**2+(g1_edgeslope[i]-g2_edgeslope[j])**2)
      edge_matching[i,j] <- edge_matching[i,j]/(0.5*g1_edgelength[i]+0.5*g2_edgelength[j])
      
    }
  }
  
  ranking <- order(edge_matching, decreasing=FALSE)
  
  #take first 50 closest edge pairs to realign and get best match score
  #================================================================
  NEdgePairs<-1000
  #NEdgePairs<-length(ranking)
  if(length(ranking)<NEdgePairs) { NEdgePairs<- length(ranking) }
  bestIndex<-numeric()
  bestMatchScore<-numeric()
  bestAlignEdge<-numeric()
  g1_shiftrot<-numeric()
  g2_shiftrot<-numeric()
  bestRank<-numeric()
  bestMCS<-numeric()
  
  
  
  
  outputlist<-numeric()
  
  outputlist<-foreach(i=1:NEdgePairs, .combine="cbind", .packages=c("network","sna", "clue") )%dopar%{
    
    source("H:/Uni Stuff/Uni Stuff/2012/PostDoc/Research Projects/ARC 2012-2015/R Code/Generic code AA/basicGraphFunctions.R")
    source("H:/Uni Stuff/Uni Stuff/2012/PostDoc/Research Projects/ARC 2012-2015/R Code/Generic code AA/graphMatchingFunctions.R")
    
    rank<-i
    
    column <- ceiling(ranking[rank]/length(edge_matching[,1]))
    row <- ranking[rank]%%length(edge_matching[,1])
    if (row==0) row <- length(edge_matching[,1])
    edgepair <- c(row, column)
    
    
    g1_temp<-shiftrotate(g1,row)
    g2_temp<-shiftrotate(g2,column)
    
    Path<-getGraphEditPath(g1_temp$g_aligned,g2_temp$g_aligned,cid,f_edgeinclude)
    MCS<-getMCS(g1_temp$g_aligned,g2_temp$g_aligned,Path$editPath)
    tempscore<-getDist(g1_temp$g_aligned,g2_temp$g_aligned,MCS,cid)$dbar_sqrt_n
    
    op<-list(score=tempscore, alignedEdgePair=edgepair, rank=rank, g1_shiftrot=g1_temp, g2_shiftrot=g2_temp, mcs=MCS)
    #op
    #outputlist<-rbind(outputlist,op)
    
  }
  
  #end parallal code
  #--
  
  #find lowest score here by converting output to dataframe, geting lowest score and corresponding features
  
  outputlist<-as.data.frame(outputlist)
  
  scorelist<-numeric()
 
 for(i in 1:dim(outputlist)[2]){
    
    scorelist<-c(scorelist,outputlist[[i]]$score)   
    
  }
  
  
  theta1_list<-numeric()
  for(i in 1:dim(outputlist)[2]){
    op<-outputlist[[i]]$g1_shiftrot$theta * (180)/pi 
    theta1_list<-c(theta1_list,op)
  }
  
  
  theta2_list<-numeric()
  for(i in 1:dim(outputlist)[2]){
    op<-outputlist[[i]]$g2_shiftrot$theta * (180)/pi
    theta2_list<-c(theta2_list,op)
  }
  
  orderIndex<-order(scorelist, decreasing=FALSE)
  thetaDiff<-theta1_list-theta2_list
  thetaSatisfy<-which(abs(thetaDiff) < maxtheta )
  
  orderIndexThetaSatisfy<-which(orderIndex %in% thetaSatisfy) #the ordered indices of orderIndex, that satisfy theta contraint
  filteredOrderIndex<-orderIndex[orderIndexThetaSatisfy] #index of edgepairs that satisfy theta contsraint, ordered in increasing order of score.
  lowestScore<-outputlist[[filteredOrderIndex[1]]]$score
  
  minIndices<-which(scorelist[filteredOrderIndex]==lowestScore) #if there are more than 1 best match score, get all of the associated edges, graphs and ranks
  
  for(i in minIndices){
    bestIndex<- rbind(bestIndex,filteredOrderIndex[i])
    bestMatchScore<-rbind(bestMatchScore,outputlist[[filteredOrderIndex[i]]]$score)
    bestAlignEdge<-rbind(bestAlignEdge,outputlist[[filteredOrderIndex[i]]]$alignedEdgePair)
    g1_shiftrot<-rbind(g1_shiftrot,outputlist[[filteredOrderIndex[i]]]$g1_shiftrot)
    g2_shiftrot<-rbind(g2_shiftrot,outputlist[[filteredOrderIndex[i]]]$g2_shiftrot)
    bestRank<-rbind(bestRank,outputlist[[filteredOrderIndex[i]]]$rank)   
    bestMCS<-rbind(bestMCS,outputlist[[filteredOrderIndex[i]]]$mcs)
    
  }
  
  
  output<-data.frame(graph1=g1_shiftrot,graph2=g2_shiftrot,bestEdgePair=bestAlignEdge,bestScore=bestMatchScore,bestRank=bestRank, bestMCS=bestMCS)
  
  
  
  if(f_plot==1){
    for(i in 1:length(minIndices)){
      windows()
      plotreggraphs(output$graph1.g_aligned[[i]],output$graph2.g_aligned[[i]])
      titlestring<-paste("Alignment on rank", i," best edge pair (# ",bestIndex[i], ")", sep="")
      title(titlestring)
      
    }#end of i loop
    
  }
  
 output 
  
}#end of function

#########################################################################################################





##########################################################################################################


#############START registerGraphs_maxtheta_lbound  , maxtheta is the maximum rotation allowed for a pair of matching edges, lbound bounds the minimum length that matching edges should have
#######################################################
registerGraphs_maxtheta_lbound<-function(g1,g2,tol,maxtheta, lbound, f_plot){
  
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
      
      if(g1_edgelength[i]<lbound || g2_edgelength[j]< lbound) {edge_matching[i,j]<-10000}
      if( sqrt((g1_edgelength[i]-g2_edgelength[j])**2)> (2*tol) ) {edge_matching[i,j]<-1000}
      
      
    }
  }
  
  ranking <- order(edge_matching, decreasing=FALSE)
  
  
  #take first 2500 closest edge pairs to realign and get best match score
  #================================================================
  NEdgePairs<-1000
  #NEdgePairs<-length(ranking)
  if(length(ranking)<NEdgePairs) { NEdgePairs<- length(ranking) }
  bestMatchScore<-numeric()
  bestAlignEdge<-numeric()
  #   g1_shiftrot<-g1
  #   g2_shiftrot<-g2
  g1_shiftrot<-numeric()
  g2_shiftrot<-numeric()
  scores<-numeric()
  alignedEdge<-numeric()
  bestRank<-numeric()
  bestIndex<-numeric()
  
  
  
  #   outputlist<-numeric()
  #   cl<-makeCluster(6, type = "SOCK")
  #   #clusterExport(cl, c("main.fun", "helper1", "helper2"))
  #   registerDoSNOW(cl)
  
  outputlist<-foreach(i=1:NEdgePairs, .combine="cbind", .packages=c("network","sna", "clue")  )%dopar%{
    source("H:/Uni Stuff/Uni Stuff/2012/PostDoc/Research Projects/ARC 2012-2015/R Code/Generic code AA/basicGraphFunctions.R")
    #for(i in 1:NEdgePairs){
    
    rank<-i
    
    column <- ceiling(ranking[rank]/length(edge_matching[,1]))
    row <- ranking[rank]%%length(edge_matching[,1])
    if (row==0) row <- length(edge_matching[,1])
    edgepair <- c(row, column)
    
    g1_temp<-shiftrotate(g1,row)
    g2_temp<-shiftrotate(g2,column)
    tempscore<-rapidscore(g1_temp$g_aligned,g2_temp$g_aligned,tol)
    
    op<-list(score=tempscore, alignedEdgePair=edgepair, rank=rank, g1_shiftrot=g1_temp, g2_shiftrot=g2_temp)
    op
    
    #     outputlist<-rbind(outputlist,op)
    #     print(cbind(op$rank,op$score))
    
    #}
    
  } # end of foreach loop
  
  #   stopCluster(cl)
  
  
  #find lowest score here by converting output to dataframe, geting lowest score and corresponding features
  
  outputlist<-as.data.frame(outputlist)
  
  
  scorelist<-numeric()
  
  for(i in 1:dim(outputlist)[2]){
    
    scorelist<-c(scorelist,outputlist[[i]]$score)   
    
  }
  
  
  theta1_list<-numeric()
  for(i in 1:dim(outputlist)[2]){
    op<-outputlist[[i]]$g1_shiftrot$theta * (180)/pi 
    theta1_list<-c(theta1_list,op)
  }
  
  
  theta2_list<-numeric()
  for(i in 1:dim(outputlist)[2]){
    op<-outputlist[[i]]$g2_shiftrot$theta * (180)/pi
    theta2_list<-c(theta2_list,op)
  }
  
  orderIndex<-order(scorelist, decreasing=FALSE)
  thetaDiff<-theta1_list-theta2_list
  thetaSatisfy<-which(abs(thetaDiff) < maxtheta )
  
  orderIndexThetaSatisfy<-which(orderIndex %in% thetaSatisfy) #the ordered indices of orderIndex, that satisfy theta contraint
  filteredOrderIndex<-orderIndex[orderIndexThetaSatisfy] #index of edgepairs that satisfy theta contsraint, ordered in increasing order of score.
  lowestScore<-outputlist[[filteredOrderIndex[1]]]$score
  
  minIndices<-which(scorelist[filteredOrderIndex]==lowestScore) #if there are more than 1 best match score, get all of the associated edges, graphs and ranks
  
  for(i in minIndices){
    bestIndex<- rbind(bestIndex,filteredOrderIndex[i])
    bestMatchScore<-rbind(bestMatchScore,outputlist[[filteredOrderIndex[i]]]$score)
    bestAlignEdge<-rbind(bestAlignEdge,outputlist[[filteredOrderIndex[i]]]$alignedEdgePair)
    g1_shiftrot<-rbind(g1_shiftrot,outputlist[[filteredOrderIndex[i]]]$g1_shiftrot)
    g2_shiftrot<-rbind(g2_shiftrot,outputlist[[filteredOrderIndex[i]]]$g2_shiftrot)
    bestRank<-rbind(bestRank,outputlist[[filteredOrderIndex[i]]]$rank)    
    
  }
  
  
  output<-data.frame(graph1=g1_shiftrot, graph2=g2_shiftrot, bestEdgePair=bestAlignEdge, bestScore=bestMatchScore, bestRank=bestRank) 
  
  
  
  if(f_plot==1){
    for(i in 1:length(minIndices)){
      windows()
      plotreggraphs(output$graph1.g_aligned[[i]],output$graph2.g_aligned[[i]])
      titlestring<-paste("Alignment on rank", i," best edge pair (# ",bestIndex[i], ")", sep="")
      title(titlestring)
      
    }#end of i loop
    
  }
  
  output
  
}# end of function
########################################################################################################

#############START registerGraphs_maxtheta  , maxtheta is the maximum rotation allowed for a pair of matching edges
#######################################################
registerGraphs_new<-function(g1,g2,tol,maxtheta, f_plot){
  
  g1_edgelength<-get.edge.value(g1,"edgelength")
  g2_edgelength<-get.edge.value(g2,"edgelength")
  g1_edgeslope<-get.edge.value(g1,"edgeslope")
  g2_edgeslope<-get.edge.value(g2,"edgeslope")
  
  g1_edgelist<-as.matrix.network(g1, matrix.type="edgelist")
  g2_edgelist<-as.matrix.network(g2, matrix.type="edgelist")
  
  g1_edgesorted<-sort(g1_edgelength, decreasing=TRUE )
  g1_edgeorder<-order(g1_edgelength, decreasing=TRUE )
  g2_edgesorted<-sort(g2_edgelength, decreasing=TRUE)
  g2_edgeorder<-order(g2_edgelength, decreasing=TRUE )
  
  # take the longest 50% of edges
  s1<-length(g1_edgesorted)/2
  s2<-length(g2_edgesorted)/2
  g1_edgesorted<-g1_edgesorted[1:s1]
  g1_edgeorder<-g1_edgeorder[1:s1]
  g2_edgesorted<- g2_edgesorted[1:s2]
  g2_edgeorder<-g2_edgeorder[1:s2]
  
  #match edges by comparing lengths and slopes
  #===========================================
  
  edge_matching <- matrix(0, length(g1_edgesorted), length( g2_edgesorted) )
  
  for (i in 1:length(g1_edgesorted)){
    for (j in 1:length(g2_edgesorted )){
      
      #      d1<-length(get.neighborhood(g1,g1_edgelist[i,1], type="combined"))+length(get.neighborhood(g1,g1_edgelist[i,2], type="combined"))
      #      d2<-length(get.neighborhood(g2,g2_edgelist[j,1], type="combined"))+length(get.neighborhood(g2,g2_edgelist[j,2], type="combined"))
      #      
      
      edge_matching[i,j] <- abs(g1_edgesorted[i]-  g2_edgesorted[j])
      
     # edge_matching[i,j] <- sqrt((g1_edgelength[i]-g2_edgelength[j])**2+(g1_edgeslope[i]-g2_edgeslope[j])**2)
      #edge_matching[i,j] <- edge_matching[i,j]/(0.5*g1_edgelength[i]+0.5*g2_edgelength[j])
      #edge_matching[i,j] <- edge_matching[i,j]/(d1+d2)
      
      
    }
  }
  
  ranking <- order(edge_matching, decreasing=FALSE)
  
  
  #take first 2500 closest edge pairs to realign and get best match score
  #================================================================
  NEdgePairs<-84 #for SFIR, NEdgePairs=1000 for SNIR
  #NEdgePairs<-length(ranking)
  if(length(ranking)<NEdgePairs) { NEdgePairs<- length(ranking) }
  bestMatchScore<-numeric()
  bestAlignEdge<-numeric()
  #   g1_shiftrot<-g1
  #   g2_shiftrot<-g2
  g1_shiftrot<-numeric()
  g2_shiftrot<-numeric()
  scores<-numeric()
  alignedEdge<-numeric()
  bestRank<-numeric()
  bestIndex<-numeric()
  
  
  
  #   outputlist<-numeric()
  #   cl<-makeCluster(6, type = "SOCK")
  #   #clusterExport(cl, c("main.fun", "helper1", "helper2"))
  #   registerDoSNOW(cl)
  
  outputlist<-foreach(i=1:NEdgePairs, .combine="cbind", .packages=c("network","sna", "clue")  )%dopar%{
    source("H:/Uni Stuff/Uni Stuff/2012/PostDoc/Research Projects/ARC 2012-2015/R Code/Generic code AA/basicGraphFunctions.R")
    #for(i in 1:NEdgePairs){
    
    rank<-i
    
    column <- ceiling(ranking[rank]/length(edge_matching[,1]))
    row <- ranking[rank]%%length(edge_matching[,1])
    if (row==0) row <- length(edge_matching[,1])
    row<-g1_edgeorder[row]
    column<-g2_edgeorder[column]
    edgepair <- c(row, column)
    
    g1_temp<-shiftrotate(g1,row)
    g2_temp<-shiftrotate(g2,column)
    tempscore<-rapidscore(g1_temp$g_aligned,g2_temp$g_aligned,tol)
    
    op<-list(score=tempscore, alignedEdgePair=edgepair, rank=rank, g1_shiftrot=g1_temp, g2_shiftrot=g2_temp)
    op
    
    #     outputlist<-rbind(outputlist,op)
    #     print(cbind(op$rank,op$score))
    
    #}
    
  } # end of foreach loop
  
  #   stopCluster(cl)
  
  
  #find lowest score here by converting output to dataframe, geting lowest score and corresponding features
  
  outputlist<-as.data.frame(outputlist)
  
  
  scorelist<-numeric()
  
  for(i in 1:dim(outputlist)[2]){
    
    scorelist<-c(scorelist,outputlist[[i]]$score)   
    
  }
  
  
  theta1_list<-numeric()
  for(i in 1:dim(outputlist)[2]){
    op<-outputlist[[i]]$g1_shiftrot$theta * (180)/pi 
    theta1_list<-c(theta1_list,op)
  }
  
  
  theta2_list<-numeric()
  for(i in 1:dim(outputlist)[2]){
    op<-outputlist[[i]]$g2_shiftrot$theta * (180)/pi
    theta2_list<-c(theta2_list,op)
  }
  
  orderIndex<-order(scorelist, decreasing=FALSE)
  thetaDiff<-theta1_list-theta2_list
  thetaSatisfy<-which(abs(thetaDiff) < maxtheta )
  
  orderIndexThetaSatisfy<-which(orderIndex %in% thetaSatisfy) #the ordered indices of orderIndex, that satisfy theta contraint
  filteredOrderIndex<-orderIndex[orderIndexThetaSatisfy] #index of edgepairs that satisfy theta contsraint, ordered in increasing order of score.
  lowestScore<-outputlist[[filteredOrderIndex[1]]]$score
  
  minIndices<-which(scorelist[filteredOrderIndex]==lowestScore) #if there are more than 1 best match score, get all of the associated edges, graphs and ranks
  
  for(i in minIndices){
    bestIndex<- rbind(bestIndex,filteredOrderIndex[i])
    bestMatchScore<-rbind(bestMatchScore,outputlist[[filteredOrderIndex[i]]]$score)
    bestAlignEdge<-rbind(bestAlignEdge,outputlist[[filteredOrderIndex[i]]]$alignedEdgePair)
    g1_shiftrot<-rbind(g1_shiftrot,outputlist[[filteredOrderIndex[i]]]$g1_shiftrot)
    g2_shiftrot<-rbind(g2_shiftrot,outputlist[[filteredOrderIndex[i]]]$g2_shiftrot)
    bestRank<-rbind(bestRank,outputlist[[filteredOrderIndex[i]]]$rank)    
    
  }
  
  
  output<-data.frame(graph1=g1_shiftrot, graph2=g2_shiftrot, bestEdgePair=bestAlignEdge, bestScore=bestMatchScore, bestRank=bestRank) 
  
  
  
  if(f_plot==1){
    for(i in 1:length(minIndices)){
      windows()
      plotreggraphs(output$graph1.g_aligned[[i]],output$graph2.g_aligned[[i]])
      #plotreggraphs_t(output$graph1.g_aligned[[i]],output$graph2.g_aligned[[i]]) #for SNIR
      titlestring<-paste("Alignment on rank", i," best edge pair (# ",bestIndex[i], ")", sep="")
      title(titlestring)
      
    }#end of i loop
    
  }
  
  output
  
}# end of function
########################################################################################################

#testhere
#RG<-registerGraphs_new(g1,g2,tol, maxtheta, 1)

getnodes <- function() {
  f <- Sys.getenv('PBS_NODEFILE')
  x <- if (nzchar(f)) readLines(f) else rep('localhost', 3)
  as.data.frame(table(x), stringsAsFactors=FALSE)
}

setcores <- function(cl, nodes) {
  f <- function(cores) assign('allocated.cores', cores, pos=.GlobalEnv)
  clusterApply(cl, nodes$Freq, f)
}
