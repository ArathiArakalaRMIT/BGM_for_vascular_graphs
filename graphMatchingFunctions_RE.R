library(splancs)
library(spatial)
library(network)
library(clue)
source("basicGraphFunctions_vpac.R")

getGraph_RE<-function(wdname,sourceid,person){
  g<-network.initialize(1)
  
  inputNodesFile<-""
  inputEdgesFile<-""
  if(person<10) personid<-paste("00", person, sep="")
  if(person>=10 && person<=99 ) personid<-paste("0", person, sep="")
  if(person >=100) personid<-person
  
  inputNodesFile<-paste(wdname,"/R",personid,"_bpoints.txt",sep="")
  inputEdgesFile<-paste(wdname,"/R",personid,"_links.txt", sep="")
  
  #checkiffileexists
  if(file.exists(inputNodesFile) && file.exists(inputEdgesFile)){
    #Create graph
    #-------------
    v <- read.table(inputNodesFile,header=FALSE, sep="\t", fill=TRUE)
    N<-length(v[,1])
    e<-read.table(inputEdgesFile,header=FALSE, sep="\t", fill=TRUE)
    e<-as.matrix(e)
    
    #arrange edges with lower vertex index first and higher vertex index second
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
    set.vertex.attribute(g,"type", v[,4]) 
    
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
        
  }
  
  
  #return graph
  #-------------
  g
  
}#end of getGraph function
###########################################################################
