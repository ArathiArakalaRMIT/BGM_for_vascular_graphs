# This code will generate the topology distributions of 4 vascular biometric graphs - PV, WV, HV and Retina.
wd<-"H:/Uni Stuff/Uni Stuff/2015/Research Projects/ARC Discovery/R Code/Generic Code AA/BGM_Vascular/"
setwd(wd)

library(splancs)
library(spatial)
library(network)
library(clue)
library(iterators)
library(foreach)
library(doSNOW)
library(sna)
library(gtools)
require(graphics)
require(utils)
#library(doMC)

# include files

source("basicGraphFunctions_vpac.R")
source("graphRegFunctions.R")
source("graphMatchingFunctions_vpac.R")
source("graphMatchingFunctions_RE.R")
source("graphMatchingFunctions_PV.R")
source("graphMatchingFunctions_WV.R")
source("graphMatchingFunctions_HV.R")
source("graphMatchingFunctions_parallal_vpac.R")
source("graphTopologicalFeatures_vpac.R")
source("graphOperations_vpac.R")


biometric<-"RE" 
sourceid <- "AutoHao"
graphList<-numeric()
outputfile<-paste("H:/Uni Stuff/Uni Stuff/2015/Research Projects/ARC Discovery/Outputs/RetinaGraphs/", sourceid, "/graphList_",sourceid,".Rdata", sep="")
load(outputfile)
graphList_RE<-graphList

biometric<-"PV" 
sourceid <- "AutoHao2L"
graphList<-numeric()
wd<-"H:/Uni Stuff/Uni Stuff/2015/Research Projects/ARC Discovery/Outputs/PalmVeinGraphs"
outputfile<-paste(wd, "/Method2Output/graphList_",sourceid,".Rdata", sep="")
load(outputfile)
graphList_PV_L<-graphList

biometric<-"PV" 
sourceid <- "AutoHao2R"
graphList<-numeric()
wd<-"H:/Uni Stuff/Uni Stuff/2015/Research Projects/ARC Discovery/Outputs/PalmVeinGraphs"
outputfile<-paste(wd, "/Method2Output/graphList_",sourceid,".Rdata", sep="")
load(outputfile)
graphList_PV_R<-graphList

biometric<-"WV" 
sourceid <- "AutoHao7L_v1"
graphList<-numeric()
wd<-"H:/Uni Stuff/Uni Stuff/2015/Research Projects/ARC Discovery/Outputs/WristVeinGraphs"
outputfile<-paste(wd, "/Method7Output/graphList_",sourceid,".Rdata", sep="")
load(outputfile)
graphList_WV_L<-graphList

biometric<-"WV" 
sourceid <- "AutoHao7R_v1"
graphList<-numeric()
wd<-"H:/Uni Stuff/Uni Stuff/2015/Research Projects/ARC Discovery/Outputs/WristVeinGraphs"
outputfile<-paste(wd, "/Method7Output/graphList_",sourceid,".Rdata", sep="")
load(outputfile)
graphList_WV_R<-graphList

biometric<-"HV" 
sourceid <- c("AutoSL/SFIR_LR", "Right")
outputfile<-paste("H:/Uni Stuff/Uni Stuff/2015/Research Projects/ARC Discovery/Outputs/HandVeinGraphs/",sourceid[1],"/graphList_",sourceid[2],".Rdata", sep="")
load(outputfile)
graphList_HV_SFIR<-graphList

biometric<-"HV" 
sourceid<-c("AutoSL/SNIR_T_LR", "Right")
outputfile<-paste("H:/Uni Stuff/Uni Stuff/2015/Research Projects/ARC Discovery/Outputs/HandVeinGraphs/",sourceid[1],"/graphList_",sourceid[2],".Rdata", sep="")
load(outputfile)
graphList_HV_SNIR<-graphList


biometric_list<-c("RE", "PV", "PV", "WV", "WV", "HV","HV")
sourceid_list<-c("AutoHao", "AutoHao2L", "AutoHao2R", "AutoHao7L_v1", "AutoHao7R_v1", "SFIR_LR", "SNIR_T_LR")
wd<-"H:/Uni Stuff/Uni Stuff/2015/Research Projects/ARC Discovery/Outputs/VascularGraphs/Topology/"

for(b in 1:7){
  biometric<-biometric_list[b]
  sourceid<-sourceid_list[b]
  if(b %in% 6:7){
    sourceid<-c(sourceid, "Right")
  }
  
  if(b==1) graphList<-graphList_RE
  if(b==2) graphList<-graphList_PV_L
  if(b==3) graphList<-graphList_PV_R
  if(b==4) graphList<-graphList_WV_L
  if(b==5) graphList<-graphList_WV_R
  if(b==6) graphList<-graphList_HV_SFIR
  if(b==7) graphList<-graphList_HV_SNIR
  
  Size<-rep(0, times=dim(graphList)[1])
  Edgecount<-rep(0, times=dim(graphList)[1])
  Edgenoderatio<-rep(0, times=dim(graphList)[1])
  C1<-rep(0, times=dim(graphList)[1])
  C2<-rep(0, times=dim(graphList)[1])
  Isolated<-rep(0, times=dim(graphList)[1])
  C1_l<-rep(0, times=dim(graphList)[1])
  C1c2<-rep(0, times=dim(graphList)[1])
  ENRc1<-rep(0, times=dim(graphList)[1]) # number of edges in c1 to size of mcs
  Dmax<-rep(0, times=dim(graphList)[1])
  starSize<-rep(0, times=dim(graphList)[1])
  twostarSize<-rep(0, times=dim(graphList)[1])
  degreedisbn<-matrix(0, nrow=dim(graphList)[1], ncol=150)
  
  
  
  for(i in 1:dim(graphList)[1]){
    
    g<-graphList$graph[[i]]
    
    
    Size[i]<-network.size(g)
    Edgecount[i]<-network.edgecount(g)
    Edgenoderatio[i]<-network.edgecount(g)/network.size(g)
    cd<-component.dist(g)
    largest_comp<-which(cd$csize==max(cd$csize))[1]
    
    C1[i]<-sort(cd$csize, decreasing=TRUE)[1]
    C2[i]<-sort(cd$csize, decreasing=TRUE)[2]
    Isolated[i]<-cd$cdist[1]
    
    n_c1<-which( cd$membership==largest_comp)
    n_todelete<-setdiff(1:network.size(g), n_c1)
    g_c1<-g
    delete.vertices(g_c1, n_todelete)
    edges_g_c1<-0
    if(network.edgecount(g_c1)>0) C1_l[i]<-sum( get.edge.value(g_c1,"edgelength") )
    C1c2[i]<-(sort(cd$csize, decreasing=TRUE)[1]+sort(cd$csize, decreasing=TRUE)[2])
    ENRc1[i]<-network.edgecount(g_c1)/network.size(g_c1)
    
    deg<-degree(g, gmode="graph", cmode="indegree")
    Dmax[i]<-max(deg)
    
    starSize[i]<-dim(get_star(g))[1]
    twostarSize[i]<-dim(get_twostar(g))[1]  
    
    degree_tmp<-rep(0, times=max(deg) )
    for(j in 1:max(deg)){
      degree_tmp[j]<-length(which(deg==j))
    }
    degree_tmp<-degree_tmp/sum(degree_tmp)
    degreedisbn[i,1:max(deg)]<-degree_tmp
    
    print(i) 
  }#end of i loop
  
  
  lastIndex<-max(Dmax)
  torem<-((lastIndex)+1):dim(degreedisbn)[2]
  degreedisbn<-degreedisbn[,-torem]
  degree_mean<-rep(0, times=lastIndex)
  degree_sd<-rep(0, times=lastIndex)
  for(j in 1:lastIndex){
    degree_mean[j]<-mean(degreedisbn[,j])
    degree_sd[j]<-sd(degreedisbn[,j])
  }
  
  #modify code and plot to file
  pdf(paste(wd, "DisbnPlots_", biometric,"_",sourceid[1],".pdf", sep=""),width=7, height=7, useDingbats=FALSE)
  #for(c in 1:12){  
  for(c in c(1,2,3,4, 10,11,12)){
    data<-numeric()
    if(c==1){ data<-Size
              titleString<-"Size"} 
    if(c==2) { data<-Edgecount
               titleString<-"Edgecount"}
    if(c==3) {data<-Edgenoderatio
              titleString<-"Edgenoderatio"} 
    if(c==4){ data<-C1
              titleString<-"C1"}
    
    if(c==5){data<-C2
             titleString<-"C2"}
    
    if(c==6){data<-Isolated
             titleString<-"Isolated"}
    
    if(c==7){data<-C1_l
             titleString<-"C1_l"}
    
    if(c==8) { data<-C1c2
               titleString<-"C1c2"}
    
    if(c==9) { data<-ENRc1
               titleString<-"ENRc1"}
    
    if(c==10){data<-Dmax
              titleString<-"Dmax"}
    
    if(c==11) {data<-starSize
               titleString<-"starSize"}
    
    if(c==12){data<-twostarSize
              titleString<-"twostarSize"}
    
    
#     h<-hist(data, breaks=10, plot=FALSE)
#     barplot(h$counts/sum(h$counts), col="grey", ylim=c(0, 1), xlim=c(0, length(h$breaks)), names.arg=h$mids, main=titleString, xlab="Data Values", ylab="Relative Frequency " )
#     #write all results to file with proper labeling of database
    outputfile<-paste(wd,titleString,"_", biometric,"_",sourceid[1],".Rdata", sep="")
    #save(data, file=outputfile)
    load(outputfile) 
    print(mean(data))
  }
  30plot(1:length(degree_mean), degree_mean, main="Degree Distribution", type='b', xlim=c(0,length(degree_mean) ), ylim=c(0,1), xlab="Degree", ylab="Probability", pch=2, cex=2, lwd=2)
  dev.off()
  outputfile<-paste(wd,"DegreeDisbn","_", biometric,"_",sourceid[1],".Rdata", sep="")
  save(degreedisbn, file=outputfile)
  load(outputfile)
  outputfile<-paste(wd,"DegreeDisbn_mean","_", biometric,"_",sourceid[1],".Rdata", sep="")
  save(degree_mean, file=outputfile)
  load(outputfile)
  outputfile<-paste(wd,"DegreeDisbn_sd","_", biometric,"_",sourceid[1],".Rdata", sep="")
  save(degree_sd, file=outputfile)
  load(outputfile)    
}#end of b loop
