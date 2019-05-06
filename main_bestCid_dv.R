# This code will determine the best cid parameters for graph matching for each of the databases, using the Bunke Shearer Metric, dv.
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

biometric_list<-c("RE", "PV", "PV", "WV", "WV", "HV","HV")
sourceid_list<-c("AutoHao", "AutoHao2L", "AutoHao2R", "AutoHao7L_v1", "AutoHao7R_v1", "SFIR_LR", "SNIR_T_LR")
folder_list<-c("RetinaGraphs","PalmVeinGraphs", "PalmVeinGraphs", "WristVeinGraphs", "WristVeinGraphs", "HandVeinGraphs", "HandVeinGraphs")
opwd<-"H:/Uni Stuff/Uni Stuff/2015/Research Projects/ARC Discovery/Outputs/VascularGraphs/bestCid/"
ipwd<-"H:/Uni Stuff/Uni Stuff/2015/Research Projects/ARC Discovery/Outputs"
bestStruct<-c("twostar", "edge", "edge", "edge", "edge", "twostar", "edge")
mcsStruct<-c("node", "edge", "star", "twostar")
for(b in c(1,2,4,6,7)){
  print("#######################################")
  print(b)
  
  #retrieve the training comparison indices and the graphs. 
  graphList<-numeric()
  outputfile<-paste(ipwd,"/",folder_list[b],"/", sourceid_list[b], "/graphList_",sourceid_list[b],".Rdata", sep="")
  #save(graphList, file=outputfile)
  load(outputfile)
  
  
  outputfile<-paste(ipwd,"/",folder_list[b],"/",sourceid_list[b],"/RegTests/trainingPersonList.Rdata", sep="")
  #save(trainingPersonList, file=outputfile)
  load(outputfile)
  
  outputfile<-paste(ipwd,"/",folder_list[b],"/",sourceid_list[b],"/RegTests/testingPersonList.Rdata", sep="")
  #save(testingPersonList, file=outputfile)
  load(outputfile)
  
  outputfile<-paste(ipwd,"/",folder_list[b],"/",sourceid_list[b],"/RegTests/gencomp_train.Rdata",sep="")
  #save(gencomp_train, file=outputfile)
  load(outputfile)
  outputfile<-paste(ipwd,"/",folder_list[b],"/",sourceid_list[b],"/RegTests/impcomp_train.Rdata", sep="")
  #save(impcomp_train, file=outputfile)
  load(outputfile)
  outputfile<-paste(ipwd,"/",folder_list[b],"/",sourceid_list[b],"/RegTests/gencomp_test.Rdata", sep="")
  #save(gencomp_test, file=outputfile)
  load(outputfile)
  outputfile<-paste(ipwd,"/",folder_list[b],"/",sourceid_list[b],"/RegTests/impcomp_test.Rdata", sep="")
  #save(impcomp_test, file=outputfile)
  load(outputfile)
  
  genComp<-gencomp_train
  impComp<-impcomp_train
  cid_n_range<-seq(from=3, to=9, by=2) # see lncs2015 paper, tuning done there
  cid_e_range<-seq(from=3, to=9, by=2)
  
  
for(m in 1:length(mcsStruct)){
  
  dv_matrix<-matrix(100, nrow=length(cid_n_range), ncol=length(cid_e_range))
  de_matrix<-matrix(100, nrow=length(cid_n_range), ncol=length(cid_e_range))
  
  
  for(cid_n in cid_n_range){
    for(cid_e in cid_e_range){
      cid_n_index<-which(cid_n_range==cid_n)
      cid_e_index<-which(cid_e_range==cid_e)
      
     
        
        #read all the files and their cid_n and cid_e
        outputfile<-paste(ipwd,"/",folder_list[b],"/",sourceid_list[b],"/BetterGraphMatching/genMCS_",mcsStruct[m],"Based_train_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tests on training data
        if(file.exists(outputfile)==FALSE) next
        load(outputfile)
        if(m==1) genMCS<-genMCS_nodeBased
        if(m==2) genMCS<-genMCS_edgeBased
        if(m==3) genMCS<-genMCS_starBased
        if(m==4) genMCS<-genMCS_twostarBased
        
      
        
        outputfile<-paste(ipwd,"/",folder_list[b],"/",sourceid_list[b],"/BetterGraphMatching/impMCS_",mcsStruct[m],"Based_train_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tests on training data
        if(file.exists(outputfile)==FALSE) next
        load(outputfile)
        if(m==1) impMCS<-impMCS_nodeBased
        if(m==2) impMCS<-impMCS_edgeBased
        if(m==3) impMCS<-impMCS_starBased
        if(m==4) impMCS<-impMCS_twostarBased
        
        
        
        
        gen_dv<-rep(0, times=dim(genMCS)[1])
        gen_de<-rep(0, times=dim(genMCS)[1])
        
        for(i in 1:dim(genMCS)[1]){
          
          mcs<-genMCS[[i]]
          
          i1<-genComp[i,1]
          i2<-genComp[i,2]
          g1<-graphList$graph[[i1]]
          g2<-graphList$graph[[i2]]
          
          #         RG3<-numeric()
          #         t3<-system.time(RG3<-registerGraphs_twostars(g1,g2,tol,maxtheta, f_plot, L))
          #         
          #         
          #         windows()
          #         par(mfrow=c(2,2))
          #         plotspatialgraph(g1)
          #         plotspatialgraph(g2)
          #         plotreggraphs(RG3$graph1.g_aligned[[1]],RG3$graph2.g_aligned[[1]] )
          #         plotspatialgraph(mcs)
          
          d_v<- 1- ( network.size(mcs)/sqrt(network.size(g1)*network.size(g2)) )
          d_e<-1- (network.edgecount(mcs)/sqrt(network.edgecount(g1)*network.edgecount(g2)))
          gen_dv[i]<-d_v
          gen_de[i]<-d_e
        }
        
        imp_dv<-rep(0, times=dim(impMCS)[1])
        imp_de<-rep(0, times=dim(impMCS)[1])
        
        for(i in 1:dim(impMCS)[1]){
          
          mcs<-impMCS[[i]]
          
          i1<-impComp[i,1]
          i2<-impComp[i,2]
          g1<-graphList$graph[[i1]]
          g2<-graphList$graph[[i2]]
          d_v<- 1- ( network.size(mcs)/sqrt(network.size(g1)*network.size(g2)) )
          d_e<-1- (network.edgecount(mcs)/sqrt(network.edgecount(g1)*network.edgecount(g2)))
          imp_dv[i]<-d_v
          imp_de[i]<-d_e
        }
        
        dv_matrix[cid_n_index,cid_e_index]<-getEER_nomx(gen_dv, imp_dv)$eer
        de_matrix[cid_n_index,cid_e_index]<-getEER_nomx(gen_de, imp_de)$eer
        
      
      
     # print(c(cid_n,cid_e))
    }#end of cid_e loop
  }#end of cid_n loop
  print(m)
  print(dv_matrix)
  print(de_matrix)
  
  outputfile<-paste(opwd,biometric_list[b],"/", sourceid_list[b],"_dv_matrix_",mcsStruct[m],".Rdata", sep="")
  save(dv_matrix, file=outputfile)
  load(outputfile)
  
  outputfile<-paste(opwd,biometric_list[b],"/", sourceid_list[b],"_de_matrix_",mcsStruct[m],".Rdata", sep="")
  save(de_matrix, file=outputfile)
  load(outputfile)
  
}#end of m loop 
  
#   which(dv_matrix==min(dv_matrix), arr.ind=TRUE)
#   which(de_matrix==min(de_matrix), arr.ind=TRUE)
#   #save the matrices in VascularGraphsOutput
#   
  
  
  
  
}#end of b loop



