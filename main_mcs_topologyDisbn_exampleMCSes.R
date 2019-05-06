# This code will generate example mcses from training and test datasets to identify what features work well by observation.
wd<-"/Users/Arathi/Documents/2015/RMIT/ARC Discovery/R Code/Generic Code AA/BGM_Vascular/"
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
library(e1071)
#require(graphics)
#require(utils)
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

biometric_list<-c("RE", "PV", "WV", "HV","HV")
#sourceid_list<-c("AutoHao", "AutoHao2L", "AutoHao2R", "AutoHao7L_v1", "AutoHao7R_v1", "SFIR_LR", "SNIR_T_LR")
sourceid_list<-c("AutoHao", "AutoHao2L", "AutoHao7L_v1", "SFIR_LR", "SNIR_T_LR")
folder_list<-c("RetinaGraphs","PalmVeinGraphs", "WristVeinGraphs", "HandVeinGraphs", "HandVeinGraphs")
opwd<-"/Users/Arathi/Documents/2015/RMIT/ARC Discovery/Outputs/VascularGraphs/bestCid/"
ipwd<-"/Users/Arathi/Documents/2015/RMIT/ARC Discovery/Outputs"
bestStruct_testdata<-c("twostar", "edge", "edge", "edge", "edge", "twostar", "edge")
mcsStruct_testdata<-c("node", "node", "node", "node", "node")
cid_n_best<-c(3,5,3,3,3)
cid_e_best<-c(9,3,7,7,5)


for(b in 1:5){
  print("#######################################")
  print(b)
  
  # get training data
  ###############################################################
  
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
  
  
  
  #read all the training files and their cid_n and cid_e
  genComp<-gencomp_train
  impComp<-impcomp_train
  cid_n<-cid_n_best[b]
  cid_e<-cid_e_best[b]
  mcsStruct<-mcsStruct_testdata[b]
  
  outputfile<-paste(ipwd,"/",folder_list[b],"/",sourceid_list[b],"/BetterGraphMatching/genMCS_",mcsStruct,"Based_train_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tests on test data
  if(file.exists(outputfile)==FALSE) next
  load(outputfile)
  genMCS<-genMCS_nodeBased
  
  outputfile<-paste(ipwd,"/",folder_list[b],"/",sourceid_list[b],"/BetterGraphMatching/impMCS_",mcsStruct,"Based_train_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tests on test data
  if(file.exists(outputfile)==FALSE) next
  load(outputfile)
  impMCS<-impMCS_nodeBased
  
  
}#end of b loop