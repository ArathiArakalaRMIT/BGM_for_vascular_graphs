# Putting all the parameters from the training experiments, run on the retina test dataset to get performance
# Output must be mcs, as we will then look at graph topological measures.

wd<-getwd()
wd<-"H:/Uni Stuff/Uni Stuff/2015/Research Projects/ARC Discovery/R Code/Generic Code AA/BGM_Retina/"
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
source("graphMatchingFunctions_parallal_vpac.R")
source("graphTopologicalFeatures_vpac.R")
source("graphOperations_vpac.R")
source("graphMatchFunc_lib.R")


biometric<-"RE" 
sourceid <- "AutoHao"

# wdname<-paste("H:/Arathi/Graph Databases/Retina Graphs/RetinaGraphs",sourceid, sep="")
# setwd(wdname)

graphList<-numeric()
wd<-"H:/Uni Stuff/Uni Stuff/2015/Research Projects/ARC Discovery/Outputs/RetinaGraphs"

outputfile<-paste(wd,"/", sourceid, "/graphList_",sourceid,".Rdata", sep="")
#save(graphList, file=outputfile)
load(outputfile)

outputfile<-paste(wd,"/",sourceid[1],"/RegTests/trainingPersonList.Rdata", sep="")
#save(trainingPersonList, file=outputfile)
load(outputfile)

outputfile<-paste(wd,"/",sourceid[1],"/RegTests/testingPersonList.Rdata", sep="")
#save(testingPersonList, file=outputfile)
load(outputfile)


outputfile<-paste(wd,"/",sourceid[1],"/RegTests/gencomp_train.Rdata",sep="")
#save(gencomp_train, file=outputfile)
load(outputfile)
outputfile<-paste(wd,"/",sourceid[1],"/RegTests/impcomp_train.Rdata", sep="")
#save(impcomp_train, file=outputfile)
load(outputfile)
outputfile<-paste(wd,"/",sourceid[1],"/RegTests/gencomp_test.Rdata", sep="")
#save(gencomp_test, file=outputfile)
load(outputfile)
outputfile<-paste(wd,"/",sourceid[1],"/RegTests/impcomp_test.Rdata", sep="")
#save(impcomp_test, file=outputfile)
load(outputfile)

tol<-10
#cid<-5
L<-100
fplot<-0
maxtheta<-360
f_edgeinclude<-1
f_plot<-0
LR<-"R"

genComp<-gencomp_test
impComp<-impcomp_test

L<-100#from registration tests
cid_n<-3 #from graph matching on training data
cid_e<-9

cl1<-makeCluster(6, type = "SOCK")
registerDoSNOW(cl1)
# stopCluster(cl1)


    genRGscores<-numeric()
    genMCS_nodeBased<-numeric()
   
    
    for ( i in 1:dim(genComp)[1]){ 
      #for ( i in randIndGen){   
      i1<-genComp[i,1]
      i2<-genComp[i,2]
      
      g1<-graphList$graph[[i1]]
      
      person1<-graphList$person[[i1]]
      LR1<-graphList$LR[[i1]]
      Session1<-graphList$Session[[i1]]
      Nr1<-graphList$Nr[[i1]]
      
      g2<-graphList$graph[[i2]]
      
      person2<-graphList$person[[i2]]
      LR2<-graphList$LR[[i2]]
      Session2<-graphList$Session[[i2]]
      Nr2<-graphList$Nr[[i2]]
      
      RG3<-numeric()
      t3<-system.time(RG3<-registerGraphs_twostars(g1,g2,tol,maxtheta, f_plot, L))
      
      genRGscores <-c(genRGscores,RG3$bestScore[[1]]) 
      outputfile<-paste(wd,"/",sourceid[1],"/BetterGraphMatching/genRGscores_test_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tetss on test data
      save(genRGscores, file=outputfile)
      load(outputfile)
      
      
      RG<-RG3 
      
      
      ############row 1
      #       cid_n<-5
      #       cid_e<-5
      mcsSize<-numeric()
      for(k in 1:dim(RG)[1]){
        Path_euc<-getGraphEditPath_nodeBased(RG$graph1.g_aligned[[k]],RG$graph2.g_aligned[[k]],cid_n, cid_e,f_edgeinclude)
        MCS_euc<-getMCS_nodeBased(RG$graph1.g_aligned[[k]],RG$graph2.g_aligned[[k]],Path_euc$editPath)
        mcsSize<-c(mcsSize,network.size(MCS_euc))
        
      }
      k_max<-which(mcsSize==max(mcsSize))[1]
      Path_euc<-getGraphEditPath_nodeBased(RG$graph1.g_aligned[[k_max]],RG$graph2.g_aligned[[k_max]],cid_n, cid_e, f_edgeinclude)
      MCS_euc<-getMCS_nodeBased(RG$graph1.g_aligned[[k_max]],RG$graph2.g_aligned[[k_max]],Path_euc$editPath)
      p1<-MCS_euc 
      tmp<-list(mcs=p1)
      genMCS_nodeBased<-rbind(genMCS_nodeBased,tmp)
      #save this to file
      outputfile<-paste(wd,"/",sourceid[1],"/BetterGraphMatching/genMCS_nodeBased_test_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tetss on training data
      save(genMCS_nodeBased, file=outputfile)
      load(outputfile)
      
      
      
    
      print(i)
      
    }# end of i loop
    
    
    
    # do imposter comparisons
    
    impRGscores<-numeric()
    impMCS_nodeBased<-numeric()
 
    
    
    for ( i in 1:dim(impComp)[1]){
      #for ( i in randIndImp){
      #for(i in 1:15){
      
      i1<-impComp[i,1]
      i2<-impComp[i,2]
      
      g1<-graphList$graph[[i1]]
      
      person1<-graphList$person[[i1]]
      LR1<-graphList$LR[[i1]]
      Session1<-graphList$Session[[i1]]
      Nr1<-graphList$Nr[[i1]]
      
      g2<-graphList$graph[[i2]]
      
      person2<-graphList$person[[i2]]
      LR2<-graphList$LR[[i2]]
      Session2<-graphList$Session[[i2]]
      Nr2<-graphList$Nr[[i2]]
      
      
      
      RG3<-numeric()
      t3<-system.time(RG3<-registerGraphs_twostars(g1,g2,tol,maxtheta, f_plot, L))
      
      impRGscores <- c(impRGscores, RG3$bestScore[[1]])
      outputfile<-paste(wd,"/",sourceid[1],"/BetterGraphMatching/impRGscores_test_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tetss on training data
      save(impRGscores, file=outputfile)
      load(outputfile)
      
      
      RG<-RG3 
      
      
      ############row 1
      #       cid_n<-5
      #       cid_e<-5
      mcsSize<-numeric()
      for(k in 1:dim(RG)[1]){
        Path_euc<-getGraphEditPath_nodeBased(RG$graph1.g_aligned[[k]],RG$graph2.g_aligned[[k]],cid_n, cid_e,f_edgeinclude)
        MCS_euc<-getMCS_nodeBased(RG$graph1.g_aligned[[k]],RG$graph2.g_aligned[[k]],Path_euc$editPath)
        mcsSize<-c(mcsSize,network.size(MCS_euc))
        
      }
      k_max<-which(mcsSize==max(mcsSize))[1]
      Path_euc<-getGraphEditPath_nodeBased(RG$graph1.g_aligned[[k_max]],RG$graph2.g_aligned[[k_max]],cid_n, cid_e,f_edgeinclude)
      MCS_euc<-getMCS_nodeBased(RG$graph1.g_aligned[[k_max]],RG$graph2.g_aligned[[k_max]],Path_euc$editPath)
      p1<-MCS_euc 
      
      tmp<-list(mcs=p1)
      impMCS_nodeBased<-rbind(impMCS_nodeBased,tmp)
      #save this to file
      outputfile<-paste(wd,"/",sourceid[1],"/BetterGraphMatching/impMCS_nodeBased_test_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tetss on training data
      save(impMCS_nodeBased, file=outputfile)
      load(outputfile)
      
      
    
      print(i)
      
      
    }# end of i loop
    
stopCluster(cl1)

# read all mcses and get ROC curves for diff distance measures and linear combinations of them. (no SVM here)

outputfile<-paste(wd,"/",sourceid[1],"/BetterGraphMatching/genMCS_nodeBased_test_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tetss on training data
load(outputfile)
genMCS<-genMCS_nodeBased

outputfile<-paste(wd,"/",sourceid[1],"/BetterGraphMatching/impMCS_nodeBased_test_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tetss on training data
load(outputfile)
impMCS<-impMCS_nodeBased


for(i in 1:dim(genComp)[1]){
  i1<-genComp[i,1]
  i2<-genComp[i,2]  
  g1<-graphList$graph[[i1]]  
  g2<-graphList$graph[[i2]]
  
  mcs<-genMCS$mcs[[i]]
  
  
  
  
  
  
}#end of i loop

for(i in 1:dim(impComp)[1]){
  
  i1<-impComp[i,1]
  i2<-impComp[i,2]  
  g1<-graphList$graph[[i1]]  
  g2<-graphList$graph[[i2]]
  
}#end of j loop


#read all mcses
#########################################################

# outputfile<-paste(wd,"/",sourceid[1],"/BetterGraphMatching/genRGscores_train_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tetss on training data
# #save(genRGscores, file=outputfile)
# load(outputfile)
# 
# outputfile<-paste(wd,"/",sourceid[1],"/BetterGraphMatching/genMCS_nodeBased_train_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tetss on training data
# #save(genMCS_nodeBased, file=outputfile)
# load(outputfile)
# genMCS_nodeBased<-as.data.frame(genMCS_nodeBased)
# 
# outputfile<-paste(wd,"/",sourceid[1],"/BetterGraphMatching/genMCS_edgeBased_train_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tetss on training data
# #save(genMCS_edgeBased, file=outputfile)
# load(outputfile)
# genMCS_edgeBased<-as.data.frame(genMCS_edgeBased)
# 
# outputfile<-paste(wd,"/",sourceid[1],"/BetterGraphMatching/genMCS_starBased_train_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tetss on training data
# #save(genMCS_starBased, file=outputfile)
# load(outputfile)
# genMCS_starBased<-as.data.frame(genMCS_starBased)
# 
# #save this to file
# outputfile<-paste(wd,"/",sourceid[1],"/BetterGraphMatching/genMCS_twostarBased_train_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tetss on training data
# #save(genMCS_twostarBased, file=outputfile)
# load(outputfile)
# genMCS_twostarBased<-as.data.frame(genMCS_twostarBased)
# 
# 
# ###############################################################
# outputfile<-paste(wd,"/",sourceid[1],"/BetterGraphMatching/impRGscores_train_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tetss on training data
# #save(impRGscores, file=outputfile)
# load(outputfile)
# 
# outputfile<-paste(wd,"/",sourceid[1],"/BetterGraphMatching/impMCS_nodeBased_train_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tetss on training data
# #save(impMCS_nodeBased, file=outputfile)
# load(outputfile)
# impMCS_nodeBased<-as.data.frame(impMCS_nodeBased)
# 
# outputfile<-paste(wd,"/",sourceid[1],"/BetterGraphMatching/impMCS_edgeBased_train_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tetss on training data
# #save(impMCS_edgeBased, file=outputfile)
# load(outputfile)
# impMCS_edgeBased<-as.data.frame(impMCS_edgeBased)
# 
# outputfile<-paste(wd,"/",sourceid[1],"/BetterGraphMatching/impMCS_starBased_train_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tetss on training data
# #save(impMCS_starBased, file=outputfile)
# load(outputfile)
# impMCS_starBased<-as.data.frame(impMCS_starBased)
# 
# #save this to file
# outputfile<-paste(wd,"/",sourceid[1],"/BetterGraphMatching/impMCS_twostarBased_train_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tetss on training data
# #save(impMCS_twostarBased, file=outputfile)
# load(outputfile)
# impMCS_twostarBased<-as.data.frame(impMCS_twostarBased)
# 
# #############################################################
# 
# 
# genscores<-matrix(0, nrow=5, ncol=dim(genMCS_edgeBased)[1])
# impscores<-matrix(0, nrow=5, ncol=dim(impMCS_edgeBased)[1])
# for(i in 1:dim(genMCS_edgeBased)[1]){
#   i1<-genComp[i,1]
#   i2<-genComp[i,2]
#   g1<-graphList$graph[[i1]]
#   g2<-graphList$graph[[i2]]
#   
#   genscores[1,i]<-genRGscores[i]
#   
#   mcs<-genMCS_nodeBased$mcs[[i]]
#   #genscores[2,i]<-1-network.edgecount(mcs)/sqrt( network.edgecount(g1)*network.edgecount(g2) )  
#   genscores[2,i]<-1-network.size(mcs)/sqrt( network.size(g1)*network.size(g2) ) 
#   
#   mcs<-genMCS_edgeBased$mcs[[i]]
#   #genscores[3,i]<-1-network.edgecount(mcs)/sqrt( network.edgecount(g1)*network.edgecount(g2) ) 
#   genscores[3,i]<-1-network.size(mcs)/sqrt( network.size(g1)*network.size(g2) ) 
#   
#   
#   mcs<-genMCS_starBased$mcs[[i]]
#   #genscores[4,i]<-1-network.edgecount(mcs)/sqrt( network.edgecount(g1)*network.edgecount(g2) ) 
#   genscores[4,i]<-1-network.size(mcs)/sqrt( network.size(g1)*network.size(g2) ) 
#   
#   
#   mcs<-genMCS_twostarBased$mcs[[i]]
#   #genscores[5,i]<-1-network.edgecount(mcs)/sqrt( network.edgecount(g1)*network.edgecount(g2) ) 
#   genscores[5,i]<-1-network.size(mcs)/sqrt( network.size(g1)*network.size(g2) ) 
#   
# }#end of i loop
# 
# for(i in 1:dim(impMCS_edgeBased)[1]){
#   i1<-impComp[i,1]
#   i2<-impComp[i,2]
#   g1<-graphList$graph[[i1]]
#   g2<-graphList$graph[[i2]]
#   impscores[1,i]<-impRGscores[i]
#   
#   mcs<-impMCS_nodeBased$mcs[[i]]
#   #impscores[2,i]<-1-network.edgecount(mcs)/sqrt( network.edgecount(g1)*network.edgecount(g2) )  
#   impscores[2,i]<-1-network.size(mcs)/sqrt( network.size(g1)*network.size(g2) )  
#   
#   mcs<-impMCS_edgeBased$mcs[[i]]
#   #impscores[3,i]<-1-network.edgecount(mcs)/sqrt( network.edgecount(g1)*network.edgecount(g2) ) 
#   impscores[3,i]<-1-network.size(mcs)/sqrt( network.size(g1)*network.size(g2) ) 
#   
#   mcs<-impMCS_starBased$mcs[[i]]
#   #impscores[4,i]<-1-network.edgecount(mcs)/sqrt( network.edgecount(g1)*network.edgecount(g2) ) 
#   impscores[4,i]<-1-network.size(mcs)/sqrt( network.size(g1)*network.size(g2) ) 
#   
#   mcs<-impMCS_twostarBased$mcs[[i]]
#   #impscores[5,i]<-1-network.edgecount(mcs)/sqrt( network.edgecount(g1)*network.edgecount(g2) ) 
#   impscores[5,i]<-1-network.size(mcs)/sqrt( network.size(g1)*network.size(g2) ) 
#   
#   
# }#end of i loop
# 
# genRange<-1:dim(genMCS_edgeBased)[1]
# impRange<-1:dim(impMCS_edgeBased)[1]
# 
# windows()
# par(mfrow=c(3,2))
# roc<-getROC_unnorm_nomx(genscores[1,genRange], impscores[1,impRange],1)
# roc<-getROC_unnorm_nomx(genscores[2,genRange], impscores[2,impRange],1)
# roc<-getROC_unnorm_nomx(genscores[3,genRange], impscores[3,impRange],1)
# roc<-getROC_unnorm_nomx(genscores[4,genRange], impscores[4,impRange],1)
# roc<-getROC_unnorm_nomx(genscores[5,genRange], impscores[5,impRange],1)
# 
# 
# windows()
# par(mfrow=c(2,2))
# plothist(genscores[1,genRange], impscores[1,impRange])
# plothist(genscores[2,genRange], impscores[2,impRange])
# plothist(genscores[3,genRange], impscores[3,impRange])
# plothist(genscores[4,genRange], impscores[4,impRange])
# windows()
# plothist(genscores[5,genRange], impscores[5,impRange])
# 
# #find the genuine comparisons that do badly
# genCompToTest<-which(genscores[1,randGenInd]>min(impscores[1,randImpInd]))
# 
# genscores[,50]
# 
# 
# c( median(genscores[1, randGenInd]), median(impscores[1,randImpInd]), median(impscores[1,randImpInd])-median(genscores[1, randGenInd]) )
# c( median(genscores[2,randGenInd]), median(impscores[2,randImpInd]), median(impscores[2,randImpInd])-median(genscores[2, randGenInd]) )
# c( median(genscores[3, randGenInd]), median(impscores[3,randImpInd]), median(impscores[3,randImpInd])-median(genscores[3,randGenInd ]) )
# c( median(genscores[4, randGenInd]), median(impscores[4,randImpInd]), median(impscores[4,randImpInd])-median(genscores[4,randGenInd ]) )
# 
# c( max(genscores[1, randGenInd]), min(impscores[1,randImpInd]), min(impscores[1,randImpInd])-max(genscores[1,randGenInd ]) )
# c(  max(genscores[2,randGenInd ]), min(impscores[2,randImpInd]), min(impscores[2,randImpInd])-max(genscores[2,randGenInd ]) )
# c(  max(genscores[3, randGenInd]), min(impscores[3,randImpInd]), min(impscores[3,randImpInd])-max(genscores[3, randGenInd]) )
# c(  max(genscores[4, randGenInd]), min(impscores[4,randImpInd]), min(impscores[4,randImpInd])-max(genscores[4,randGenInd ]) )
