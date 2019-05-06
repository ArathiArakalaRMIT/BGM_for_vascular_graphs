# This code will build a model using feature selection and parameters selection done on training data and tested on the testing data.
#wd<-"H:/Uni Stuff/Uni Stuff/2015/Research Projects/ARC Discovery/R Code/Generic Code AA/BGM_Vascular/"
wd<-"/Users/Arathi/Documents/2016/RMIT/ARC Discovery/R Code/Generic Code AA/BGM_Vascular"
setwd(wd)

library(splancs)
library(spatial)
library(network)
library(clue)
library(iterators)
library(foreach)
library(doSNOW)
library(statnet.common)
library(utils)
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
opwd<-"/Users/Arathi/Documents/2016/RMIT/ARC Discovery/Outputs/VascularGraphs/bestCid/"
ipwd<-"/Users/Arathi/Documents/2016/RMIT/ARC Discovery/Outputs"
bestStruct_testdata<-c("twostar", "edge", "edge", "edge", "edge", "twostar", "edge")
mcsStruct_testdata<-c("node", "node", "node", "node", "node")
cid_n_best<-c(3,5,3,3,3)
cid_e_best<-c(9,3,7,7,5)


for(b in c(2,3,4,5,1)){
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
  
  
 
  
  #Get the vector for training data
  Size_gen<-rep(0, times=dim(genMCS)[1])
  Edgecount_gen<-rep(0, times=dim(genMCS)[1])
  Edgenoderatio_gen<-rep(0, times=dim(genMCS)[1])
  C1_gen<-rep(0, times=dim(genMCS)[1])
  C2_gen<-rep(0, times=dim(genMCS)[1])
  Isolated_gen<-rep(0, times=dim(genMCS)[1])
  C1_l_gen<-rep(0, times=dim(genMCS)[1])
  C1c2_gen<-rep(0, times=dim(genMCS)[1])
  ENRc1_gen<-rep(0, times=dim(genMCS)[1]) # number of edges in c1 to size of mcs
  Dmax_gen<-rep(0, times=dim(genMCS)[1])
  starSize_gen<-rep(0, times=dim(genMCS)[1])
  twostarSize_gen<-rep(0, times=dim(genMCS)[1])
  degreedisbn_gen<-matrix(0, nrow=dim(genMCS)[1], ncol=150)
  SpatialSize_gen<-rep(0, times=dim(genMCS)[1])
  
  
  for(i in 1:dim(genMCS)[1]){
    
   
    g<-genMCS[[i]]          
    i1<-genComp[i,1]
    i2<-genComp[i,2]
    g1<-graphList$graph[[i1]]
    g2<-graphList$graph[[i2]]
    
    Size_gen[i]<-1- (network.size(g)/sqrt(network.size(g1)*network.size(g2)) )
    Edgecount_gen[i]<-1- ( network.edgecount(g)/sqrt(network.edgecount(g1)*network.edgecount(g2)) )
    
    Edgenoderatio_gen[i]<-network.edgecount(g)/network.size(g)
    
    cd<-component.dist(g)
    largest_comp<-which(cd$csize==max(cd$csize))[1]
   
   
    
     C1_gen[i]<-1- ( (sort(cd$csize, decreasing=TRUE)[1]) /sqrt(network.size(g1)*network.size(g2)) )
     C2_gen[i]<-1- ( sort(cd$csize, decreasing=TRUE)[2]/ sqrt(network.size(g1)*network.size(g2)) )
     Isolated_gen[i]<-1- ( cd$cdist[1]/ sqrt(network.size(g1)*network.size(g2)) )
    
    n_c1<-which( cd$membership==largest_comp)
    n_todelete<-setdiff(1:network.size(g), n_c1)
    g_c1<-g
    delete.vertices(g_c1, n_todelete)
    edges_g_c1<-0
    if(network.edgecount(g_c1)>0) C1_l_gen[i]<-1-( sum( get.edge.value(g_c1,"edgelength") )/sqrt(sum(get.edge.value(g1,"edgelength"))*sum(get.edge.value(g2,"edgelength")) ) )
    C1c2_gen[i]<-1- ( (sort(cd$csize, decreasing=TRUE)[1]+sort(cd$csize, decreasing=TRUE)[2])/sqrt(network.size(g1)*network.size(g2)) )
    ENRc1_gen[i]<-network.edgecount(g_c1)/network.size(g_c1)
    
    deg<-degree(g, gmode="graph", cmode="indegree")
    Dmax_gen[i]<-max(deg)
    
    
    starSize_gen[i]<-1- ( dim(get_star(g))[1]/sqrt(dim(get_star(g1))[1]*dim(get_star(g2))[1])  )
    twostarSize_gen[i]<-1- ( dim(get_twostar(g))[1]/sqrt(dim(get_twostar(g1))[1]*dim(get_twostar(g2))[1])  )
    
    degree_tmp<-rep(0, times=max(deg) )
    for(j in 1:max(deg)){
      degree_tmp[j]<-length(which(deg==j))
    }
    degree_tmp<-degree_tmp/sum(degree_tmp)
    degreedisbn_gen[i,1:max(deg)]<-degree_tmp
    
    #to get the spatial size, get smallest and largest x and y coordinates and get area of rectangle
    #of mcs
    x_minmax<-range(get.vertex.attribute(g, "xcoord"))
    y_minmax<-range(get.vertex.attribute(g, "ycoord"))
    area_g<-(x_minmax[2]-x_minmax[1])*(y_minmax[2]-y_minmax[1])
    if(area_g<=0.0067) area_g<-exp(-5)
    #of g1
    x_minmax<-range(get.vertex.attribute(g1, "xcoord"))
    y_minmax<-range(get.vertex.attribute(g1, "ycoord"))
    area_g1<-(x_minmax[2]-x_minmax[1])*(y_minmax[2]-y_minmax[1])
    #of g2
    x_minmax<-range(get.vertex.attribute(g2, "xcoord"))
    y_minmax<-range(get.vertex.attribute(g2, "ycoord"))
    area_g2<-(x_minmax[2]-x_minmax[1])*(y_minmax[2]-y_minmax[1])
    
    SpatialSize_gen[i]<-1-( ( area_g)/(sqrt( area_g1)* sqrt(area_g2))  )
    #nodes per sq area
  #  SpatialSize_gen[i]<-(  (Size_gen[i]/area_g)/(sqrt( network.size(g1)/area_g1)* sqrt(network.size(g2)/area_g2))  )
    
    
    
#     windows()
#     par(mfrow=c(2,2), oma=c(0,0,2,0))
#     plotspatialgraph(g1)
#     plotspatialgraph(g2)
#     plotspatialgraph(g)
#     titlestring<-paste("MCS of graph ", i1, "and ", i2, sep=" ")
#     title(titlestring, outer=TRUE)
#       
  
    
  } # end of i loop
  
  
  Size_imp<-rep(0, times=dim(impMCS)[1])
  Edgecount_imp<-rep(0, times=dim(impMCS)[1])
  Edgenoderatio_imp<-rep(0, times=dim(impMCS)[1])
  C1_imp<-rep(0, times=dim(impMCS)[1])
  C2_imp<-rep(0, times=dim(impMCS)[1])
  Isolated_imp<-rep(0, times=dim(impMCS)[1])
  C1_l_imp<-rep(0, times=dim(impMCS)[1])
  C1c2_imp<-rep(0, times=dim(impMCS)[1])
  ENRc1_imp<-rep(0, times=dim(impMCS)[1]) # number of edges in c1 to size of mcs
  Dmax_imp<-rep(0, times=dim(impMCS)[1])
  starSize_imp<-rep(0, times=dim(impMCS)[1])
  twostarSize_imp<-rep(0, times=dim(impMCS)[1])
  degreedisbn_imp<-matrix(0, nrow=dim(impMCS)[1], ncol=150)
  SpatialSize_imp<-rep(0, times=dim(impMCS)[1])
  
  for(i in 1:dim(impMCS)[1]){
    
    g<-impMCS[[i]]
    
    i1<-impComp[i,1]
    i2<-impComp[i,2]
    g1<-graphList$graph[[i1]]
    g2<-graphList$graph[[i2]]
    
    
    Size_imp[i]<-1- (network.size(g)/sqrt(network.size(g1)*network.size(g2)) )
    Edgecount_imp[i]<-1- ( network.edgecount(g)/sqrt(network.edgecount(g1)*network.edgecount(g2)) )
    Edgenoderatio_imp[i]<-network.edgecount(g)/network.size(g)
    cd<-component.dist(g)
    largest_comp<-which(cd$csize==max(cd$csize))[1]
    
    C1_imp[i]<-1- ( (sort(cd$csize, decreasing=TRUE)[1]) /sqrt(network.size(g1)*network.size(g2)) )
    C2_imp[i]<-1- ( sort(cd$csize, decreasing=TRUE)[2]/ sqrt(network.size(g1)*network.size(g2)) )
    
    Isolated_imp[i]<-1- ( cd$cdist[1]/ sqrt(network.size(g1)*network.size(g2)) )
    
    n_c1<-which( cd$membership==largest_comp)
    n_todelete<-setdiff(1:network.size(g), n_c1)
    g_c1<-g
    delete.vertices(g_c1, n_todelete)
    edges_g_c1<-0
    if(network.edgecount(g_c1)>0) C1_l_imp[i]<-1-( sum( get.edge.value(g_c1,"edgelength") )/sqrt(sum(get.edge.value(g1,"edgelength"))*sum(get.edge.value(g2,"edgelength")) ) )
    C1c2_imp[i]<-1-( (sort(cd$csize, decreasing=TRUE)[1]+sort(cd$csize, decreasing=TRUE)[2])/ sqrt(network.size(g1)*network.size(g2)) )
    ENRc1_imp[i]<-network.edgecount(g_c1)/network.size(g_c1)
    
    deg<-degree(g, gmode="graph", cmode="indegree")
    Dmax_imp[i]<-max(deg)
    
    starSize_imp[i]<-1- ( dim(get_star(g))[1]/sqrt(dim(get_star(g1))[1]*dim(get_star(g2))[1])  )
    twostarSize_imp[i]<-1- ( dim(get_twostar(g))[1]/sqrt(dim(get_twostar(g1))[1]*dim(get_twostar(g2))[1])  )
    
    degree_tmp<-rep(0, times=max(deg) )
    for(j in 1:max(deg)){
      degree_tmp[j]<-length(which(deg==j))
    }
    degree_tmp<-degree_tmp/sum(degree_tmp)
    degreedisbn_imp[i,1:max(deg)]<-degree_tmp
    
    #to get the spatial size, get smallest and largest x and y coordinates and get area of rectangle
    #of mcs
    x_minmax<-range(get.vertex.attribute(g, "xcoord"))
    y_minmax<-range(get.vertex.attribute(g, "ycoord"))
    area_g<-(x_minmax[2]-x_minmax[1])*(y_minmax[2]-y_minmax[1])
    if(area_g<=0.0067) area_g<-exp(-5)
    #of g1
    x_minmax<-range(get.vertex.attribute(g1, "xcoord"))
    y_minmax<-range(get.vertex.attribute(g1, "ycoord"))
    area_g1<-(x_minmax[2]-x_minmax[1])*(y_minmax[2]-y_minmax[1])
    #of g2
    x_minmax<-range(get.vertex.attribute(g2, "xcoord"))
    y_minmax<-range(get.vertex.attribute(g2, "ycoord"))
    area_g2<-(x_minmax[2]-x_minmax[1])*(y_minmax[2]-y_minmax[1])
    
    SpatialSize_imp[i]<-1-( ( area_g)/(sqrt( area_g1)* sqrt(area_g2))  )
  #  SpatialSize_imp[i]<-(  (Size_imp[i]/area_g)/(sqrt( network.size(g1)/area_g1)* sqrt(network.size(g2)/area_g2))  )
    
#     
#         windows()
#         par(mfrow=c(2,2), oma=c(0,0,2,0))
#         plotspatialgraph(g1)
#         plotspatialgraph(g2)
#         plotspatialgraph(g)
#         titlestring<-paste("MCS of graph ", i1, "and ", i2, sep=" ")
#         title(titlestring, outer=TRUE)
#           
    
  }
  
  quartz()
  plot(1:length(Size_gen), Edgecount_gen, type="p", col=3, xlim=c(0,length(Size_imp)), ylim=c(0,1))
  points(1:length(Size_imp), Edgecount_imp, col=2 )
  
  C2_imp[which(is.na(C2_imp)==TRUE)]<-1
  C1c2_imp[which(is.na(C1c2_imp)==TRUE)]<-1
  C2_gen[which(is.na(C2_gen)==TRUE)]<-1
  C1c2_gen[which(is.na(C1c2_gen)==TRUE)]<-1
  
#   
   gen_data_train<-cbind( Size_gen,  Edgecount_gen, Edgenoderatio_gen, C1_gen,  C2_gen,  Isolated_gen, C1_l_gen, C1c2_gen,  ENRc1_gen,  Dmax_gen, starSize_gen, twostarSize_gen, SpatialSize_gen)
#   #pairs(gen_data)
   imp_data_train<-cbind( Size_imp,  Edgecount_imp, Edgenoderatio_imp, C1_imp,  C2_imp,  Isolated_imp, C1_l_imp, C1c2_imp,  ENRc1_imp,  Dmax_imp, starSize_imp, twostarSize_imp, SpatialSize_imp)
  class<-rep(1, times=length(Size_gen))
  genscores_all<-cbind(gen_data_train, class)
  class<-rep(2, times=length(Size_imp))
  impscores_all<-cbind(imp_data_train, class)
  
#   # print out a selection of 16 (4x4) genuine and imposter mcses with their vectors
#   printIndex_gen<-sample(1:dim(gen_data_train)[1], size=16, replace=FALSE)
#   printIndex_imp<-sample(1:dim(imp_data_train)[1], size=16, replace=FALSE)
#   windows()
#   par(mfrow=c(2,2))
#   for(i in printIndex_gen[1:4]){
#     g<-genMCS[[i]]   
#     i1<-genComp[i,1]
#     i2<-genComp[i,2]
#     titlestring<-paste("MCS of graph ", i1, "and ", i2, sep=" ")
#     plotspatialgraph(g)
#     title(titlestring)
#   
#   }
#   windows()
#   par(mfrow=c(2,2))
#   for(i in printIndex_gen[5:8]){
#     g<-genMCS[[i]]   
#     i1<-genComp[i,1]
#     i2<-genComp[i,2]
#     titlestring<-paste("MCS of graph ", i1, "and ", i2, sep=" ")
#     plotspatialgraph(g)
#     title(titlestring)
#     
#   }
#   
#   windows()
#   par(mfrow=c(2,2))
#   for(i in printIndex_gen[9:12]){
#     g<-genMCS[[i]]   
#     i1<-genComp[i,1]
#     i2<-genComp[i,2]
#     titlestring<-paste("MCS of graph ", i1, "and ", i2, sep=" ")
#     plotspatialgraph(g)
#     title(titlestring)
#     
#   }
#   
# 
# genD<-gen_data_train[printIndex_gen,]
#   
#   windows()
#   par(mfrow=c(2,2))
#   for(i in printIndex_imp[1:4]){
#     g<-impMCS[[i]]   
#     i1<-impComp[i,1]
#     i2<-impComp[i,2]
#     titlestring<-paste("MCS of graph ", i1, "and ", i2, sep=" ")
#     plotspatialgraph(g)
#     title(titlestring)
#     
#   }
#   
#   windows()
#   par(mfrow=c(2,2))
#   for(i in printIndex_imp[5:8]){
#     g<-impMCS[[i]]   
#     i1<-impComp[i,1]
#     i2<-impComp[i,2]
#     titlestring<-paste("MCS of graph ", i1, "and ", i2, sep=" ")
#     plotspatialgraph(g)
#     title(titlestring)
#     
#   }
#   windows()
#   par(mfrow=c(2,2))
#   for(i in printIndex_imp[9:12]){
#     g<-impMCS[[i]]   
#     i1<-impComp[i,1]
#     i2<-impComp[i,2]
#     titlestring<-paste("MCS of graph ", i1, "and ", i2, sep=" ")
#     plotspatialgraph(g)
#     title(titlestring)
#     
#   }
#   impD<-imp_data_train[printIndex_imp,]
#   
#   genscores_all<-cbind(Size=Size_gen, Edgecount=Edgecount_gen, Dmax=Dmax_gen, C1c2=C1c2_gen, C1_l=C1_l_gen, starSize=starSize_gen)
#   impscores_all<-cbind(Size=Size_imp, Edgecount=Edgecount_imp, Dmax=Dmax_imp, C1c2=C1c2_imp, C1_l=C1_l_imp, starSize=starSize_imp)
#   class<-rep(1, times=length(Size_gen))
#   genscores_all<-cbind(genscores_all, class)
#   class<-rep(2, times=length(Size_imp))
#   impscores_all<-cbind(impscores_all, class)
#   
#     
#   scores_train<-rbind(genscores_all, impscores_all)
#   scores_train<-as.data.frame(scores_train)
#   
#   model<-svm(class~.,data=scores, cost=as.numeric(model00[[1]]),kernel="linear" )
#   summary(model)
#   
#   x<-subset(scores, select=-class)
#   y<-subset(scores, select=class)
#   
#   
#   # compute decision values and probabilities:
#   pred<-predict(model, x, decision.values = TRUE)
#   svmScores<-attr(pred, "decision.values")
#   svmScores_gen<-svmScores[1:dim(genscores)[1]]
#   svmScores_imp<-svmScores[-(1:dim(genscores)[1])]
#   
#   roc_svm<-getROC_sim_nomx( svmScores_gen, svmScores_imp,0)
#   roc_Size<-getROC_unnorm_nomx(Size_gen, Size_imp,0)
#   roc_EdgeCount<-getROC_unnorm_nomx(Edgecount_gen, Edgecount_imp,0)
#   
#   legendtext<-c(paste("SVM scores, EER=", round(roc_svm$eer,4), sep=""),paste("Size, EER=", round(roc_Size$eer,4), sep=""),paste("Edgecount, EER=", round(roc_EdgeCount$eer,4), sep=""))
#   
#   quartz()
#   #plot(roc_svm$x, roc_svm$y, col=1, lwd=2, type='l', lty=1, ylim=c(0,1), xlim=c(0,1))
#   plot(roc_Size$x, roc_Size$y, col=2, lwd=2, type='l', lty=2, ylim=c(0,1), xlim=c(0,1))
#   points(roc_EdgeCount$x, roc_EdgeCount$y, col=3, lwd=2, type='l', lty=3)
#   legend(0.5,1, legendtext, col=1:3, lwd=2, lty=1:3)
# 
  
  ##########################################################################################
  # Get the vector for test data
  genComp<-gencomp_test
  impComp<-impcomp_test
  cid_n<-cid_n_best[b]
  cid_e<-cid_e_best[b]
  mcsStruct<-mcsStruct_testdata[b]
  
  #read all the test files and their cid_n and cid_e
  outputfile<-paste(ipwd,"/",folder_list[b],"/",sourceid_list[b],"/BetterGraphMatching/genMCS_",mcsStruct,"Based_test_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tests on test data
  if(file.exists(outputfile)==FALSE) next
  load(outputfile)
  genMCS<-genMCS_nodeBased
  
  outputfile<-paste(ipwd,"/",folder_list[b],"/",sourceid_list[b],"/BetterGraphMatching/impMCS_",mcsStruct,"Based_test_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tests on test data
  if(file.exists(outputfile)==FALSE) next
  load(outputfile)
  impMCS<-impMCS_nodeBased
  
  if(b==2){
    impMCS<-numeric()
    outputfile<-paste(ipwd,"/",folder_list[b],"/",sourceid_list[b],"/BetterGraphMatching/impMCS_",mcsStruct,"Based_test_",cid_n,"_",cid_e,"_1to3386.Rdata", sep="") #cost function tests on test data
    if(file.exists(outputfile)==FALSE) next
    load(outputfile)
    impMCS<-rbind(impMCS, impMCS_nodeBased)
    
    outputfile<-paste(ipwd,"/",folder_list[b],"/",sourceid_list[b],"/BetterGraphMatching/impMCS_",mcsStruct,"Based_test_",cid_n,"_",cid_e,"_3387to9601.Rdata", sep="") #cost function tests on test data
    if(file.exists(outputfile)==FALSE) next
    load(outputfile)
    impMCS<-rbind(impMCS, impMCS_nodeBased)
    
  }
          
  
  Size_gen<-rep(0, times=dim(genMCS)[1])
  Edgecount_gen<-rep(0, times=dim(genMCS)[1])
  Edgenoderatio_gen<-rep(0, times=dim(genMCS)[1])
  C1_gen<-rep(0, times=dim(genMCS)[1])
  C2_gen<-rep(0, times=dim(genMCS)[1])
  Isolated_gen<-rep(0, times=dim(genMCS)[1])
  C1_l_gen<-rep(0, times=dim(genMCS)[1])
  C1c2_gen<-rep(0, times=dim(genMCS)[1])
  ENRc1_gen<-rep(0, times=dim(genMCS)[1]) # number of edges in c1 to size of mcs
  Dmax_gen<-rep(0, times=dim(genMCS)[1])
  starSize_gen<-rep(0, times=dim(genMCS)[1])
  twostarSize_gen<-rep(0, times=dim(genMCS)[1])
  degreedisbn_gen<-matrix(0, nrow=dim(genMCS)[1], ncol=150)
  SpatialSize_gen<-rep(0, times=dim(genMCS)[1])
  
  
  for(i in 1:dim(genMCS)[1]){
        
   
          g<-genMCS[[i]]          
          i1<-genComp[i,1]
          i2<-genComp[i,2]
          g1<-graphList$graph[[i1]]
          g2<-graphList$graph[[i2]]
          
          Size_gen[i]<-1- (network.size(g)/sqrt(network.size(g1)*network.size(g2)) )
          Edgecount_gen[i]<-1- ( network.edgecount(g)/sqrt(network.edgecount(g1)*network.edgecount(g2)) )
          
          Edgenoderatio_gen[i]<-network.edgecount(g)/network.size(g)
          
          cd<-component.dist(g)
          largest_comp<-which(cd$csize==max(cd$csize))[1]
          
          C1_gen[i]<-1- ( (sort(cd$csize, decreasing=TRUE)[1]) /sqrt(network.size(g1)*network.size(g2)) )
          C2_gen[i]<-1- ( sort(cd$csize, decreasing=TRUE)[2]/ sqrt(network.size(g1)*network.size(g2)) )
          Isolated_gen[i]<-1- ( cd$cdist[1]/ sqrt(network.size(g1)*network.size(g2)) )
          
          n_c1<-which( cd$membership==largest_comp)
          n_todelete<-setdiff(1:network.size(g), n_c1)
          g_c1<-g
          delete.vertices(g_c1, n_todelete)
          edges_g_c1<-0
          if(network.edgecount(g_c1)>0) C1_l_gen[i]<-1-( sum( get.edge.value(g_c1,"edgelength") )/sqrt(sum(get.edge.value(g1,"edgelength"))*sum(get.edge.value(g2,"edgelength")) ) )
          C1c2_gen[i]<-1- ( (sort(cd$csize, decreasing=TRUE)[1]+sort(cd$csize, decreasing=TRUE)[2])/sqrt(network.size(g1)*network.size(g2)) )
          ENRc1_gen[i]<-network.edgecount(g_c1)/network.size(g_c1)
          
          deg<-degree(g, gmode="graph", cmode="indegree")
          Dmax_gen[i]<-max(deg)
          
          starSize_gen[i]<-1- ( dim(get_star(g))[1]/sqrt(dim(get_star(g1))[1]*dim(get_star(g2))[1])  )
          twostarSize_gen[i]<-1- ( dim(get_twostar(g))[1]/sqrt(dim(get_twostar(g1))[1]*dim(get_twostar(g2))[1])  )
          
          degree_tmp<-rep(0, times=max(deg) )
          for(j in 1:max(deg)){
            degree_tmp[j]<-length(which(deg==j))
          }
          degree_tmp<-degree_tmp/sum(degree_tmp)
          degreedisbn_gen[i,1:max(deg)]<-degree_tmp
          
          
          #to get the spatial size, get smallest and largest x and y coordinates and get area of rectangle
          #of mcs
          x_minmax<-range(get.vertex.attribute(g, "xcoord"))
          y_minmax<-range(get.vertex.attribute(g, "ycoord"))
          area_g<-(x_minmax[2]-x_minmax[1])*(y_minmax[2]-y_minmax[1])
          if(area_g<=0.0067) area_g<-exp(-5)
          #of g1
          x_minmax<-range(get.vertex.attribute(g1, "xcoord"))
          y_minmax<-range(get.vertex.attribute(g1, "ycoord"))
          area_g1<-(x_minmax[2]-x_minmax[1])*(y_minmax[2]-y_minmax[1])
          #of g2
          x_minmax<-range(get.vertex.attribute(g2, "xcoord"))
          y_minmax<-range(get.vertex.attribute(g2, "ycoord"))
          area_g2<-(x_minmax[2]-x_minmax[1])*(y_minmax[2]-y_minmax[1])
          
          SpatialSize_gen[i]<-1-( ( area_g)/(sqrt( area_g1)* sqrt(area_g2))  )
          #nodes per sq area
 #         SpatialSize_gen[i]<-(  (Size_gen[i]/area_g)/(sqrt( network.size(g1)/area_g1)* sqrt(network.size(g2)/area_g2))  )
          
  

        } # end of i loop
        
  
  Size_imp<-rep(0, times=dim(impMCS)[1])
  Edgecount_imp<-rep(0, times=dim(impMCS)[1])
  Edgenoderatio_imp<-rep(0, times=dim(impMCS)[1])
  C1_imp<-rep(0, times=dim(impMCS)[1])
  C2_imp<-rep(0, times=dim(impMCS)[1])
  Isolated_imp<-rep(0, times=dim(impMCS)[1])
  C1_l_imp<-rep(0, times=dim(impMCS)[1])
  C1c2_imp<-rep(0, times=dim(impMCS)[1])
  ENRc1_imp<-rep(0, times=dim(impMCS)[1]) # number of edges in c1 to size of mcs
  Dmax_imp<-rep(0, times=dim(impMCS)[1])
  starSize_imp<-rep(0, times=dim(impMCS)[1])
  twostarSize_imp<-rep(0, times=dim(impMCS)[1])
  degreedisbn_imp<-matrix(0, nrow=dim(impMCS)[1], ncol=150)
  SpatialSize_imp<-rep(0, times=dim(impMCS)[1])
  
        for(i in 1:dim(impMCS)[1]){
     
         g<-impMCS[[i]]
          
          i1<-impComp[i,1]
          i2<-impComp[i,2]
          g1<-graphList$graph[[i1]]
          g2<-graphList$graph[[i2]]
          
          
          Size_imp[i]<-1- (network.size(g)/sqrt(network.size(g1)*network.size(g2)) )
          Edgecount_imp[i]<-1- ( network.edgecount(g)/sqrt(network.edgecount(g1)*network.edgecount(g2)) )
          Edgenoderatio_imp[i]<-network.edgecount(g)/network.size(g)
          cd<-component.dist(g)
          largest_comp<-which(cd$csize==max(cd$csize))[1]
          
          C1_imp[i]<-1- ( (sort(cd$csize, decreasing=TRUE)[1]) /sqrt(network.size(g1)*network.size(g2)) )
          C2_imp[i]<-1- ( sort(cd$csize, decreasing=TRUE)[2]/ sqrt(network.size(g1)*network.size(g2)) )
         
          Isolated_imp[i]<-1- ( cd$cdist[1]/ sqrt(network.size(g1)*network.size(g2)) )
          
          n_c1<-which( cd$membership==largest_comp)
          n_todelete<-setdiff(1:network.size(g), n_c1)
          g_c1<-g
          delete.vertices(g_c1, n_todelete)
          edges_g_c1<-0
          if(network.edgecount(g_c1)>0) C1_l_imp[i]<-1-( sum( get.edge.value(g_c1,"edgelength") )/sqrt(sum(get.edge.value(g1,"edgelength"))*sum(get.edge.value(g2,"edgelength")) ) )
          C1c2_imp[i]<-1-( (sort(cd$csize, decreasing=TRUE)[1]+sort(cd$csize, decreasing=TRUE)[2])/ sqrt(network.size(g1)*network.size(g2)) )
          ENRc1_imp[i]<-network.edgecount(g_c1)/network.size(g_c1)
          
          deg<-degree(g, gmode="graph", cmode="indegree")
          Dmax_imp[i]<-max(deg)
          
         starSize_imp[i]<-1- ( dim(get_star(g))[1]/sqrt(dim(get_star(g1))[1]*dim(get_star(g2))[1])  )
         twostarSize_imp[i]<-1- ( dim(get_twostar(g))[1]/sqrt(dim(get_twostar(g1))[1]*dim(get_twostar(g2))[1])  )
         
          degree_tmp<-rep(0, times=max(deg) )
          for(j in 1:max(deg)){
            degree_tmp[j]<-length(which(deg==j))
          }
          degree_tmp<-degree_tmp/sum(degree_tmp)
          degreedisbn_imp[i,1:max(deg)]<-degree_tmp
          
         
         #to get the spatial size, get smallest and largest x and y coordinates and get area of rectangle
         #of mcs
         x_minmax<-range(get.vertex.attribute(g, "xcoord"))
         y_minmax<-range(get.vertex.attribute(g, "ycoord"))
         area_g<-(x_minmax[2]-x_minmax[1])*(y_minmax[2]-y_minmax[1])
         if(area_g<=0.0067) area_g<-exp(-5)
         #of g1
         x_minmax<-range(get.vertex.attribute(g1, "xcoord"))
         y_minmax<-range(get.vertex.attribute(g1, "ycoord"))
         area_g1<-(x_minmax[2]-x_minmax[1])*(y_minmax[2]-y_minmax[1])
         #of g2
         x_minmax<-range(get.vertex.attribute(g2, "xcoord"))
         y_minmax<-range(get.vertex.attribute(g2, "ycoord"))
         area_g2<-(x_minmax[2]-x_minmax[1])*(y_minmax[2]-y_minmax[1])
         
         SpatialSize_imp[i]<-1-( ( area_g)/(sqrt( area_g1)* sqrt(area_g2))  )
#         SpatialSize_imp[i]<-(  (Size_imp[i]/area_g)/(sqrt( network.size(g1)/area_g1)* sqrt(network.size(g2)/area_g2))  )
                  
          
          
        }
  
  quartz()
  plot(1:length(Size_gen), Size_gen, type="p", col=3, xlim=c(0,length(Size_imp)), ylim=c(0,1))
  points(1:length(Size_imp), Size_imp, col=2 )
  
 indToChk<-which(Size_gen>0.7)
  
 C2_imp[which(is.na(C2_imp)==TRUE)]<-1
 C1c2_imp[which(is.na(C1c2_imp)==TRUE)]<-1
  C2_gen[which(is.na(C2_gen)==TRUE)]<-1
  C1c2_gen[which(is.na(C1c2_gen)==TRUE)]<-1
  
#   # print out a selection of 16 (4x4) genuine and imposter mcses with their vectors
#   printIndex_gen<-sample(1:dim(gen_data_test)[1], size=16, replace=FALSE)
#   printIndex_imp<-sample(1:dim(imp_data_test)[1], size=16, replace=FALSE)
#   windows()
#   par(mfrow=c(2,2))
#   for(i in printIndex_gen[9:12]){
#     g<-genMCS[[i]]   
#     i1<-genComp[i,1]
#     i2<-genComp[i,2]
#     titlestring<-paste("MCS of graph ", i1, "and ", i2, sep=" ")
#     plotspatialgraph(g)
#     title(titlestring)
#     
#   }
#   
#   genD<-gen_data_test[printIndex_gen,]
#   
#   windows()
#   par(mfrow=c(2,2))
#   for(i in printIndex_imp[1:4]){
#     g<-impMCS[[i]]   
#     i1<-impComp[i,1]
#     i2<-impComp[i,2]
#     titlestring<-paste("MCS of graph ", i1, "and ", i2, sep=" ")
#     plotspatialgraph(g)
#     title(titlestring)
#     
#   }
#   
#   impD<-imp_data_test[printIndex_imp,]
#   
  
# ############################################  
#   windows()
#   par(mfrow=c(2,2))
#   plothist(Size_gen, Size_imp)
#   title(sub="size")
#   plothist(Edgecount_gen,Edgecount_imp)
#   title(sub="Edgecount")
#   plothist(Edgenoderatio_gen, Edgenoderatio_imp)
#   title(sub="Edgenoderatio")
#   plothist(C1_gen, C1_imp)
#   title(sub="C1")
#   
#   windows()
#   par(mfrow=c(2,2))
#   plothist(C2_gen, C2_imp)
#   title(sub="C2")
#   plothist(Isolated_gen,Isolated_imp)
#   title(sub="Isolated")
#   plothist(starSize_gen, starSize_imp)
#   title(sub="starSize")
#   plothist(twostarSize_gen, twostarSize_imp)
#   title(sub="twostarSize")
  
#   windows()
#   par(mfrow=c(2,2))
#   plothist(C1c2_gen, C1c2_imp)
#   title(sub="C1c2")
#   plothist(ENRc1_gen,ENRc1_imp)
#   title(sub="ENRc1")
#   hist(C1_l_gen)
#   hist(C1_l_imp)
#   
#    
#   windows()
#   par(mfrow=c(1,2))
#   hist(Dmax_gen)
#   title(sub=paste("#degree 3 graphs= ", length(which(Dmax_gen==3)) , sep="") )
#   hist(Dmax_imp)
#   title(sub=paste("#degree 3 graphs= ", length(which(Dmax_imp==3)) , sep="") )
  
    
#   windows()
#   plot(1:5, (colSums(degreedisbn_gen )/sum(degreedisbn_gen))[1:5], type='b', lwd=2, col=3, ylim=c(0,1))
#   lines( 1:5, (colSums(degreedisbn_imp)/sum(degreedisbn_imp))[1:5], type='b', lwd=2, col=2)
##################################################
#   gen_data<-cbind( Size_gen,  Edgecount_gen, Edgenoderatio_gen, C1_gen,  C2_gen,  Isolated_gen, C1_l_gen, C1c2_gen,  ENRc1_gen,  Dmax_gen, starSize_gen, twostarSize_gen)
# #pairs(gen_data)
#   imp_data<-cbind( Size_imp,  Edgecount_imp, Edgenoderatio_imp, C1_imp,  C2_imp,  Isolated_imp, C1_l_imp, C1c2_imp,  ENRc1_imp,  Dmax_imp, starSize_imp, twostarSize_imp)
  
#############################   
  #1
#  
#   gen<-cbind(Edgecount_gen,C1_l_gen)
#   imp<-cbind(Edgecount_imp,C1_l_imp)
#   
# windows()
# plot(gen, pch=1, col=Dmax_gen, xlim=c(0,1), ylim=c(0,1), cex=1)
#   points(imp, pch=2, col=Dmax_imp)
#   title(main= cor(Edgecount_gen,C1_l_gen, method="spearman"), sub= cor(Edgecount_gen,Dmax_gen, method="spearman"))
#   
#   
#   #2
#   
#   gen<-cbind(Size_gen,Edgecount_gen)
#   imp<-cbind(Size_imp,Edgecount_imp)
#   
#   windows()
#   plot(gen, pch=1, col=Dmax_gen, xlim=c(0,1), ylim=c(0,2), cex=1)
#   points(imp, pch=2, col=Dmax_imp)
#   title(main= cor(Size_gen,Edgecount_gen, method="spearman"))
#   
#   #3
#   gen<-cbind(Size_gen,C1c2_gen)
#   imp<-cbind(Size_imp,C1c2_imp)
#   
#   windows()
#   plot(gen, pch=1, col=Dmax_gen, xlim=c(0,1), ylim=c(0,2), cex=1)
#   points(imp, pch=2, col=Dmax_imp)
#   title(main= cor(Size_gen,C1c2_gen, method="spearman"))
#   
#   
#   
#   #4
#   gen<-cbind(Size_gen,C1_l_gen)
#   imp<-cbind(Size_imp,C1_l_imp)
#   
#   windows()
#   plot(gen, pch=1, col=Dmax_gen, xlim=c(0,1), ylim=c(0,2), cex=1)
#   points(imp, pch=2, col=Dmax_imp)
#   title(main= cor(Size_gen,C1_l_gen, method="spearman"))
#   
#   
#   #5
#   gen<-cbind(Size_gen,starSize_gen)
#   imp<-cbind(Size_imp,starSize_imp)
#   
#   windows()
#   plot(gen, pch=1, col=Dmax_gen, xlim=c(0,1), ylim=c(0,2), cex=1)
#   points(imp, pch=2, col=Dmax_imp)
#   title(main= cor(Size_gen,C1_l_gen, method="spearman"))
#   
#   
# 
#   ############
#   ### Now try classification using decision tree with thresholds from training data
#   class<-rep(0, times=length(Size_gen))
#   if(b==1){
#     genscores<-cbind(Size=Size_gen, Dmax=Dmax_gen, C1_l=C1_l_gen, starSize=starSize_gen, C1c2=C1c2_gen)
#   }  
#   if(b==2){
#     genscores<-cbind(Size=Size_gen, Dmax=Dmax_gen, C1_l=C1_l_gen, Edgecount=Edgecount_gen)
#   }
#   if(b==3){
#     genscores<-cbind(Edgecount=Edgecount_gen, C1c2=C1c2_gen, C1_l=C1_l_gen, Dmax=Dmax_gen )
#   }
#   if(b==4){
#     genscores<-cbind(Size=Size_gen, Edgecount=Edgecount_gen,Dmax=Dmax_gen  )#, C1c2=C1c2_gen, starSize=starSize_gen, 
#   }
#   if(b==5){
#     genscores<-cbind(Size=Size_gen, Edgecount=Edgecount_gen, Dmax=Dmax_gen  ) #, C1c2=C1c2_gen, C1=C1_gen
#   }
#   genscores<-cbind(genscores, class)
#   
#   class<-rep(0, times=length(Size_imp))
#   if(b==1){
#     impscores<-cbind(Size=Size_imp, Dmax=Dmax_imp, C1_l=C1_l_imp, starSize=starSize_imp, C1c2=C1c2_imp)
#   }  
#   if(b==2){
#     impscores<-cbind(Size=Size_imp, Dmax=Dmax_imp, C1_l=C1_l_imp, Edgecount=Edgecount_imp)
#   }
#   if(b==3){
#     impscores<-cbind(Edgecount=Edgecount_imp, C1c2=C1c2_imp, C1_l=C1_l_imp, Dmax=Dmax_imp )
#   }
#   if(b==4){
#     impscores<-cbind(Size=Size_imp, Edgecount=Edgecount_imp, Dmax=Dmax_imp  )#, C1c2=C1c2_imp, starSize=starSize_imp
#   }
#   if(b==5){
#     impscores<-cbind(Size=Size_imp, Edgecount=Edgecount_imp, Dmax=Dmax_imp  ) #, C1c2=C1c2_imp, C1=C1_imp
#   }
#   impscores<-cbind(impscores, class)
#   
#   
#   # Size_gen, Edgecount_gen, Edgenoderatio_gen, C1_gen, C2_gen, Isolated_gen, C1_l_gen, C1c2_gen, ENRc1_gen, Dmax_gen
#   # starSize_gen, twostarSize_gen
#   
#   
#   #we need to get these onto the genuine side soemhow using other measures
#   genIndToCheck<-which(Edgecount_gen>min(Edgecount_imp))
#   impIndToCheck<-which(Edgecount_imp<max(Edgecount_gen))
#   y_gen<-twostarSize_gen
#   y_imp<-twostarSize_imp
#   windows()
#   plot(x=Edgecount_imp[impIndToCheck], y=y_imp[impIndToCheck], pch=1, col=2, xlim=c(0,1), ylim=c(0,1))
#   points(x=Edgecount_gen[genIndToCheck], y=y_gen[genIndToCheck], pch=15, col=3)
# #look at the graphs compared
#   for(i in genIndToCheck[2:10]){
#     g<-genMCS[[i]]          
#     i1<-genComp[i,1]
#     i2<-genComp[i,2]
#     g1<-graphList$graph[[i1]]
#     g2<-graphList$graph[[i2]]
#     windows()
#     par(mfrow=c(2,2))
#     plotspatialgraph(g)
#     title(paste("gen MCS#", i, sep=""))
#     plotspatialgraph(g1)
#     title(paste("graph #", i1, sep=""))
#     plotspatialgraph(g2)
#     title(paste("graph #", i2, sep=""))
#     
#   }
#   for(i in impIndToCheck[1:5]){
#     g<-impMCS[[i]]          
#     i1<-impComp[i,1]
#     i2<-impComp[i,2]
#     g1<-graphList$graph[[i1]]
#     g2<-graphList$graph[[i2]]
#     windows()
#     par(mfrow=c(2,2))
#     plotspatialgraph(g)
#     title(paste("imp MCS#", i, sep=""))
#     plotspatialgraph(g1)
#     title(paste("graph #", i1, sep=""))
#     plotspatialgraph(g2)
#     title(paste("graph #", i2, sep=""))
#     
#   }
#   
# 
#   #############################################################
#   ### Let us use 3 (col=green) for genuine and 2(col=red) for imposter
#   
#   
#   th_Dmax<-3
#   th_C1c2<- min(C1c2_imp) #if scores is below this it is genuine
#   th_C1_l<-min(C1_l_imp)
#   th_starSize<-min(starSize_imp)
#   th_Size<-min(Size_imp)
#   inputscores<-numeric()
#   FMR<-numeric()
#   FNMR<-numeric()
#   GAR<-numeric()
#   GRR<-numeric()
#   for(s in 1:2){
#     if(s==1) inputscores<-genscores
#     if(s==2) inputscores<-impscores
#     for(i in 1:dim(inputscores)[1]){
#       scores<-inputscores[i,]
#       if(as.numeric(scores["Dmax"])<th_Dmax ) {scores["class"]<-2}
#       if(as.numeric(scores["Dmax"])>th_Dmax ) {scores["class"]<-3}
#       if(as.numeric(scores["Dmax"])==th_Dmax ) {
#         if(as.numeric(scores["C1c2"])<th_C1c2) scores["class"]<-3
#         if(as.numeric(scores["C1_l"])<th_C1_l) scores["class"]<-3
#         if(as.numeric(scores["starSize"])<th_starSize) scores["class"]<-3
#         if(as.numeric(scores["Size"])<th_Size) scores["class"]<-3 
#         if(as.numeric(scores["Size"])>=th_Size) scores["class"]<-2 
#       }
#       
#       inputscores[i,]<-scores
#     }#end of i loop
#     if(s==1){
#       FNMR<-(length(which(inputscores[,6]==2) )/length(inputscores[,1]))*100 #falsesly calling it imposter
#       GAR<-(length(which(inputscores[,6]==3) )/length(inputscores[,1]))*100 
#     }
#     if(s==2){
#       FMR<-(length(which(inputscores[,6]==3) )/length(inputscores[,1]))*100 #falsesly calling it genuine
#       GRR<-(length(which(inputscores[,6]==2) )/length(inputscores[,1]))*100 
#     }
#     
#   }
#  print(c(FMR, FNMR))
  ############################################################################
  
  gen_data_test<-cbind( Size_gen,  Edgecount_gen, Edgenoderatio_gen, C1_gen,  C2_gen,  Isolated_gen, C1_l_gen, C1c2_gen,  ENRc1_gen,  Dmax_gen, starSize_gen, twostarSize_gen, SpatialSize_gen)
  #   #pairs(gen_data)
  imp_data_test<-cbind( Size_imp,  Edgecount_imp, Edgenoderatio_imp, C1_imp,  C2_imp,  Isolated_imp, C1_l_imp, C1c2_imp,  ENRc1_imp,  Dmax_imp, starSize_imp, twostarSize_imp, SpatialSize_imp)
  class<-rep(1, times=length(Size_gen))
  genscores_all_test<-cbind(gen_data_test, class)
  class<-rep(2, times=length(Size_imp))
  impscores_all_test<-cbind(imp_data_test, class)
  
  
  
#   genscores_all_test<-cbind(Size=Size_gen, Edgecount=Edgecount_gen, Dmax=Dmax_gen, C1c2=C1c2_gen, C1_l=C1_l_gen, starSize=starSize_gen)
#   impscores_all_test<-cbind(Size=Size_imp, Edgecount=Edgecount_imp, Dmax=Dmax_imp, C1c2=C1c2_imp, C1_l=C1_l_imp, starSize=starSize_imp)
#   class<-rep(1, times=length(Size_gen))
#   genscores_all_test<-cbind(genscores_all_test, class)
#   class<-rep(2, times=length(Size_imp))
#   impscores_all_test<-cbind(impscores_all_test, class)
  
####plot the distribution of genuine and imposter points in training and testing data
  
#   windows()
#   par(mfrow=c(2,2))
#   plot(genscores_all[,1], genscores_all[,2], col=genscores_all[,10], pch=15, xlim=c(0,1), ylim=c(0,1), main="train gen" )
#   legend(0.8, 0.5, legend=c(1,2,3,4), col=c(1,2,3,4) , pch=15)
#   plot(impscores_all[,1], impscores_all[,2], col=impscores_all[,10], pch=16,, xlim=c(0,1), ylim=c(0,1), main="train imp")
#   legend(0.8, 0.5, legend=c(1,2,3,4), col=c(1,2,3,4) , pch=16)
#   plot(genscores_all_test[,1], genscores_all_test[,2], col=genscores_all_test[,10], pch=15, xlim=c(0,1), ylim=c(0,1), main="test gen" )
#   legend(0.8, 0.5, legend=c(1,2,3,4), col=c(1,2,3,4) , pch=15)
#   plot(impscores_all_test[,1], impscores_all_test[,2], col=impscores_all_test[,10], pch=16, xlim=c(0,1), ylim=c(0,1), main="test imp")
#   legend(0.8, 0.5, legend=c(1,2,3,4), col=c(1,2,3,4) , pch=16)

  trainData<-rbind(genscores_all, impscores_all)
  testData<-rbind(genscores_all_test, impscores_all_test)
  
#   #for the last column imposters that are above 20, just cap it off there
#   trainData[which(trainData[,13]>20),13]<-20
#   testData[which(testData[,13]>20),13]<-20
#   
  testData_test<-testData[,-dim(testData)[2]]
  
#   windows()
#   par(mfrow=c(2,1))
#   plot(trainData[,1], trainData[,2], col=trainData[,10], pch=trainData[,14], main="train")
#   plot(testData[,1], testData[,2], col=testData[,10], pch=testData[,14], main="test")
  
  trainData<-as.data.frame(trainData)
  testData<-as.data.frame(testData)
  testData_test<-as.data.frame(testData_test)
  
  require(MASS)
  priorProb<-c(dim(genscores_all)[1]/sum(dim(genscores_all)[1],dim(impscores_all)[1] ), dim(impscores_all)[1]/sum(dim(genscores_all)[1],dim(impscores_all)[1] ) )
  #z<-lda(class ~., trainData, prior=c(1,1)/2  )
  z<-lda(class ~., trainData, prior=priorProb  )
  train_lda<-predict(z, trainData)
  train_gen<- train_lda$x[1:dim(genscores_all)[1]]
  train_imp<-train_lda$x[-(1:dim(genscores_all)[1])]
  
  
  ###### Let us plot ROCS of all the training data and get the EER threshold value. #####
  ###### Then get the FMR and FNMR of the test set ######################################
  ####################################################################################################################
#   quartz()
#   plot(1:length(train_gen), train_gen, xlim=c(0,length(train_imp)), pch=1, col="green", ylim=c(min(train_gen), max(train_imp)) , type='p')
#   points(1:length(train_imp), train_imp, pch=2, col="red")
  
  roc_lda_train<-getROC_unnorm_nomx(train_gen, train_imp,0)
  
  roc_Size_train<-getROC_unnorm_nomx(genscores_all[,1], impscores_all[,1],0)
  roc_EdgeCount_train<-getROC_unnorm_nomx(genscores_all[,2], impscores_all[,2],0)
  roc_Edgenoderatio_train<-getROC_sim_nomx(genscores_all[,3], impscores_all[,3], 0)
  roc_C1_train<-getROC_unnorm_nomx(genscores_all[,4], impscores_all[,4], 0)
  roc_C2_train<-getROC_unnorm_nomx(genscores_all[,5], impscores_all[,5], 0)
  roc_Isolated_train<-getROC_unnorm_nomx(genscores_all[,6], impscores_all[,6], 0)
  roc_C1_l_train<-getROC_unnorm_nomx(genscores_all[,7], impscores_all[,7], 0)
  roc_C1c2_train<-getROC_unnorm_nomx(genscores_all[,8], impscores_all[,8], 0)
  roc_ENRc1_train<-getROC_sim_nomx(genscores_all[,9], impscores_all[,9], 0)
  roc_Dmax_train<-getROC_sim_nomx(genscores_all[,10], impscores_all[,10], 0)
  roc_starSize_train<-getROC_unnorm_nomx(genscores_all[,11], impscores_all[,11], 0)
  roc_twostarSize_train<-getROC_unnorm_nomx(genscores_all[,12], impscores_all[,12], 0)
  roc_SpatialSize_train<-getROC_sim_nomx(genscores_all[,13], impscores_all[,13], 0)
  
  
  roc_list<-c( round(roc_Size_train$eer,4), round(roc_EdgeCount_train$eer,4), round(roc_Edgenoderatio_train$eer,4), round(roc_C1_train$eer,4), round(roc_C2_train$eer,4), round(roc_Isolated_train$eer,4), round(roc_C1_l_train$eer,4), round(roc_C1c2_train$eer,4), round(roc_ENRc1_train$eer,4), round(roc_Dmax_train$eer,4), round(roc_starSize_train$eer,4), round(roc_twostarSize_train$eer,4) , round(roc_SpatialSize_train$eer,4)  )
  bestEERInd<-order(roc_list, decreasing=FALSE)
  
  dbType<-""
  if(b==1) dbType<-"Retina"
  if(b==2) dbType<-"Palmvein"
  if(b==3) dbType<-"Wristvein"
  if(b==4) dbType<-"Handvein - SFIR"
  if(b==5) dbType<-"Handvein - SNIR"
  titleText<-paste("Training Set DET curves for ", dbType, sep = "")
  
  EER_list<-c(round(roc_lda_train$eer,4), roc_list[bestEERInd[1:3]])
  feature_list<-c("LDA", "", "", "")
  legendtext1<-paste(feature_list[1], " EER = ", EER_list[1], sep="")
  for(ind in 1:3){
    i<-bestEERInd[ind]
    if(i==1) feature_list[ind+1]<-"d_v"
    if(i==2) feature_list[ind+1]<-"d_e"
    if(i==3) feature_list[ind+1]<-"rho"
    if(i==4) feature_list[ind+1]<-"d_C1"
    if(i==5) feature_list[ind+1]<-"d_C2"
    if(i==6) feature_list[ind+1]<-"d_isolated"
    if(i==7) feature_list[ind+1]<-"d_C1C2"
    if(i==8) feature_list[ind+1]<-"d_rho_C1"
    if(i==9) feature_list[ind+1]<-"d_max"
    if(i==10) feature_list[ind+1]<-"d_S"
    if(i==11) feature_list[ind+1]<-"d_S^2"
    if(i==12) feature_list[ind+1]<-"d_a"
    
  }
  legendtext1<-paste(feature_list, "EER = ", EER_list,sep=" ")
  #feature_list1<-expression("LDA ", d[e], d[v], d[S^2]) #HV
  #feature_list1<-expression("LDA ", d[e], d[v], d[rho[C1]]) #PV, WV
  #feature_list1<-expression("LDA ", d[v], d[e], d[S^2]) #HV-SFIR
  #feature_list1<-expression("LDA ", d[v], d[e], d[rho[C1]]) #HV-SNIR
  feature_list1<-expression("LDA ", d[v], d[e], d[C1C2]) #HV-SNIR
  
  quartz()
  plot(roc_lda_train$x, roc_lda_train$y, log="xy", col=1, lwd=2, type='l', lty=1, main=titleText, xlab="FMR (log scale)", ylab="FNMR (log scale)")
 #, ylim=c(0,1), xlim=c(0,1)
  for(ind in 1:3){
   i<-bestEERInd[ind]
   roc_train<-getROC_unnorm_nomx(genscores_all[,i], impscores_all[,i], 0)
   lines(roc_train$x, roc_train$y, col=ind+1, lwd=2, lty=ind+1)
 }
    
  legend(0.3,1, legend = feature_list1, col=1:4, lwd=3, lty=1:4)
  
   # get the three thresholds for the best performing feature EER, FMR100 and FMR1000.
    bestIndex<-bestEERInd[1]
    roc_best<-getROC_unnorm_nomx(genscores_all[,bestIndex], impscores_all[,bestIndex],0)
    Thresholds<-c(0,0,0) # EER, FMR100 and FMR1000 thresholds
    FNMR_tr<-c(0,0,0)
    Thresholds[1]<-roc_best$th
    FNMR_tr[1]<-roc_best$eer
    # Finding the threshold for FMR100
    z1=abs(0.01-roc_best$x)
    FMR100_index<-order(z1)[1]
    roc_best$y[FMR100_index]
    Thresholds[2]<-roc_best$w[FMR100_index]
    FNMR_tr[2]<-roc_best$y[FMR100_index]
    
    # Finding the threshold for FMR1000
    z1=abs(0.001-roc_best$x)
    FMR1000_index<-order(z1)[1]
    roc_best$y[FMR1000_index]
    Thresholds[3]<-roc_best$w[FMR1000_index]
    FNMR_tr[3]<-roc_best$y[FMR1000_index]
    
    print(Thresholds)
    print(FNMR_tr)
    
    

#   legendtext1<-c(paste("LDA scores, EER=", round(roc_lda_train$eer,4), sep=""),paste("Size, EER=", round(roc_Size_train$eer,4), sep=""),paste("Edgecount, EER=", round(roc_EdgeCount_train$eer,4), sep=""),paste("Edgenoderatio, EER=", round(roc_Edgenoderatio_train$eer,4), sep=""),paste("C1, EER=", round(roc_C1_train$eer,4), sep=""),paste("C2, EER=", round(roc_C2_train$eer,4), sep=""))
#   legendtext2<-c(paste("Isolated, EER=", round(roc_Isolated_train$eer,4), sep=""),paste("C1_l, EER=", round(roc_C1_l_train$eer,4), sep=""),paste("C1c2, EER=", round(roc_C1c2_train$eer,4), sep=""),paste("ENRc1, EER=", round(roc_ENRc1_train$eer,4), sep=""),paste("Dmax, EER=", round(roc_Dmax_train$eer,4), sep=""),paste("starSize, EER=", round(roc_starSize_train$eer,4), sep=""),paste("twostarSize, EER=", round(roc_twostarSize_train$eer,4), sep=""))
#   windows()
#   par(mfrow=c(2,1))
#   
#   plot(roc_lda_train$x, roc_lda_train$y, col=1, lwd=2, type='l', lty=1, ylim=c(0,1), xlim=c(0,1), main=titleText)
#   lines(roc_Size_train$x, roc_Size_train$y, col=2, lwd=2, lty=2)
#   lines(roc_EdgeCount_train$x, roc_EdgeCount_train$y, col=3, lty=3)
#   lines(roc_Edgenoderatio_train$x, roc_Edgenoderatio_train$y, col=4, lwd=2, lty=4)
#   lines(roc_C1_train$x, roc_C1_train$y, col=5, lwd=2, lty=5)
#   lines(roc_C2_train$x, roc_C2_train$y, col=6, lwd=2, lty=6)
#   legend(0.6,1, legendtext1, col=1:6, lwd=3, lty=1:6)
#   
#   plot(roc_Isolated_train$x, roc_Isolated_train$y, col=7, lwd=2, lty=7, type='l', ylim=c(0,1), xlim=c(0,1), main=titleText)
#   lines(roc_C1_l_train$x, roc_C1_l_train$y, col=8, lwd=2,  lty=8)
#   lines(roc_C1c2_train$x, roc_C1c2_train$y, col=9, lwd=2,  lty=9)
#   lines(roc_ENRc1_train$x, roc_ENRc1_train$y, col=10, lwd=2,  lty=10)
#   lines(roc_Dmax_train$x, roc_Dmax_train$y, col=11, lwd=2,  lty=11)
#   lines(roc_starSize_train$x, roc_starSize_train$y, col=12, lwd=2,  lty=12)
#   lines(roc_twostarSize_train$x, roc_twostarSize_train$y, col=13, lwd=2,  lty=13)    
#   legend(0.6,1, legendtext2, col=7:13, lwd=3, lty=7:13)
#   
 
  #############################################
  ### get a list of thresholds for EER ########
  #############################################
  
  training_Th<-c(round(roc_lda_train$th,4), round(roc_Size_train$th,4), round(roc_EdgeCount_train$th,4), round(roc_Edgenoderatio_train$th,4), round(roc_C1_train$th,4), round(roc_C2_train$th,4), round(roc_Isolated_train$th,4), round(roc_C1_l_train$th,4), round(roc_C1c2_train$th,4), round(roc_ENRc1_train$th,4), round(roc_Dmax_train$th,4), round(roc_starSize_train$th,4), round(roc_twostarSize_train$th,4), round(roc_SpatialSize_train$th,4)   )

  
  
  
  
  ########### generate the lda test data
  
  Test_lda<-predict(z, testData_test)
  
 op<-cbind(Test_lda$class,testData[,14])
op[which(op[,1]!=op[,2]),]
  head(Test_lda$posterior)
  
  test_gen<- Test_lda$x[1:dim(genscores_all_test)[1]]
  test_imp<-Test_lda$x[-(1:dim(genscores_all_test)[1])]
#  windows()
#   plot(1:length(test_gen), test_gen, xlim=c(0,length(test_imp)), pch=1, col="green", ylim=c(min(test_gen), max(test_imp)) , type='p')
#   points(1:length(test_imp), test_imp, pch=2, col="red")
  
 
  ############# Plot the testdata and show the threshold ################################
#   windows()
#   par(mfrow=c(2,2), oma=c(0,0,2,0))
#   for(i in c(1,(bestEERInd[1:3]+1) )){
#     genData<-numeric()
#     impData<-numeric()
#     titlestring<-""
#    if(i==1){
#      genData<-test_gen
#      impData<-test_imp
#      titlestring<-"Test scores from LDA model output"
#    }
#     if(i>1){
#       genData<-genscores_all_test[,i-1]
#       impData<-impscores_all_test[,i-1]
#       titlestring<-paste("Test scores from ", colnames(genscores_all_test)[i-1], sep="" )
#     }
#      
#    
#     plot(1:length(genData), genData, xlim=c(0,length(impData)), pch=1, col="green", ylim=c(min(genData), max(impData)) , type='p', main=titlestring)
#     points(1:length(impData), impData, pch=2, col="red")
#     abline(h=training_Th[i], col="blue")
#     
#   }
#   
#   
  
#   # calculate the FMR and FNMR at the training threshold for each feature
#   
#   FMR_list<-rep(0, times=14)
#   FNMR_list<-rep(0, times=14)
#   
#   for(i in c(1,(bestEERInd[1:3]+1) )){
#     genData<-numeric()
#     impData<-numeric()
#     index<-0
#     if(i==1){
#       genData<-test_gen
#       impData<-test_imp
#      
#     }
#     if(i>1){
#       genData<-genscores_all_test[,i-1]
#       impData<-impscores_all_test[,i-1]
#       
#     }
#     
#     FMR_list[i]<-length(which(genData>=training_Th[i])) / length(genData)  
#     FNMR_list[i]<- length(which(impData<training_Th[i])) / length(impData)  
#     
#   }
#   cbind(FMR_list,FNMR_list)
  
  roc_lda<-getROC_unnorm_nomx(test_gen, test_imp,0)
  #roc_lda<-getROC_unnorm_nomx(train_gen, train_imp,0)
  roc_Size_test<-getROC_unnorm_nomx(genscores_all_test[,1], impscores_all_test[,1],0)
  roc_EdgeCount_test<-getROC_unnorm_nomx(genscores_all_test[,2], impscores_all_test[,2],0)
  roc_Edgenoderatio_test<-getROC_sim_nomx(genscores_all_test[,3], impscores_all_test[,3], 0)
  roc_C1_test<-getROC_unnorm_nomx(genscores_all_test[,4], impscores_all_test[,4], 0)
  roc_C2_test<-getROC_unnorm_nomx(genscores_all_test[,5], impscores_all_test[,5], 0)
  roc_Isolated_test<-getROC_unnorm_nomx(genscores_all_test[,6], impscores_all_test[,6], 0)
  roc_C1_l_test<-getROC_unnorm_nomx(genscores_all_test[,7], impscores_all_test[,7], 0)
  roc_C1c2_test<-getROC_unnorm_nomx(genscores_all_test[,8], impscores_all_test[,8], 0)
  roc_ENRc1_test<-getROC_sim_nomx(genscores_all_test[,9], impscores_all_test[,9], 0)
  roc_Dmax_test<-getROC_sim_nomx(genscores_all_test[,10], impscores_all_test[,10], 0)
  roc_starSize_test<-getROC_unnorm_nomx(genscores_all_test[,11], impscores_all_test[,11], 0)
  roc_twostarSize_test<-getROC_unnorm_nomx(genscores_all_test[,12], impscores_all_test[,12], 0)
  roc_SpatialSize_test<-getROC_sim_nomx(genscores_all_test[,13], impscores_all_test[,13], 0)
  
  titleText<-"ROCs for graph features and LDA measure"
  
  #####get FMR and FNMRs for test data
  FMR_test<-c(0,0,0)
  FNMR_test<-c(0,0,0)
  roc_test<-getROC_unnorm_nomx(genscores_all_test[,bestIndex], impscores_all_test[,bestIndex],0)
  #get at EERTh
  z2<-abs(Thresholds[1]-roc_test$w)
  EER_index<-order(z2)[1]
  FMR_test[1]<-roc_test$x[EER_index]
  FNMR_test[1]<-roc_test$y[EER_index]
  
  #get at FMR100 Th
  z2<-abs(Thresholds[2]-roc_test$w)
  FMR100_index<-order(z2)[1]
  FMR_test[2]<-roc_test$x[FMR100_index]
  FNMR_test[2]<-roc_test$y[FMR100_index]
  
  #get at FMR1000 Th
  z2<-abs(Thresholds[3]-roc_test$w)
  FMR1000_index<-order(z2)[1]
  FMR_test[3]<-roc_test$x[FMR1000_index]
  FNMR_test[3]<-roc_test$y[FMR1000_index]
  
  print(FMR_test)
  print(FNMR_test)
  
#   legendtext1<-c(paste("LDA scores, EER=", round(roc_lda$eer,4), sep=""),paste("Size, EER=", round(roc_Size$eer,4), sep=""),paste("Edgecount, EER=", round(roc_EdgeCount$eer,4), sep=""),paste("Edgenoderatio, EER=", round(roc_Edgenoderatio$eer,4), sep=""),paste("C1, EER=", round(roc_C1$eer,4), sep=""),paste("C2, EER=", round(roc_C2$eer,4), sep=""))
#   legendtext2<-c(paste("Isolated, EER=", round(roc_Isolated$eer,4), sep=""),paste("C1_l, EER=", round(roc_C1_l$eer,4), sep=""),paste("C1c2, EER=", round(roc_C1c2$eer,4), sep=""),paste("ENRc1, EER=", round(roc_ENRc1$eer,4), sep=""),paste("Dmax, EER=", round(roc_Dmax$eer,4), sep=""),paste("starSize, EER=", round(roc_starSize$eer,4), sep=""),paste("twostarSize, EER=", round(roc_twostarSize$eer,4), sep=""))
#   windows()
#   par(mfrow=c(2,1))
#   
#   plot(roc_lda$x, roc_lda$y, col=1, lwd=2, type='l', lty=1, ylim=c(0,1), xlim=c(0,1), main=titleText)
#   lines(roc_Size$x, roc_Size$y, col=2, lwd=2, lty=2)
#   lines(roc_EdgeCount$x, roc_EdgeCount$y, col=3, lty=3)
#   lines(roc_Edgenoderatio$x, roc_Edgenoderatio$y, col=4, lwd=2, lty=4)
#   lines(roc_C1$x, roc_C1$y, col=5, lwd=2, lty=5)
#   lines(roc_C2$x, roc_C2$y, col=6, lwd=2, lty=6)
#   legend(0.6,1, legendtext1, col=1:6, lwd=3, lty=1:6)
#   
#   plot(roc_Isolated$x, roc_Isolated$y, col=7, lwd=2, lty=7, type='l', ylim=c(0,1), xlim=c(0,1), main=titleText)
#   lines(roc_C1_l$x, roc_C1_l$y, col=8, lwd=2,  lty=8)
#   lines(roc_C1c2$x, roc_C1c2$y, col=9, lwd=2,  lty=9)
#   lines(roc_ENRc1$x, roc_ENRc1$y, col=10, lwd=2,  lty=10)
#   lines(roc_Dmax$x, roc_Dmax$y, col=11, lwd=2,  lty=11)
#   lines(roc_starSize$x, roc_starSize$y, col=12, lwd=2,  lty=12)
#   lines(roc_twostarSize$x, roc_twostarSize$y, col=13, lwd=2,  lty=13)    
#   legend(0.6,1, legendtext2, col=7:13, lwd=3, lty=7:13)
  
  roc_list<-c( round(roc_Size$eer,4), round(roc_EdgeCount$eer,4), round(roc_Edgenoderatio$eer,4), round(roc_C1$eer,4), round(roc_C2$eer,4), round(roc_Isolated$eer,4), round(roc_C1_l$eer,4), round(roc_C1c2$eer,4), round(roc_ENRc1$eer,4), round(roc_Dmax$eer,4), round(roc_starSize$eer,4), round(roc_twostarSize$eer,4)   )
  order(roc_list, decreasing=FALSE)
  
  
  ############################################################################################
  ########## TRY some thing else #############################################################
  ########## Tune using N-1 sample points and test the last one. #############################
  
  
  
  
  
  
  
#   ##########################
  #svm_range<-1:12
  svm_range<-order(roc_list)[1:2]
#   genscores_train<-genscores_all[, c(svm_range,dim(genscores_all)[2])]
#   impscores_train<-impscores_all[,c(svm_range,dim(impscores_all)[2])]
  
  genscores_train<-cbind( train_gen, genscores_all[, c(svm_range,dim(genscores_all)[2])] )
  impscores_train<-cbind(train_imp, impscores_all[,c(svm_range,dim(impscores_all)[2])] )
  ###Now build an svm classifier
  genscores_train[,dim(genscores_train)[2]]<-1
  impscores_train[,dim(impscores_train)[2]]<-2
  
  gen_tune<-sample(1:dim(genscores_train)[1], size=10)  
  imp_tune<-sample(1:dim(impscores_train)[1], size=10)  
  scores_tune<-as.data.frame(rbind(genscores_train[gen_tune,], impscores_train[imp_tune,]))
  #   model00<-tune.svm(class~.,data=scores_tune, cost=10^(-2:2), kernel="linear")
  #   summary(model00)
  #   
  gen_train<-1:dim(genscores_train)[1]  
 # imp_train<-sample(1:dim(impscores_train)[1], size=dim(genscores_train)[1], replace=FALSE)
  #instead of random selection, take the worst performing imposter scores for nodescore.
  imp_train<-order(impscores_train[,1])[1:dim(genscores_train)[1]]
  scores_train<-as.data.frame(rbind(genscores_train[gen_train,], impscores_train[imp_train,]))
  
  
 # genscores_test<-genscores_all_test[, c(svm_range,dim(genscores_all_test)[2])]
#  impscores_test<-impscores_all_test[,c(svm_range,dim(impscores_all_test)[2])]
  
  genscores_test<-cbind( test_gen,genscores_all_test[, c(svm_range,dim(genscores_all_test)[2])])
  impscores_test<-cbind( test_imp,impscores_all_test[,c(svm_range,dim(impscores_all_test)[2])])
  #   ###Now build an svm classifier
  genscores_test[,dim(genscores_test)[2]]<-1
  impscores_test[,dim(impscores_test)[2]]<-2
  #   gen_tune<-sample(1:dim(genscores_test)[1], size=10)  
  #   imp_tune<-sample(1:dim(impscores_test)[1], size=10)  
  #   scores_tune<-as.data.frame(rbind(genscores_test[gen_tune,], impscores_test[imp_tune,]))
  #   model00<-tune.svm(class~.,data=scores_tune, cost=10^(-2:2), kernel="linear")
  #   summary(model00)
  
  scores_test<-rbind(genscores_test, impscores_test)
  scores_test<-as.data.frame(scores_test)
  
  trainSize<-dim(scores_train)[1]
  testSize<-dim(scores_test)[1]
  
  colnames(scores_train)<-colnames(scores_test)
  scores<-rbind(scores_train,scores_test)
  op_class<-scores$class
  trainIndex<-1:trainSize
  testIndex<-(trainSize+1):(trainSize+testSize)
  
  # Tune svm model
    model00<-tune.svm(class~.,data=scores_tune, cost=10^(-2:2), gamma=10^(-4:-1), kernel="linear")
    summary(model00)
  plot(model00)

   # Build svm model
    # model<-svm(class~.,data=scores[testIndex[c(gen_tune, imp_tune)],], cost=as.numeric(model00[[1]]),kernel="linear", type="C-classification" )
    model<-svm(class~.,data=scores[trainIndex,], gamma=as.numeric(model00[[1]][1]) , cost=as.numeric(model00[[1]][2]),kernel="linear", type="C-classification" )
  
  summary(model)

  y<-scores[testIndex,-dim(scores)[2]]
  
  pred<-predict(model, y )
   table(pred, op_class[testIndex])
  
  
  
  # TEST SVM MODEL
   
  # compute decision values and probabilities:
  pred<-predict(model, y, decision.values = TRUE)
 
  svmScores<-attr(pred, "decision.values")
    
  svmScores_gen<-svmScores[1:dim(genscores_test)[1]]
  svmScores_imp<-svmScores[-(1:dim(genscores_test)[1])]
  
#   windows()
#   plot(1:length( svmScores_gen),  svmScores_gen, pch=15, col="green", xlim=c(0, max(length( svmScores_gen), length( svmScores_imp))), ylim=c(-2, 2) )
#   points(1:length( svmScores_imp),  svmScores_imp, pch=16, col="red")
  
  #getROC_unnorm_nomx
  roc_svm<-getROC_sim_nomx( svmScores_gen, svmScores_imp,0)
  roc_Size<-getROC_unnorm_nomx(Size_gen, Size_imp,0)
  roc_EdgeCount<-getROC_unnorm_nomx(Edgecount_gen, Edgecount_imp,0)
  
  legendtext<-c(paste("SVM scores, EER=", round(roc_svm$eer,4), sep=""),paste("Size, EER=", round(roc_Size$eer,4), sep=""),paste("Edgecount, EER=", round(roc_EdgeCount$eer,4), sep=""))

  windows()
  plot(roc_svm$x, roc_svm$y, col=1, lwd=2, type='l', lty=1, ylim=c(0,1), xlim=c(0,1))
  points(roc_Size$x, roc_Size$y, col=2, lwd=2, type='l', lty=2)
  points(roc_EdgeCount$x, roc_EdgeCount$y, col=3, lwd=2, type='l', lty=3)
  legend(0.5,1, legendtext, col=1:3, lwd=2, lty=1:3)
#   
#   #plot scores
#   minimum<-min(range(svmScores_gen)[1], range(svmScores_imp)[1] )
#   maximum<-max(range(svmScores_gen)[2], range(svmScores_imp)[2] )
#   windows()
#   plot(1:length(svmScores_gen), svmScores_gen, col="green", pch=1, ylim=c(minimum, maximum), xlim=c(1, max(length(svmScores_gen), length(svmScores_imp))))
#   points(1:length(svmScores_imp), svmScores_imp, col="red", pch=2)
#   
#   minimum<-min(range(Size_gen)[1], range(Size_imp)[1] )
#   maximum<-max(range(Size_gen)[2], range(Size_imp)[2] )
#   windows()
#   plot(1:length(Size_gen), Size_gen, col="green", pch=1, ylim=c(minimum, maximum), xlim=c(1, max(length(Size_gen), length(Size_imp))))
#   points(1:length(Size_imp), Size_imp, col="red", pch=2)
#   
#   minimum<-min(range(Edgecount_gen)[1], range(Edgecount_imp)[1] )
#   maximum<-max(range(Edgecount_gen)[2], range(Edgecount_imp)[2] )
#   windows()
#   plot(1:length(Edgecount_gen), Edgecount_gen, col="green", pch=1, ylim=c(minimum, maximum), xlim=c(1, max(length(Edgecount_gen), length(Edgecount_imp))))
#   points(1:length(Edgecount_imp), Edgecount_imp, col="red", pch=2)
  
}#end of b loop

