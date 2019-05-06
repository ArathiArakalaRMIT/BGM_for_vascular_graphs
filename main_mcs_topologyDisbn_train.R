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
opwd<-"H:/Uni Stuff/Uni Stuff/2015/Research Projects/ARC Discovery/Outputs/VascularGraphs/bestCid/"
ipwd<-"H:/Uni Stuff/Uni Stuff/2015/Research Projects/ARC Discovery/Outputs"
bestStruct_testdata<-c("twostar", "edge", "edge", "edge", "edge", "twostar", "edge")
mcsStruct_testdata<-c("node", "node", "node", "node", "node")
cid_n_best<-c(3,5,3,3,3)
cid_e_best<-c(9,3,7,7,5)


for(b in 1:5){
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
  cid_n<-cid_n_best[b]
  cid_e<-cid_e_best[b]
  mcsStruct<-mcsStruct_testdata[b]
  
  #read all the files and their cid_n and cid_e
  outputfile<-paste(ipwd,"/",folder_list[b],"/",sourceid_list[b],"/BetterGraphMatching/genMCS_",mcsStruct,"Based_train_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tests on training data
  if(file.exists(outputfile)==FALSE) next
  load(outputfile)
  genMCS<-genMCS_nodeBased
  
  outputfile<-paste(ipwd,"/",folder_list[b],"/",sourceid_list[b],"/BetterGraphMatching/impMCS_",mcsStruct,"Based_train_",cid_n,"_",cid_e,".Rdata", sep="") #cost function tests on training data
  if(file.exists(outputfile)==FALSE) next
  load(outputfile)
  impMCS<-impMCS_nodeBased
          
  
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
          
          starSize_gen[i]<-1- ( dim(get_star(g))[1]/sqrt(network.size(g1)*network.size(g2))  )
          twostarSize_gen[i]<-1- ( dim(get_twostar(g))[1]/sqrt(network.size(g1)*network.size(g2))   )
          
          degree_tmp<-rep(0, times=max(deg) )
          for(j in 1:max(deg)){
            degree_tmp[j]<-length(which(deg==j))
          }
          degree_tmp<-degree_tmp/sum(degree_tmp)
          degreedisbn_gen[i,1:max(deg)]<-degree_tmp
          
  

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
          
          starSize_imp[i]<-1- ( dim(get_star(g))[1]/sqrt(network.size(g1)*network.size(g2))  )
          twostarSize_imp[i]<-1- ( dim(get_twostar(g))[1]/sqrt(network.size(g1)*network.size(g2))  )
          
          degree_tmp<-rep(0, times=max(deg) )
          for(j in 1:max(deg)){
            degree_tmp[j]<-length(which(deg==j))
          }
          degree_tmp<-degree_tmp/sum(degree_tmp)
          degreedisbn_imp[i,1:max(deg)]<-degree_tmp
          
          
          
          
        }
  
 C2_imp[which(is.na(C2_imp)==TRUE)]<-1
 C1c2_imp[which(is.na(C1c2_imp)==TRUE)]<-1
  C2_gen[which(is.na(C2_gen)==TRUE)]<-1
  C1c2_gen[which(is.na(C1c2_gen)==TRUE)]<-1
#   
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
#   
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

  gen_data<-cbind( Size_gen,  Edgecount_gen, Edgenoderatio_gen, C1_gen,  C2_gen,  Isolated_gen, C1_l_gen, C1c2_gen,  ENRc1_gen,  Dmax_gen, starSize_gen, twostarSize_gen)
#pairs(gen_data)
  imp_data<-cbind( Size_imp,  Edgecount_imp, Edgenoderatio_imp, C1_imp,  C2_imp,  Isolated_imp, C1_l_imp, C1c2_imp,  ENRc1_imp,  Dmax_imp, starSize_imp, twostarSize_imp)
  
#   #1
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
  
  genscores_all<-cbind(Size=Size_gen, Edgecount=Edgecount_gen, Dmax=Dmax_gen, C1c2=C1c2_gen, C1_l=C1_l_gen, starSize=starSize_gen)
  impscores_all<-cbind(Size=Size_imp, Edgecount=Edgecount_imp, Dmax=Dmax_imp, C1c2=C1c2_imp, C1_l=C1_l_imp, starSize=starSize_imp)
  class<-rep(0, times=length(Size_gen))
  genscores_all<-cbind(genscores_all, class)
  class<-rep(0, times=length(Size_imp))
  impscores_all<-cbind(impscores_all, class)
  
  svm_range<-1:2
  genscores<-genscores_all[, c(svm_range,dim(genscores_all)[2])]
  impscores<-impscores_all[,c(svm_range,dim(impscores_all)[2])]
  ###Now build an svm classifier
  genscores[,dim(genscores)[2]]<-"G"
  impscores[,dim(impscores)[2]]<-"I"
  gen_tune<-sample(1:dim(genscores)[1], size=10)  
  imp_tune<-sample(1:dim(impscores)[1], size=10)  
  scores_tune<-as.data.frame(rbind(genscores[gen_tune,], impscores[imp_tune,]))
  model00<-tune.svm(class~.,data=scores_tune, cost=10^(-2:2), kernel="linear")
  summary(model00)
  
  scores<-rbind(genscores, impscores)
  scores<-as.data.frame(scores)
  
  
  
  
  model<-svm(class~.,data=scores, cost=as.numeric(model00[[1]]),kernel="linear" )
  summary(model)
  x<-subset(scores, select=-class)
  y<-subset(scores, select=class)
  pred<-as.vector( predict(model, x, decision.values = TRUE) )

  
  
  GAR<-length(which(pred[1:length(genscores[,1])]=="G"))/length(genscores[,1]) * 100
  FNMR<-100-GAR
  GRR<-length(which(pred[(length(genscores[,1])+1):(length(genscores[,1])+length(impscores[,1]) )]=="I"))/length(impscores[,1]) * 100
  FMR<-100-GRR
  
  print(c(FMR, FNMR))
  
  
  # compute decision values and probabilities:
  pred<-predict(model, x, decision.values = TRUE)
  svmScores<-attr(pred, "decision.values")
  svmScores_gen<-svmScores[1:dim(genscores)[1]]
  svmScores_imp<-svmScores[-(1:dim(genscores)[1])]
  
  roc_svm<-getROC_sim_nomx( svmScores_gen, svmScores_imp,0)
  roc_Size<-getROC_unnorm_nomx(Size_gen, Size_imp,0)
  roc_EdgeCount<-getROC_unnorm_nomx(Edgecount_gen, Edgecount_imp,0)
  
  legendtext<-c(paste("SVM scores, EER=", round(roc_svm$eer,4), sep=""),paste("Size, EER=", round(roc_Size$eer,4), sep=""),paste("Edgecount, EER=", round(roc_EdgeCount$eer,4), sep=""))

  windows()
  plot(roc_svm$x, roc_svm$y, col=1, lwd=2, type='l', lty=1, ylim=c(0,1), xlim=c(0,1))
  points(roc_Size$x, roc_Size$y, col=2, lwd=2, type='l', lty=2)
  points(roc_EdgeCount$x, roc_EdgeCount$y, col=3, lwd=2, type='l', lty=3)
  legend(0.5,1, legendtext, col=1:3, lwd=2, lty=1:3)
}#end of b loop


windows()
plot(1:length(svmScores_gen), svmScores_gen, col="green", pch=1, ylim=c(-1,1))
points(1:length(svmScores_imp), svmScores_imp, col="red", pch=2)