library(igraph)
library(ROntoTools)
library(graph)
library(e1071)
library(NetComp)


setwd("D:/BreastCancerResearch/DrugRepositioning/R/DNN")
source("NetUnion2.r")

getKEGG<-function(){
  kpg <- keggPathwayGraphs("hsa", relPercThresh=0, updateCache = FALSE, verbose = TRUE)
  kpg <- setEdgeWeights(kpg, edgeTypeAttr = "subtype", edgeWeightByType = list(activation = 1, inhibition = -1, expression = 1, repression = -1), defaultWeight = 0)
  save(kpg,file = "KEGG.RData")
}

CheckAllGenesInKegg<-function(kpg)
{
  L=NULL
  for (i in kpg) {
    M=data.frame(row.names(as(i,"matrix")))
    L=rbind(L,M)
  }
  
  L=unique(L)
  return(L)
}  

mergeMatrix <- function()
{
  load("KEGG.RData")
  
  #M=as(as_adj(igraph.from.graphNEL(kpg[[1]],weight = TRUE)),"matrix")
  M=as(kpg[[1]],"matrix")
  for (i in 2:50) {
    print(i)
    M2=as(kpg[[i]],"matrix")
    M=NetUnion2(M,M2)
    #M=as(NetComp::netUnion(M,M2),"matrix")
  }
  M_1=M
  
  M=as(kpg[[51]],"matrix")
  for (i in 52:100) {
    print(i)
    M2=as(kpg[[i]],"matrix")
    M=NetUnion2(M,M2)
    #M=as(NetComp::netUnion(M,M2),"matrix")
  }
  M_2=M
  
  M=as(kpg[[101]],"matrix")
  for (i in 102:150) {
    print(i)
    M2=as(kpg[[i]],"matrix")
    M=NetUnion2(M,M2)
    #M=as(NetComp::netUnion(M,M2),"matrix")
  }
  M_3=M
  
  M=as(kpg[[151]],"matrix")
  for (i in 152:200) {
    print(i)
    M2=as(kpg[[i]],"matrix")
    M=NetUnion2(M,M2)
    #M=as(NetComp::netUnion(M,M2),"matrix")
  }
  M_4=M
  
  M=as(kpg[[201]],"matrix")
  for (i in 202:250) {
    print(i)
    M2=as(kpg[[i]],"matrix")
    M=NetUnion2(M,M2)
    #M=as(NetComp::netUnion(M,M2),"matrix")
  }
  M_5=M
  
  M=as(kpg[[251]],"matrix")
  for (i in 252:length(kpg)) {
    print(i)
    M2=as(kpg[[i]],"matrix")
    M=NetUnion2(M,M2)
    #M=as(NetComp::netUnion(M,M2),"matrix")
  }
  M_6=M
  
  #----------------
  M=NetUnion2(M_1,M_2)
  M=NetUnion2(M,M_3)
  M=NetUnion2(M,M_4)
  M=NetUnion2(M,M_5)
  M=NetUnion2(M,M_6)
  
  return(M)
}

createKEGG_Matrix <-function(){
  M=mergeMatrix()
  #Genes=read.csv("Kegg_Genes.txt")
  #G=Genes$Gene.Names
  #row.names(M)<-colnames(M)<-G
  M_original=M
  
  for (i in 1:dim(M)[1])
    for (j in 1:dim(M)[1]){
      if(M[i,j]==0)
        M[i,j]=NA
      else
        M[i,j]=abs(M[i,j])
    }
  ShortestPath=allShortestPaths(M)
  
  save(file="KEGG_Matrix.RData",M,M_original,ShortestPath)
  
  
  # example -----
  which(row.names(M)=="TP53")
}

extendMatrix<-function(M){
  Updated=0
  for (i in 1:dim(M)[1]) {
    for (j in 1:dim(M)[1])
      if(M[i,j]!=0)
        for (k in 1:dim(M)[1]) 
          if(M[i,k]!=0 & M[j,k]!=0)
          {
            if(M[i,j]==1 & M[j,k]==1 & i!=k)
            {
              Updated=Updated+1
              M[i,k]=1
            }
            else if (M[i,j]==-1 & M[j,k]==-1 & i!=k)
            {
              Updated=Updated+1
              M[i,k]=-1
            }
          }
    print(paste0(i," -> ",Updated))
    
  }
  return(list(M,Updated))  
}
  
# Res=extendMatrix(M)
# M=Res[[1]]
# Updated=Res[[2]]
#
#start=20
#end=84
#P=extractPath(AS, start, end)

runAll<-function()
{
  #getKEGG()
  createKEGG_Matrix()
}




  