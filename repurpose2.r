rm(list=ls(all=TRUE))




#==========================
Cutted=1    # 1206  4946 5876
Subtype=1
Cutoff=8  # cutoff for drug expression genes (between 0 and 10) - the higher the value, the algorithm keeps less genes

BidirectionalPath=FALSE
VisualizeGraph=FALSE
#==========================

library(e1071)
library(igraph)
library(visNetwork)
library(igraph)

setwd("D:/BreastCancerResearch/DrugRepositioning/R/DNN/")

#load("DrugGene.RData")

source("NetUnion2.r")

#--------------------------------------
fetchLincsExp <- function()
{
  GG=list()
  for(i in 1:length(Lincs))
  {
    #print(i)
    
    list=data.frame(Lincs[,i])
    list$Name=rownames(Lincs)
    names(list)=c("DrugExp","GeneName")
    Data=list[which(abs(list$DrugExp)>=Cutoff),]
    
    # print(Data$GeneName)
    GG[[i]]<-Data
  }
  return(GG)
}


fetchDrugExp <- function(ind, Candidates)
{
    list=data.frame(Lincs[,ind])
    list$Gene=rownames(Lincs)
    names(list)=c("DrugExp","GeneName")
    Data=list[which(list$GeneName %in% Candidates),]
    
    return(Data)
}

fetchCandExp <- function(Candidates)
{
  Data=cancerExp[which(cancerExp$Gene %in% Candidates),]
  names(Data)=c("GeneName","CancerExp")
  return(Data)
}


Visualize<- function()
{
  library(visNetwork)
  
  V(G)$color="yellow"
  V(G)[which(V(G)$name %in% DrugGenes)]$color="blue"
  V(G)[which(V(G)$name %in% CancerGenes)]$color="red"
  V(G)[which(!V(G)$name %in% GeneNames$Name)]$size=8
  
  E(G)[which(E(G)$weight==-1)]$dashes=TRUE
  E(G)$color="black"
  #visIgraph(G,physics = TRUE, smooth = TRUE)
  
  net <- visIgraph(G,physics = TRUE,smooth = TRUE) %>% 
    visOptions(highlightNearest = TRUE) %>%
    visNodes(shadow = FALSE)
  
  net$width=1400
  net$height=800
  visSave(net, file = paste0(getwd(),"/Graphs/",netcount,".html") ,background = "white",selfcontained = TRUE)
  #visSave(net, file = paste0(getwd(),"/Graphs/test.html"), background = "white")
  
}

upStream<-function(Node)
{
  v=which(DDM[,which(colnames(DDM)==Node)]!=0)
  if(length(v)==0)
    return(NULL)
  else
    return(v)
}

DownStreamNodes <- function(Node)
{
  length(ego(G,order = 1000,mode = "out",nodes = Node)[[1]])-1
}


### TODO: write a recursive using ComputePerturbation & upStream
ComputeDrugPerturbation<-function(ind)
{
  #st=names(which(colSums(DDM)==0))
  #G[which(G$GeneName %in% st),]$Score=1
  if(DPerto[ind]!=0)
    Pert=DPerto[ind]
  else
  {
    cn=colnames(DDM)[ind]
    #print(cn)
    
    S=which(DDM[,ind]!=0)
    if(length(S)==0)
      Pert=AllExp$DrugExp[ind]
    else
    {
      temp=0
      for(gg in S)
      {
        cmp <- ComputeDrugPerturbation(gg)
        dst <- DownStreamNodes(colnames(DDM)[gg])
        
        temp <- temp+(DDM[gg,ind]*cmp/dst)
      }
      
      Pert <- AllExp$DrugExp[ind]+temp
    }
    DPerto[ind]=Pert
    
  }
  return(as.numeric(Pert))
}

ComputeCancerPerturbation<-function(ind)
{
  #st=names(which(colSums(DDM)==0))
  #G[which(G$GeneName %in% st),]$Score=1
  
  if(CPerto[ind]!=0)
    Pert=CPerto[ind]
  else
  {
    cn=colnames(DDM)[ind]
    #print(cn)
    
    S=which(DDM[,ind]!=0)
    if(length(S)==0)
      Pert=AllExp$CancerExp[ind]
    else
    {
      temp=0
      for(gg in S)
      {
        cmp <- ComputeCancerPerturbation(gg)
        dst <- DownStreamNodes(colnames(DDM)[gg])
        
        temp <- temp+(DDM[gg,ind]*cmp/dst)
      }
      Pert <- AllExp$CancerExp[ind]+temp
    } 
    CPerto[ind]=Pert
    
  }
  return(Pert)
}








#========================


#--------- Load Lincs Data ---------
load('Lincs.RData')
AllDrugGenes <- rownames(Lincs) <- Genes$V1
rm(Genes)
LincsGenes <- fetchLincsExp()
#-----------------------------------


#----------------- load BC scores for a given subtype --------------
load("BC_Scores.RData")
cancerExp=BC_Scores[,Subtype]
AllCancerGenes=BC_Scores$GeneName
cancerExp=cbind(data.frame(BC_Scores$GeneName),data.frame(cancerExp))
names(cancerExp)=c("Gene","Exp")
#-------------------------------------------------------------------


#--------- load Candidate genes of a given subtype ---------
load("UnifiedCandidateGenes.RData")
Candidates=AllCandidateGenes[which(AllCandidateGenes$Subtype==Subtype),]$Gene
CancerGenes=as.character(Candidates)
#-----------------------------------------------------------


#--------------- load Kegg Data ---------------------
load("KEGG_Matrix.RData")


#xx=which(GeneNames$Name %in% AllCancerGenes & GeneNames$Name %in% AllDrugGenes)

#M=M[xx,xx]
#M_original=M_original[xx,xx]
#GeneNames=data.frame(GeneNames[xx,])
#names(GeneNames)="Name"
#----------------------------------------------------

#Chem=unique(DrugGene$ChemicalID)
#NotFound=0
#Found=0



Pert=data.frame()
library(tictoc)
cat("\014")

for(netcount in Cutted:length(LincsGenes))
{
  print(netcount)
  
  DrugGenes=LincsGenes[[netcount]]$GeneName

  UnifiedGenes=union(CancerGenes,DrugGenes)
  DrugExpr <- fetchDrugExp(netcount,UnifiedGenes)
  CanExpr <- fetchCandExp(UnifiedGenes)
  AllExp=merge(CanExpr,DrugExpr,by='GeneName')
  
  x=merge(data.frame(DrugGenes),data.frame(AllExp$GeneName),by.x='DrugGenes',by.y='AllExp.GeneName')
  DrugGenes=x$DrugGenes
  
  x=merge(data.frame(Candidates),data.frame(AllExp$GeneName),by.x='Candidates',by.y='AllExp.GeneName')
  CancerGenes=x$Candidates
  
  
  #cval=Chem[netcount]
  
  #cval="C521487"  #ATTENTION ----------------------------
  
  #genes=DrugGene[which(DrugGene$ChemicalID==cval),]$GeneSymbol

  
  DDM=matrix(0,nrow(AllExp),nrow(AllExp))
  rownames(DDM) <- colnames(DDM) <- AllExp$GeneName
  
  if(VisualizeGraph==TRUE)
  {
    cat("\014")
    #print(genes)
    
    cat(netcount)
    cat(" ----- ")
    cat(length(DrugGenes))
    cat(" <> ")
    cat(length(CancerGenes))
    cat("\n")
  }
  
  #ll=DrugGenes[which(DrugGenes %in% GeneNames$Name)]
  
  tic('create DNN...')
  
  for(dg in DrugGenes)
  {
    start=which(GeneNames$Name==dg)
    if(length(start)>0)
    {
      AllMGenes=data.frame()
      DDM0=matrix(0,1,1)
      rownames(DDM0) <- colnames(DDM0)<- dg
      
      for(g in CancerGenes)
      {

        
        end=which(GeneNames$Name==g)
        
        if(length(end)>0)
        {
          
          #if(!is.na(ShortestPath$length[start,end]))
          if(dg %in% AllDrugGenes & g %in% AllCancerGenes )
          {
            
            
            if(start!=end)
            {
              if(BidirectionalPath==TRUE & !is.na(ShortestPath$length[end,start]) & is.na(ShortestPath$length[start,end]))
                P=extractPath(ShortestPath, end, start)
              else
                P=extractPath(ShortestPath, start, end)
              NewM=matrix(data = 0,nrow = length(P),ncol = length(P))
              
              MGenes=data.frame()
              
              Counter=1
              
              for(gg in 1:(length(P)-1))
              {
                if(row.names(M_original)[gg]%in% AllCancerGenes & row.names(M_original)[gg]%in% AllDrugGenes )
                {
                  Counter=Counter+1
                  
                  CG=which( AllCancerGenes== row.names(M_original)[gg] )
                  DG=which( AllDrugGenes== row.names(M_original)[gg] )
                  
                  CMexp=cancerExp[CG,]$Exp
                  DMexp=Lincs[DG,netcount]
                  MG=row.names(M_original)[gg]
                  
                  xx=cbind(MG,CMexp,DMexp)
                  MGenes=rbind(MGenes,xx)
                  
                  NewM[gg,gg+1]=M_original[P[gg],P[gg+1]]
                }
              }
              if (Counter==length(P))
              {
                rownames(NewM) <- colnames(NewM) <- row.names(M_original[P,P])
                
                DDM0=NetUnion2(DDM0,NewM)
                
                AllMGenes=rbind(AllMGenes,MGenes)
              }
            }
            
            else
            {
              NewM=matrix(data=M_original[start,end],1,1)
              rownames(NewM) <- colnames(NewM)<- g
              
              DDM0=NetUnion2(DDM0,NewM)
            }
            
            
          }  
        }
        
      }
      DDM=NetUnion2(DDM,DDM0)
      
      names(AllMGenes)=c('GeneName','CancerExp','DrugExp')
      AllExp=rbind(AllExp,AllMGenes)
    }
  }
  
  
  toc()
  
  for (res in 1:nrow(DDM))
    DDM[res,res]=0  
  

  G=graph_from_adjacency_matrix(DDM,mode = "directed",weighted = TRUE)
  
  if(VisualizeGraph==TRUE)
    Visualize()
  
  
  
  AllExp=unique(AllExp[order(AllExp$GeneName),])
  
  
  if(netcount==56)
    browser()
  
  DrugPert=numeric(nrow(DDM))  
  CancerPert=numeric(nrow(DDM))
  
  
  DPerto=numeric(nrow(DDM))
  CPerto=numeric(nrow(DDM))
  

  for(g in 1:nrow(DDM))
  {
    #print(g)
    DrugPert[g] <- as.numeric(ComputeDrugPerturbation(g))
    CancerPert[g] <- as.numeric(ComputeCancerPerturbation(g))
    #DrugPert[g] <- CancerPert[g] <- 0
  }
  toc()
  tic('statistical testing...')
  KS <- ks.test(DrugPert,CancerPert)$p.value
  Pearson <- cor(DrugPert,CancerPert)
  Dname=paste0(Drugs$pert_iname[netcount],'_',Drugs$pert_idose[netcount],'_',Drugs$pert_itime[netcount])
  
  td=data.frame(Dname,KS,Pearson)
  names(td)=c('Drug','KS','Pearson')
  
  Pert=rbind(Pert,td)
  
  
  #for(g in Candidates$GeneName)
  #{
  #   target=which(colnames(DDM)==g)
  #   if(length(target)>0)
  #   Pert=ComputePerturbation(target)
  # }
  # G=data.frame(row.names(DDM))
  # names(G)="GeneName"
  # G$Score=0
  
  
}

