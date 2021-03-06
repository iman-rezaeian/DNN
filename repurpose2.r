rm(list=ls(all=TRUE))

cat("\014")



#==========================
Cutted=1   # 5467 12017
Cutoff=9  # cutoff for drug expression genes (between 0 and 10) - the higher the value, the algorithm keeps less genes
Cutoff2=50

BidirectionalPath=FALSE
VisualizeGraph=FALSE
DEBUG=FALSE
#==========================
library(tictoc)
library(e1071)
library(igraph)
library(visNetwork)
library(igraph)
library(data.table)

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
    
    xx=list[order(abs(list$DrugExp),decreasing = TRUE),]
    
    Data2 = xx[1:Cutoff2,]
    Data=list[which(abs(list$DrugExp)>=Cutoff),]
    
    # print(Data$GeneName)
    GG[[i]]<-Data2
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
  V(G)[which(!V(G)$name %in% GeneNames)]$size=8
  
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

  if(DPerto[ind]!=0)
  {
    Pert=DPerto[ind]
  }
  else
  {
    DinProcess[ind]<<-1
    
    cn=colnames(DDM)[ind]
    #print(cn)
    
    S=which(DDM[,ind]!=0)
    if(length(S)==0)
      Pert=as.numeric(AllExp$DrugExp[ind])
    else
    {
      temp=0
      for(gg in S)
      {
        if(DinProcess[gg]==0)
        {
          cmp <- ComputeDrugPerturbation(gg)
          dst <- DownStreamNodes(colnames(DDM)[gg])
          
          temp <- temp+(DDM[gg,ind]*cmp/dst)
          
        }
      }

      Pert = as.numeric(AllExp$DrugExp[ind])+temp
    }
    DPerto[ind]<<-Pert
    DinProcess[ind]<<-0
    
  }
  
  return(Pert)
}

ComputeCancerPerturbation<-function(ind)
{
  
  if(CPerto[ind]!=0)
  {
    Pert=CPerto[ind]
  }
  else
  {
    CinProcess[ind]<<-1
    
    cn=colnames(DDM)[ind]
    #print(cn)
    
    S=which(DDM[,ind]!=0)
    if(length(S)==0)
      Pert=as.numeric(AllExp$CancerExp[ind])
    else
    {
      temp=0
      for(gg in S)
      {
        if(CinProcess[gg]==0)
        {
          cmp <- ComputeCancerPerturbation(gg)
          dst <- DownStreamNodes(colnames(DDM)[gg])
          
          temp <- temp+(DDM[gg,ind]*cmp/dst)
          
        }
      }
      
      Pert = as.numeric(AllExp$CancerExp[ind])+temp
    }
    CPerto[ind]<<-Pert
    CinProcess[ind]<<-0
    
  }
  
  return(Pert)
  
}










#========================
for(Subtype in 1:10)
{
    print(Subtype)
  
    SubT0=as.character(floor((Subtype+1)/2))
    Stag0=as.character((Subtype-1)%%2)
    ff=paste0('FinalResults_Subtype',Subtype,'.RData')
    
    #--------- Load Lincs Data ---------
    load('Lincs.RData')
    AllDrugGenes <- rownames(Lincs) <- Genes$V1
    rm(Genes)
    Drugs$index=1:nrow(Drugs)
    Drugs$pert_idose = as.numeric(gsub("[a-z A-Z]", "", Drugs$pert_idose))
    Drugs = data.table(Drugs)
    fmin = Drugs[ , .SD[which.min(pert_idose)], by = list(pert_iname,pert_itime)]
    fmax = Drugs[ , .SD[which.max(pert_idose)], by = list(pert_iname,pert_itime)]
    f = unique(rbind(fmin,fmax))
    I = sort(f$index)
    Drugs = Drugs[I,]
    Lincs = Lincs[,I]
    LincsGenes <- fetchLincsExp()
    #-----------------------------------
    
    
    #----------------- load BC scores for a given subtype --------------
    load("BC_Scores_10Subtypes.RData")
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
    
    
    #xx=which(GeneNames %in% AllCancerGenes & GeneNames %in% AllDrugGenes)
    
    #M=M[xx,xx]
    #M_original=M_original[xx,xx]
    #GeneNames=data.frame(GeneNames[xx,])
    #names(GeneNames)="Name"
    #----------------------------------------------------
    
    #Chem=unique(DrugGene$ChemicalID)
    #NotFound=0
    #Found=0
    
    
    
    AllPert=data.frame()
    
    for(netcount in Cutted:length(LincsGenes))
    {
      #print(netcount)
      
    
      
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
      
      #ll=DrugGenes[which(DrugGenes %in% GeneNames)]
      
      #tic('create DNN...')
      
      for(dg in DrugGenes)
      {
    
        start=which(GeneNames==dg)
        
    
        
        if(length(start)>0)
        {
          #if(start==2826)
           # browser()
          
          AllMGenes=data.frame()
          
          DDM0=matrix(0,1,1)
          rownames(DDM0) <- colnames(DDM0)<- dg
          
          
          for(g in CancerGenes)
          {
    
            end=which(GeneNames==g)
            
            if(length(end)>0)
            {
              
                if(start!=end)
                {
    
                  
                  if(BidirectionalPath==TRUE & !is.na(ShortestPath$length[end,start]) & is.na(ShortestPath$length[start,end]))
                    P=extractPath(ShortestPath, end, start)
                  else
                    P=extractPath(ShortestPath, start, end)
    
                  MGenes=data.frame()
                  
                  Counter=1
                  
                  for(gg in 1:length(P))
                  {
                      Counter=Counter+1
                      
                      CG=which( AllCancerGenes== row.names(M_original)[P[gg]])
                      DG=which( AllDrugGenes== row.names(M_original)[P[gg]])
                      
                      CMexp=cancerExp[CG,]$Exp
                      DMexp=Lincs[DG,netcount]
                      MG=row.names(M_original)[P[gg]]
                      
                      xx=cbind(MG,CMexp,DMexp)
                      MGenes=rbind(MGenes,xx)
                      
                  }
                  NewM=M_original[P,P]
                  rownames(NewM) <- colnames(NewM) <- row.names(M_original[P,P])
                  
                  DDM0=NetUnion2(DDM0,NewM)
                  
                  AllMGenes=rbind(AllMGenes,MGenes)
                }
                
                else
                {
                  NewM=matrix(data=M_original[start,end],1,1)
                  rownames(NewM) <- colnames(NewM)<- g
                  
                  DDM0=NetUnion2(DDM0,NewM)
                }
                
                
               
            }
            
          }
          DDM=NetUnion2(DDM,DDM0)
          
          names(AllMGenes)=c('GeneName','CancerExp','DrugExp')
          AllExp=unique(rbind(AllExp,AllMGenes))
          
          
        }
      }
      
      
      #toc()
      
      for (res in 1:nrow(DDM))
        DDM[res,res]=0  
      
    
      G=graph_from_adjacency_matrix(DDM,mode = "directed",weighted = TRUE)
      
      if(VisualizeGraph==TRUE)
        Visualize()
      
      
      
      AllExp=unique(AllExp[order(AllExp$GeneName),])
      
    
      
      DrugPert=numeric(nrow(DDM))  
      CancerPert=numeric(nrow(DDM))
      
      
      DPerto=numeric(nrow(DDM))
      CPerto=numeric(nrow(DDM))
      
      DinProcess=numeric(nrow(DDM))
      CinProcess=numeric(nrow(DDM))
      
      if(netcount==56 & DEBUG==TRUE)
        browser()
      
    
      for(g in 1:nrow(DDM))
      {
        #print(g)
        DrugPert[g] <- as.numeric(ComputeDrugPerturbation(g))
        CancerPert[g] <- as.numeric(ComputeCancerPerturbation(g))
        #DrugPert[g] <- CancerPert[g] <- 0
      }
      #toc()
    
      KS <- ks.test(DrugPert,CancerPert)
      KS_pval = KS$p.value
      KS_D = KS$statistic
      Pearson <- cor(DrugPert,CancerPert,method = 'pearson')
      Kendal <- cor(DrugPert,CancerPert,method = 'kendall')
      Spearman <- cor(DrugPert,CancerPert,method = 'spearman')
    
      Dname=paste0(Drugs$pert_iname[netcount],'_',Drugs$pert_idose[netcount],'_',Drugs$pert_itime[netcount])
      
      td=data.frame(Dname,KS_pval,KS_D,Pearson,Kendal,Spearman)
      names(td)=c('Drug','KS(pval)','KS(D)','Pearson','Kendal','Spearman')
      
      AllPert=rbind(AllPert,td)
      
      
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
    
    save(AllPert,file = ff)

}