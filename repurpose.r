library(e1071)
library(igraph)
library(visNetwork)
library(igraph)

#setwd("D:/BreastCancerResearch/DrugRepositioning/R/DNN/")

source("NetUnion2.r")


load("DrugGene.RData")
load("KEGG_Matrix.RData")
load("CandidateGenes.RData")

Candidates=Subtype1


#---------------------------------------

upStream<-function(DDM,Node)
{
  v=which(DDM[,which(colnames(DDM)==Node)]!=0)
  if(length(v)==0)
    return(NULL)
  else
    return(v)
}

DownStreamNodes <- function(Network , Node)
{
  length(ego(Network,order = 1000,mode = "out",nodes = Node)[[1]])-1
}


### TODO: write a recursive using ComputePerturbation & upStream
ComputePerturbation<-function(ind)
{
  #st=names(which(colSums(DDM)==0))
  #G[which(G$GeneName %in% st),]$Score=1
  
  S=which(DDM[,ind]!=0)
  if(length(S)==0)
    Pert[ind]=ExpressionChange(g)
  else
  {
    temp=0
    for(gg in S)
    {
      temp=temp+(DDM[gg,ind]*ComputePerturbation(gg)/DownStreamNodes(gg))
    }
    Pert[ind]=ExpressionChange(g)+temp
  } 
  return(Pert)
}


#========================


GeneNames=data.frame(read.csv("Genes.csv"))
names(GeneNames)="Name"


row.names(M)<-colnames(M)<-GeneNames$Name
row.names(M_original)<-colnames(M_original)<-GeneNames$Name

Chem=unique(DrugGene$ChemicalName)
#NotFound=0
#Found=0

for(netcount in 1:length(Chem))
{
  
  

  cval=Chem[netcount]
  
  Counter=0

  genes=DrugGene[which(DrugGene$ChemicalName==cval),]$GeneSymbol

  Dmatrix=matrix(0,length(genes),length(genes))
  rownames(Dmatrix) <- colnames(Dmatrix) <- genes
  
  Cmatrix=matrix(0,length(Candidates$GeneName),length(Candidates$GeneName))
  rownames(Cmatrix) <- colnames(Cmatrix) <- Candidates$GeneName
  
  DDM=NetUnion2(Dmatrix,Cmatrix)
  
  cat("\014")
  cat(netcount)
  cat(" ----- ")
  cat(length(genes))
  cat(" <> ")
  cat(length(Candidates$GeneName))
  cat("\n")
  
  for(dg in genes)
  {
    
    start=which(GeneNames$Name==dg)
    for(g in Candidates$GeneName)
    {

      end=which(GeneNames$Name==g)
        
        if(length(start)>0 & length(end)>0)
        {
          #if(!is.na(ShortestPath$length[start,end]))
          {
          
            Counter=Counter+1
            
            if(start!=end)
            {
              P=extractPath(ShortestPath, start, end)
              NewM=matrix(data = 0,nrow = length(P),ncol = length(P))
              for(gg in 1:(length(P)-1))
              {
                NewM[gg,gg+1]=M_original[P[gg],P[gg+1]]
              }
              rownames(NewM) <- colnames(NewM) <- row.names(M_original[P,P])
              
              DDM=NetUnion2(DDM,NewM)
            }
            
            else
            {
              NewM=matrix(data=M_original[start,end],1,1)
              rownames(NewM) <- colnames(NewM)<- g
              DDM=NetUnion2(DDM,NewM)
              
            }
            
            
          }  
        }
        
    }
  }
  
  for (res in 1:nrow(DDM))
    DDM[res,res]=0  
  
  G=graph_from_adjacency_matrix(DDM,mode = "directed",weighted = TRUE)
  V(G)$color="yellow"
  V(G)[which(V(G)$name %in% genes)]$color="blue"
  V(G)[which(V(G)$name %in% Candidates$GeneName)]$color="red"
  
  E(G)[which(E(G)$weight==-1)]$dashes=TRUE
  E(G)$color="black"
  #visIgraph(G,physics = TRUE, smooth = TRUE)
  net <- visIgraph(G,physics = TRUE,smooth = TRUE) %>% visOptions(highlightNearest = FALSE) 
  visSave(net, file = paste0(getwd(),"/Graphs/",cval,".html") ,background = "white")
  #visSave(net, file = paste0(getwd(),"/Graphs/test.html"), background = "white")
  
  
  #print(c)
  rr=0

}  
  # for(g in Candidates$GeneName)
  # {
  #   target=which(colnames(DDM)==g)
  #   if(length(target)>0)
  #     
  #     Pert=ComputePerturbation(target)
  # 
  #   }
  # 
  # G=data.frame(row.names(DDM))
  # names(G)="GeneName"
  # G$Score=0
  # 
  
#}

