
library("illuminaHumanv4.db") 

  setwd("D:/BreastCancerResearch/DrugRepositioning/data/BC expression")
  
  Label=read.table('classlabels.txt',header = TRUE,sep = "\t")
  Label=Label$IntClustMemb
  
  NExp=read.table('normals_ExpressionMatrix.txt',header = TRUE,sep = "\t")
  
  probes=rownames(NExp)
  
  S <- data.frame(unlist(mapIds(illuminaHumanv4.db, probes, "SYMBOL","PROBEID")))
  
  names(S)="GeneName"
  S$ProbeID=probes
  S2=data.frame(S[complete.cases(S), ])
  row.names(S2)=NULL
  
  RowNames=data.frame(rownames(NExp))
  names(RowNames)='ID'
  RowNames$Index=as.numeric(rownames(RowNames))
  
  B=merge(RowNames,S2,by.x='ID',by.y='ProbeID',all=TRUE)
  B2 = B[order(B$Index),]
  
  Genes=unique(B$GeneName)
  Genes=Genes[which(!is.na(Genes))]

  
  finalScores=matrix(0,nrow=length(Genes),ncol=10)
  
  CExp_o=read.table('discovery_ExpressionMatrix.txt',header = TRUE,sep = "\t")
  
  for (Subtype in 1:10) 
  {
    print(Subtype)
    
    SType= which(Label==Subtype)
    
    CExp=CExp_o[,SType]
    
    Gval=array(0,dim=length(Genes))
    
    CancerM=array(0,dim=length(Genes))
    NormalM=array(0,dim=length(Genes))
    
    for (i in 1:length(Genes)) {
      L=as.numeric(B[which(B$GeneName==Genes[i]),]$Index)
      # L now contains the row numbers of matrix containing a gene
      Cancer=unlist(CExp[L,])
      Normal=unlist(NExp[L,])
      #scale(unlist(NExp[L,]),center = TRUE, scale = TRUE)
      CancerM[i]=median(Cancer)
      NormalM[i]=median(Normal)
      
      Gval[i]=median(Cancer)/median(Normal)
      
      
      #How to model Level5 on our data here?
  
    }
    CancerM=CancerM[which(!is.na(CancerM))]
    NormalM=NormalM[which(!is.na(NormalM))]
    
    cx=scale(CancerM,center = TRUE, scale = TRUE)
    nx=scale(NormalM,center = TRUE, scale = TRUE)
    f=cx-nx
    sc=10/max(abs(f))
    finalScores[,Subtype]=f*sc
    
  }
  
  BC_Scores=data.frame(finalScores)
  BC_Scores$GeneName=Genes
  
  save(BC_Scores,file="BC_Scores.RData")








