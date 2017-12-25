load('Lincs.RData')
Genes=data.frame(Genes)
names(Genes)="Gene"


load('AllCandidateGenes.RData')

GG=list()
for(i in 1:length(Lincs))
{
  print(i)
  
  list=data.frame(sort(Lincs[,i],index.return=TRUE,decreasing = TRUE))
  list[,3]=Genes$Gene
  names(list)=c("Exp","Ind","GeneName")
  Data=list[which(abs(list$Exp)>2),]
  
  print(Data$GeneName)
  GG[[i]]<-Data
}
LincsGenes=GG
save(LincsGenes,file = "LincsGenes.RData")


