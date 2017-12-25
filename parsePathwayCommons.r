library(data.table)
library(Matrix)

setwd("D:/BreastCancerResearch/DrugRepositioning/R/DNN")
#---------------------------------

# setWeigth(Pathway)
# {
#   if(Pathway$INTERACTION_TYPE==)
#   
# }

#---------------------------
load("PathwayCommons.RData")

#Pathway <- setWeigth(Pathway)

v=unique(union(Pathway$PARTICIPANT_A,Pathway$PARTICIPANT_B))
Data=Matrix(0,nrow=length(v),ncol = length(v),sparse = TRUE)

row.names(Data) <- colnames(Data) <- v
prog=0
for (i in 1:nrow(Pathway))
{
  if (prog<round(i/nrow(Pathway),2))
  {
    prog=round(i/nrow(Pathway),2)
    print(prog)
  }
  
  Data[which(row.names(Data)==Pathway$PARTICIPANT_A[i]),which(row.names(Data)==Pathway$PARTICIPANT_B[i])]=1
}



