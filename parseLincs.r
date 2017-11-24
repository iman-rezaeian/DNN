library(data.table)

setwd("D:/BreastCancerResearch/DrugRepositioning/LINC_interpreter")


parseLincs <- function()
{
  L=fread("Lincs.csv",sep = ",")
  
  
  
}

parseLincsMetaData <- function()
{
  Drugs=fread("Lincs_drugs_All.txt",sep = "\t")
  return(Drugs)
}

#---------------------------

Drugs=parseLincsMetaData()

Drugs2=Drugs[which(Drugs$pert_idose!=-666),]


cellInfo=fread("LINCS_cell_info.txt",sep = "\t")
sigInfo=fread("LINCS_sig_info.txt",sep = "\t")

br=cellInfo[which(cellInfo$primary_site=="breast"),]$cell_id

xx=which((sigInfo$cell_id %in% br) & (sigInfo$pert_idose!=-666))
xx2=(sigInfo$cell_id %in% br) & (sigInfo$pert_idose!=-666)

xx2[xx2==TRUE]=NA
xx2[xx2==FALSE]="NULL"


Lincs=fread("Lincs.csv",sep=",", colClasses=xx2)  # nrows = 10, 
names(Lincs)=as.character(xx)
Drugs=Drugs[xx,]

save(Lincs,Drugs,file = "Lincs.RData")
#lincs2=Lincs[,..xx]



