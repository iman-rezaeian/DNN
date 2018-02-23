
load('FinalResults_AllSubtypes.RData')

Subtype1=Subtype1[order(Subtype1$Pearson),]
Subtype2=Subtype2[order(Subtype2$Pearson),]
Subtype3=Subtype3[order(Subtype3$Pearson),]
Subtype4=Subtype4[order(Subtype4$Pearson),]
Subtype5=Subtype5[order(Subtype5$Pearson),]
Subtype6=Subtype6[order(Subtype6$Pearson),]
Subtype7=Subtype7[order(Subtype7$Pearson),]
Subtype8=Subtype8[order(Subtype8$Pearson),]
Subtype9=Subtype9[order(Subtype9$Pearson),]
Subtype10=Subtype10[order(Subtype10$Pearson),]

AllDrugs=cbind(Subtype1[,c('Drug','Pearson')],
               Subtype2[,c('Drug','Pearson')],
               Subtype3[,c('Drug','Pearson')],
               Subtype4[,c('Drug','Pearson')],
               Subtype5[,c('Drug','Pearson')],
               Subtype6[,c('Drug','Pearson')],
               Subtype7[,c('Drug','Pearson')],
               Subtype8[,c('Drug','Pearson')],
               Subtype9[,c('Drug','Pearson')],
               Subtype10[,c('Drug','Pearson')])

N1=paste0('Subtype ',1:10)
N2=paste0('Corr ',1:10)
colnames(AllDrugs)=rbind(N1,N2)
row.names(AllDrugs)=1:nrow(AllDrugs)


AllCorr=AllDrugs[,seq(2,20,2)]
AllD=AllDrugs[,seq(1,19,2)]


DD=unique(AllD)
for (i in 1:nrow(DD)) {
  which(Subtype2$Drug==Subtype1[2,]$Drug)
}


write.csv(AllDrugs,'Top20Drugs.csv')

row.names(Subtype1)=1:nrow(Subtype1)
row.names(Subtype2)=1:nrow(Subtype2)
row.names(Subtype3)=1:nrow(Subtype3)
row.names(Subtype4)=1:nrow(Subtype4)
row.names(Subtype5)=1:nrow(Subtype5)
row.names(Subtype6)=1:nrow(Subtype6)
row.names(Subtype7)=1:nrow(Subtype7)
row.names(Subtype8)=1:nrow(Subtype8)
row.names(Subtype9)=1:nrow(Subtype9)
row.names(Subtype10)=1:nrow(Subtype10)


which(Subtype2$Drug==Subtype1[2,]$Drug)




