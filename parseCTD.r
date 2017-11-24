parseCTD <- function()
{
  CTD = read.csv("CTD_chem_gene_ixns.csv",stringsAsFactors = FALSE)
  H=CTD[which(CTD$OrganismID==9606 | is.na(CTD$OrganismID)),]
  DrugGene= H[,c("ChemicalName","ChemicalID","GeneSymbol")]
  DrugGene=unique(DrugGene)
  
  #--------------------------------------
  save(DrugGene,file = "DrugGene.RData")
}
