setwd("C://Users//86189//Desktop//neo_project//GEM")

genes_GEM=read.table("genes.tsv",header = T,sep = "\t",quote = "",stringsAsFactors = F)
genes_GEM=subset(genes_GEM,select = c("genes","geneSymbols"))
genes_GEM$genes=gsub("\"","",genes_GEM$genes);genes_GEM$geneSymbols=gsub("\"","",genes_GEM$geneSymbols)

GEM=read.table("Human-GEM.tsv",header = T,sep = "\t",quote = "",stringsAsFactors = F)

reactions_GEM=read.table("reactions.tsv",header = T,sep = "\t",quote = "",stringsAsFactors = F)
reactions_GEM=subset(reactions_GEM,select = c("rxns","rxnHMR2ID"))
reactions_GEM$rxns=gsub("\"","",reactions_GEM$rxns);reactions_GEM$rxnHMR2ID=gsub("\"","",reactions_GEM$rxnHMR2ID)

######################
experi_A549=as.data.frame(experi_A549)

## find common reactions in experi_A549$Reaction_ID
#reactions_GEM=reactions_GEM[reactions_GEM$rxnHMR2ID %in% experi_A549$Reaction_ID,]
#GEM=GEM[GEM$Rxn.name %in% reactions_GEM$rxns,]

experi_A549_GEM=merge(reactions_GEM,experi_A549,by.x = "rxnHMR2ID",by.y = "Reaction_ID")
experi_A549_GEM=subset(experi_A549_GEM,select = -c(rxnHMR2ID))

# save into .csv file
write.csv(experi_A549_GEM,"experi_A549_GEM.csv",row.names = F,col.names = T,quote = F)
