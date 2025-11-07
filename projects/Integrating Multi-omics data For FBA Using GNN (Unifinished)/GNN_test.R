
test_GNN <- function(gnn1){
gnn1=as.data.frame(gnn1)
if (ncol(gnn1)==2){
  # gnn1$X=rownames(gnn1);colnames(gnn1)=c("X","abundence")
  rownames(gnn1)=gnn1$X;gnn1=subset(gnn1, select=-c(X))
}
colnames(gnn1)="abundence"
# min-max normalization
gnn1$abundence=(gnn1$abundence-min(gnn1$abundence))/(max(gnn1$abundence)-min(gnn1$abundence))
library(METAFlux)
#gnn1=rna_A549_list$Sample3[,c(1,2)];rownames(gnn1)=gnn1[,1];gnn1=subset(gnn1, select=-c(SYMBOL))
scores<-calculate_reaction_score(gnn1)
temp=compute_flux(mras=scores,medium=human_blood)
temp=as.data.frame(temp)
temp$Reaction_ID=rownames(temp)
temp=merge(temp,experi_A549,by="Reaction_ID")
rownames(temp)=temp$Reaction_ID;temp=temp[,-1]
gc()
correlation=cor(as.double(temp$abundence),as.double(temp$experimental_flux),method="spearman")
print(correlation)
gc()
}

# import .csv
gnn1=as.data.frame(read.csv("stable_activities.csv", header=TRUE, sep=","))
test_GNN(gnn1)


gnn2=as.data.frame(read.csv("stable_activities_highD.csv", header=TRUE, sep=","))
test_GNN(gnn2)

gnn3=input_lung
colnames(gnn3)="abundence"
scores<-calculate_reaction_score(gnn3)
temp=compute_flux(mras=scores,medium=human_blood)
temp=as.data.frame(temp)
temp$Reaction_ID=rownames(temp)
temp=merge(temp,experi_A549,by="Reaction_ID")
rownames(temp)=temp$Reaction_ID;temp=temp[,-1]
gc()
correlation=cor(as.double(temp$abundence),as.double(temp$experimental_flux),method="spearman")
print(correlation)
gc()
