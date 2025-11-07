setwd("C:/Users/86189/Desktop/neo_project")
gc()

View(A549_flux_total_total)
stab_A549=A549_flux_total_total
colnames(stab_A549)=c("study1","study2","study3","study4","study5","study6","study7","study8")

stab_A549$Reaction_ID=rownames(stab_A549)
stab_A549=merge(stab_A549,experi_A549,by="Reaction_ID")
rownames(stab_A549)=stab_A549$Reaction_ID
stab_A549=subset(stab_A549,select=-c(Reaction_ID))

correlation=rep(0,length(colnames(stab_A549)))
for (i in 1:length(colnames(stab_A549))){
  correlation[i]=cor(as.double(stab_A549[,i]),as.double(stab_A549$experimental_flux),method="spearman")
}
stab_A549_correlation=rbind(stab_A549,correlation)


stab_A549=A549_flux_total_total
colnames(stab_A549)=c("study1","study2","study3","study4","study5","study6","study7","study8")

stab_A549$Reaction_ID=rownames(stab_A549)
stab_A549=merge(stab_A549,experi_A549,by="Reaction_ID")
rownames(stab_A549)=stab_A549$Reaction_ID
stab_A549=subset(stab_A549,select=-c(Reaction_ID))

correlation=rep(0,length(colnames(stab_A549)))
for (i in 1:length(colnames(stab_A549))){
  correlation[i]=cor(as.double(stab_A549[,i]),as.double(stab_A549$experimental_flux),method="spearman")
}
stab_A549_correlation=rbind(stab_A549,correlation)

######

stab_K562=K562_flux_total
colnames(stab_K562)=c("study1","study2","study3","study4","study5")

stab_K562$Reaction_ID=rownames(stab_K562)
stab_K562=merge(stab_K562,experi_K562,by="Reaction_ID")
rownames(stab_K562)=stab_K562$Reaction_ID
stab_K562=subset(stab_K562,select=-c(Reaction_ID))

correlation=rep(0,length(colnames(stab_K562)))
for (i in 1:length(colnames(stab_K562))){
  correlation[i]=cor(as.double(stab_K562[,i]),as.double(stab_K562$experimental_flux),method="spearman")
}
stab_K562_correlation=rbind(stab_K562,correlation)


## Bar plot
stab_weighted_total <- cbind(
  t(stab_A549_correlation[nrow(stab_A549_correlation), 1:(ncol(stab_A549_correlation) - 1)]),
  t(stab_K562_correlation[nrow(stab_K562_correlation), 1:(ncol(stab_K562_correlation) - 1)])
)


colnames(stab_weighted_total)=c("A549","K562");
stab_weighted_total=as.data.frame(stab_weighted_total)
stab_weighted_total$weight=as.numeric(rownames(stab_weighted_total))
stab_weighted_total_long=pivot_longer(as.data.frame(stab_weighted_total),cols=1:2,names_to="Cell_line",values_to="Correlation")

ggplot(data=stab_weighted_total_long,aes(x=weight,y=Correlation,fill=Cell_line)) + geom_bar(stat="identity",position="dodge") + theme_minimal()+labs(title="Weighted mean")






# ################## check the pathway activity variation
# colnames(A549_flux_total_total_total)=c("sample1","sample2","sample3","sample4","sample5","sample6","sample7","sample8","sample9","sample10","sample11","sample12","sample13","sample14","sample16","sample17","sample18","sample19")
# colnames(K562_flux_total_total)=c("sample1","sample2","sample3","sample4","sample5","sample6","sample7","sample8","sample9","sample10","sample11","sample12","sample13","sample14","sample16","sample17","sample18","sample19","sample20","sample21","sample22")
# score_A549=pathway_comp(human_gem,A549_flux_total_total_total,"A549")
# score_K562=pathway_comp(human_gem,K562_flux_total_total,"K562")
# ## compute the pathway activity variation for each row
# score_A549$mean=apply(score_A549,1,mean);score_A549$sd=apply(score_A549,1,sd)
# score_A549$cv=score_A549$sd/score_A549$mean
# 
# score_K562$mean=apply(score_K562,1,mean);score_K562$sd=apply(score_K562,1,sd)
# score_K562$cv=score_K562$sd/score_K562$mean
# 
# score_A549_cv=data.frame(pathway=rownames(score_A549),cv=score_A549$cv)
# score_K562_cv=data.frame(pathway=rownames(score_K562),cv=score_K562$cv)
# 
# 
# draw_cv_sorted(score_A549_cv,"A549",score_K562_cv,"K562",4)
# 
# 
# 
# ##############
# 
# rownames(flux_A549_all)=flux_A549_all$reaction
# flux_A549_all=subset(flux_A549_all,select=-c(reaction))
# rownames(flux_K562_all)=flux_K562_all$reaction
# flux_K562_all=subset(flux_K562_all,select=-c(reaction))
# 
# 
# score_A549_all=pathway_comp(human_gem,flux_A549_all,"A549_all")
# score_K562=pathway_comp(human_gem,flux_K562_all,"K562")
# ## compute the pathway activity variation for each row
# score_A549_all$mean=apply(score_A549_all,1,mean);score_A549_all$sd=apply(score_A549_all,1,sd)
# score_A549_all$cv=score_A549_all$sd/score_A549_all$mean
# 
# score_K562$mean=apply(score_K562,1,mean);score_K562$sd=apply(score_K562,1,sd)
# score_K562$cv=score_K562$sd/score_K562$mean
# 
# score_A549_all_cv=data.frame(pathway=rownames(score_A549_all),cv=score_A549_all$cv)
# score_K562_cv=data.frame(pathway=rownames(score_K562),cv=score_K562$cv)
# 
# 
# draw_cv_sorted(score_A549_all_cv,"A549_all",score_K562_cv,"K562",4)





