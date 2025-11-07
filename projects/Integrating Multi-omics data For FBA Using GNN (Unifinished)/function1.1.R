pathway_show_mutated <- function(human_gem,flux_leukemia_mutated,name){
  pathway<-unique(unlist(human_gem$SUBSYSTEM));
  pathway_score<-list();
  for (i in pathway){
    path=i;
    activity_score<-c();
    for (d in 1:ncol(flux_leukemia_mutated)){
      activity_score[d]<-mean(abs(flux_leukemia_mutated[which(unlist(human_gem$SUBSYSTEM)==i),d]));
    } 
    pathway_score[[i]]<-activity_score;
  }
  all_pathway_score<-as.data.frame(do.call(rbind,pathway_score));
  colnames(all_pathway_score)=colnames(flux_leukemia_mutated);
  diff_pathway_score <- data.frame(pathway = rownames(all_pathway_score))
  i=1
 while (i <= ncol(all_pathway_score)) {
    temp=(all_pathway_score[,i+1]-all_pathway_score[,i])/abs(all_pathway_score[,i]);
    diff_pathway_score[[colnames(all_pathway_score)[i]]] = temp;
    i=i+2;
 }
  rownames(diff_pathway_score) <- diff_pathway_score$pathway;diff_pathway_score <- diff_pathway_score[,-1]
  mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256);
  pheatmap::pheatmap(diff_pathway_score,cluster_cols = F,cluster_rows = T,border_color = NA,color = rev(mapal),scale = "column",fontsize = 3,main = name)
  #return(all_pathway_score)
} 
