pathway_show <- function(human_gem,flux,name){
  pathway<-unique(unlist(human_gem$SUBSYSTEM));
  pathway_score<-list();
  for (i in pathway){
    path=i;
    activity_score<-c();
    for (d in 1:ncol(flux)){
      activity_score[d]<-mean(abs(flux[which(unlist(human_gem$SUBSYSTEM)==i),d]));
    } 
    pathway_score[[i]]<-activity_score;
  }
  all_pathway_score<-as.data.frame(do.call(rbind,pathway_score));
  colnames(all_pathway_score)=colnames(flux);
  mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256);
  pheatmap::pheatmap(all_pathway_score,cluster_cols = T,border_color = NA,
                     treeheight_row = 0,treeheight_col = 4,
                     fontsize_col  = 5,fontsize_row = 5,cellwidth = 18,cellheight = 5,
                     color = rev(mapal),scale = "row",fontsize = 3,main = name,
                    )
  #return(all_pathway_score)
} 
  