## let's get into annotation files
setwd('D:/westlake_summer/database/');
getwd()

## Complex Portal
library(readr)
ComplexPortal=read_tsv("2025-07-24-09-25.tsv",col_names = TRUE)
colnames(ComplexPortal)=gsub(" ","_",colnames(ComplexPortal))
colnames(ComplexPortal)=gsub("#","",colnames(ComplexPortal))
length(unique(ComplexPortal[[1]]))
# filter
evi_codes=c("ECO:0000353","ECO:0005543","ECO:0005610","ECO:0005544","ECO:0005546")  # 2 tier highest confidence score
ComplexPortal$Evidence_Code=gsub("\\(.*\\)","",ComplexPortal$Evidence_Code)
ComplexPortal=ComplexPortal[ComplexPortal$Evidence_Code %in% evi_codes,]
ComplexPortal=as_tibble(ComplexPortal)
ComplexPortal$proteins=strsplit(ComplexPortal$Expanded_participant_list,"[|]")
ComplexPortal_slim=ComplexPortal[,c("Complex_ac","Recommended_name","proteins")]
for (i in 1:length(ComplexPortal_slim$proteins)){
  ComplexPortal_slim$proteins[[i]]=gsub("\\(.*\\)","",ComplexPortal_slim$proteins[[i]])
} 

#################### MAPPING GENES IN COMPLEX LEVEL 
## refer uniprot ID to gene symbol
library(openxlsx)
up2gene=read.xlsx("uniprotkb_taxonomy_id_559292_2025_07_28.xlsx",colNames = TRUE)
up2gene$gene=gsub(" .*","",up2gene$Gene.Names)
map_up2gene <- function(up_id,up2gene){
  # up_id: a vector of uniprot IDs
  # up2gene: a data frame with columns 'Entry' and 'gene'
  matches <- match(up_id, up2gene$Entry)
  genes <- up2gene$gene[matches]
  return(genes)
}
# apply the mapping function to the proteins column
ComplexPortal_slim$genes <- lapply(ComplexPortal_slim$proteins, map_up2gene, up2gene=up2gene)

# count #genes in each complex
ComplexPortal_slim$num_all_genes=lapply(ComplexPortal_slim$genes, function(x){length(unlist(x))})



# map genes to complex
genes2complex <- map_df(unlist(unique(yeast_roi_alpha$set_genotype)), function(x) {
  complex_list <- list()
  for (i in 1:length(ComplexPortal_slim$genes)) {
    if (x %in% ComplexPortal_slim$genes[[i]]) {
      complex_list[[length(complex_list) + 1]] <- ComplexPortal_slim$Complex_ac[i]
    }
  }
  tibble(gene = x, complex = list(complex_list))
})


genes2complex$num=sapply(genes2complex$complex, length)
# filter out genes with no complex
genes2complex <- genes2complex[genes2complex$num > 0, ]

# build dataframe for plotting
lifes_set_alpha_complex_df <- data.frame()
for (i in 1:length(genes2complex$gene)){
  g=genes2complex$gene[i]
  lifespan=lifes_set_alpha_df$mean_rela[lifes_set_alpha_df$gene==g]
  steep=lifes_set_alpha_df$steepness_rela[lifes_set_alpha_df$gene==g]
  complex_vector=unlist(genes2complex$complex[[i]])
  for (j in 1:length(complex_vector)){
    complex=complex_vector[j]
    # add new row to the dataframe
    lifes_set_alpha_complex_df <- rbind(lifes_set_alpha_complex_df, 
                                        data.frame(gene=g, lifespan_rela=lifespan, steepness_rela=steep, complex=complex))
  }
}
temp=counts_alpha[,c("gene","num_of_experi")]
# map num_of_experi to the new dataframe
lifes_set_alpha_complex_df$num_of_experi <- sapply(lifes_set_alpha_complex_df$gene, function(g) {
  if (g %in% temp$gene) {
    return(temp$num_of_experi[temp$gene == g])
  } else {
    return(NA)
  }
})
rm(temp)

# count the number of genes in each complex
library(tidyverse)
counts_complex <- lifes_set_alpha_complex_df %>%
  group_by(complex) %>%
  summarise(num_genes = n_distinct(gene), .groups = 'drop')


################!!! FILTER CRETERIA I !!!################
good_complex=unlist(counts_complex[counts_complex$num_genes > 1, "complex"])


## zoom into best complexes
selected_complexes=counts_complex$complex[order(-counts_complex$num_genes)][1:4]
lifes_set_alpha_complex_df$description <- sapply(lifes_set_alpha_complex_df$complex, function(x) {
  if (x %in% ComplexPortal_slim$Complex_ac) {
    return(ComplexPortal_slim$Recommended_name[ComplexPortal_slim$Complex_ac == x])
  } else {
    return(NA)
  }
})


### prepare for enrichment
complex_df=ComplexPortal_slim[ComplexPortal_slim$Complex_ac %in% good_complex,]
complex_df=merge(complex_df,counts_complex,by.x="Complex_ac",by.y="complex",all=FALSE)
complex_df$ratio_genes_exist=as.numeric(complex_df$num_genes)/as.numeric(complex_df$num_all_genes)
################!!! FILTER CRETERIA II !!!################
complex_df=complex_df[complex_df$ratio_genes_exist>0.5,]


############### SOFT enrichment
# NEW NEW: use new weights calculated from SR_weight.r
distance_df <- read.csv("./data/distance_df_new_weight.csv", header = TRUE, stringsAsFactors = FALSE)
# head(distance_df)
colnames(distance_df)
# check NA of distance_df$weight_new
anyNA(distance_df$weight_new)
print(distance_df$weight_new)
distance_weight_df=subset(distance_df,select=c(eta_rela, epsilon_rela, xc_rela, scaling_rela))
rownames(distance_weight_df)=distance_df$gene
# 
# 
# distance_weight_df$xc_epsilon_rela=pmin(distance_weight_df$xc_rela,distance_weight_df$epsilon_rela)
# 
# between_xc_epsilon <- lifes_set_alpha_bp_short_df
# between_xc_epsilon$xc_x=xc_x(between_xc_epsilon$mean_rela)
# between_xc_epsilon$epsilon_x=epsilon_x(between_xc_epsilon$mean_rela)
# # find rows with steepness_rela(y) between xc_x & epsilon_x
# between_xc_epsilon=subset(between_xc_epsilon, (xc_x<steepness_rela & steepness_rela<epsilon_x)| (xc_x>steepness_rela & steepness_rela > epsilon_x))
# between_xc_epsilon=as.character(between_xc_epsilon$gene);length(between_xc_epsilon)
# # set rows in between_xc_epsilon value 0
# distance_weight_df[between_xc_epsilon, "xc_epsilon_rela"] <- 0.001
# rm(between_xc_epsilon)
# 
# distance_df$xc_epsilon_rela=distance_weight_df$xc_epsilon_rela

# # see the propotion
# temp=colSums(distance_weight_df);temp=temp/sum(temp);label=paste(names(temp),round(100*temp,2),"%")
# pie(temp,labels = label,main="percentage of sum of weights of genes")

### calculate p-value for each complex
## IMPROVED: sample with probability proportional to weight_new
generate_null_dist <- function(n, N=1000, weight_df, weight_vector=NULL){
  # weight_df: matrix/df with genes as rows, parameters as columns
  # weight_vector: sampling probability for each gene (should sum to 1)
  #               if NULL, use uniform (original behavior)
  # n: number of genes to sample
  # N: number of simulations
  
  gene_list <- as.character(rownames(weight_df))
  result <- matrix(NA, nrow = N, ncol = ncol(weight_df))
  
  if (is.null(weight_vector)) {
    # Original: uniform sampling
    sampling_prob <- rep(1/nrow(weight_df), nrow(weight_df))
  } else {
    # New: weighted sampling by weight_new
    sampling_prob <- weight_vector / sum(weight_vector, na.rm=TRUE)
  }
  
  for (i in 1:N){
    # Sample n genes with probability proportional to weights
    # (without replacement within each simulation)
    gene_null <- sample(gene_list, size = n, replace = FALSE, prob = sampling_prob)
    temp <- weight_df[gene_null, ]
    scores <- colSums(temp)
    result[i, ] <- scores
  }
  
  colnames(result) <- colnames(weight_df)
  result <- as.data.frame(result)
  return(result)
}


for (i in 1:length(complex_df$Complex_ac)){
  n=complex_df$num_genes[i]
  genes=unlist(complex_df$genes[i])
  temp=na.omit(distance_weight_df[genes,])
  scores=colSums(temp)
  N_sample=10000
  
  # WEIGHTED SAMPLING: use weight_new as sampling probability
  all_genes <- rownames(distance_weight_df)
  
  # Match gene names in distance_df to get correct weights
  # (safer than using row indices, avoids NA from mismatched names)
  idx_in_df <- match(all_genes, distance_df$gene)
  sampling_weights <- distance_df$weight_new[idx_in_df]
  
  # Ensure no NA and all weights are positive
  sampling_weights <- pmax(sampling_weights, 1e-10)
  # sampling_weights <- sampling_weights / sum(sampling_weights)
  
  scores_null=generate_null_dist(n, N_sample, distance_weight_df, weight_vector=sampling_weights)
  count_vector <- lapply(1:length(scores), function(j) {
    sum(scores_null[[j]] < scores[j])
  })
  count_vector <- unlist(count_vector);names(count_vector)=names(temp)
  # complex_df$SOFT_eta_p[i]=(count_vector["eta_rela"]+1)/(N_sample+1)
  # complex_df$SOFT_epsilon_p[i]=(count_vector["epsilon_rela"]+1)/(N_sample+1)
  # complex_df$SOFT_xc_p[i]=(count_vector["xc_rela"]+1)/(N_sample+1)
  # complex_df$SOFT_scaling_p[i]=(count_vector["scaling_rela"]+1)/(N_sample+1)
  # complex_df$SOFT_xc_epsilon_p[i]=(count_vector["xc_epsilon_rela"]+1)/(N_sample+1)
  # NEW
  complex_df$SOFT_eta_p[i]=(count_vector["eta_rela"]+1)/(N_sample+1)
  complex_df$SOFT_epsilon_p[i]=(count_vector["epsilon_rela"]+1)/(N_sample+1)
  complex_df$SOFT_xc_p[i]=(count_vector["xc_rela"]+1)/(N_sample+1)
  complex_df$SOFT_scaling_p[i]=(count_vector["scaling_rela"]+1)/(N_sample+1)
  if (i %% 5 == 0) {
    cat(sprintf("✓ Processed %d/%d complexes\n", i, length(complex_df$Complex_ac)))
  }
}

# remove proteins and genes columns
# complex_df=subset(complex_df,select=-c(proteins,genes))

# save into .xlsx
# setwd("./data/")
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "Complexes")
writeData(wb, "Complexes", complex_df)
saveWorkbook(wb, "./data/complex_df_with_p_values_new_new.xlsx", overwrite = TRUE)

# 
# ## ITS GREAT TO SET p critical = 0.05  :>
# ################!!! FILTER CRETERIA III !!!################
# p_critical=0.05
# complex_list=list(
# "eta"=subset(complex_df, SOFT_eta_p<p_critical),
# "epsilon"=subset(complex_df, SOFT_epsilon_p<p_critical),
# "xc"=subset(complex_df, SOFT_xc_p<p_critical),
# "scaling"=subset(complex_df, SOFT_scaling_p<p_critical),
# "xc_epsilon"=subset(complex_df, SOFT_xc_epsilon_p<p_critical)
# )
# 
# 
# 
# ## plot
# library(ggplot2)
# library(ggforce)
# library(ggnewscale)
# plot_complex <-function(complex_selected,name,p_cri){
# ggplot(data.frame(x = seq(0, 2, length.out = 1000)),aes(x=x))+
#     geom_hline(yintercept=1, color = "grey",linewidth = 0.3)+
#     geom_vline(xintercept=1, color = "grey",linewidth = 0.3)+
#     geom_line(aes(y = eta_x(x)),size = 0.3) +
#     geom_line(aes(y = epsilon_x(x)),size = 0.3) +
#     geom_line(aes(y = xc_x(x)),size = 0.3) +
#     geom_point(data=simu_long,aes(x=lifespan_rela,y=steepness_rela,shape=para),
#                size=0.5,alpha=0.5)+
#     scale_shape_manual(values=c(1,2,3), name="Parameters") +
#     new_scale_color() +
#     geom_mark_hull(data=lifes_set_alpha_complex_df[lifes_set_alpha_complex_df$complex %in% complex_selected,],
#                    aes(x=lifespan_rela, y=steepness_rela, 
#                        color=complex,group=complex,fill=complex),
#                    alpha=0.2, colour=NA, expand=0.05,show.legend = FALSE) +
#     geom_point(data=lifes_set_alpha_complex_df[lifes_set_alpha_complex_df$complex %in% complex_selected,],
#                aes(x=lifespan_rela, y=steepness_rela, 
#                    color=complex), alpha=1)+
#     xlim(0,2)+ylim(0.0,2)+
#     coord_fixed(ratio=1)+
#     # theme(legend.position = "none")+
#     labs(x="Lifespan Relative to WT", y="Steepness Relative to WT",
#         title=paste(name,": complexes with p-value < ",p_cri,sep = ""))+
#     theme(
#       # legend.position = "none",
#       panel.grid.major = element_blank(),
#       # panel.grid.minor = element_blank(),
#       panel.background = element_blank(),
#       axis.line = element_line(colour = "black"),
#       plot.title = element_text(size = 10),     # 调整主标题大小
#       plot.subtitle = element_text(size = 8),  # 调整副标题大小
#       axis.title = element_text(size = 8),     # 调整所有坐标轴标题大小
#       axis.text = element_text(size = 5)       # 调整所有坐标轴标签大小
#     )
#   
# }
# 
# plot_complex(complex_list$eta$Complex_ac,"eta",p_critical)
# plot_complex(complex_list$epsilon$Complex_ac,"epsilon",p_critical)
# plot_complex(complex_list$xc$Complex_ac,"Xc",p_critical)
# plot_complex(complex_list$scaling$Complex_ac,"scaling",p_critical)
# # new
# plot_complex(complex_list$xc_epsilon$Complex_ac,"Xc_Epsilon",p_critical)
# 
# 
# 
# 






