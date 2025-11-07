## let's get into annotation files
setwd('C:/Users/20145/Desktop/westlake_summer/database/');
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

## plotting
library(ggplot2)
library(ggforce)
ggplot(lifes_set_alpha_complex_df[lifes_set_alpha_complex_df$complex %in% good_complex,],aes(x=lifespan_rela, y=steepness_rela, 
      color=complex, size=num_of_experi,group=complex,fill=complex)) +
  geom_mark_hull(alpha=0.2, colour=NA, expand=0.02,show.legend = FALSE) +
  geom_point(alpha=1,show.legend=FALSE)+
  labs(x="Lifespan Relative to WT", y="Steepness Relative to WT") +
  # add a title
  ggtitle("Lifespan and Steepness of Genes in Complexes") +
  xlim(0,2)+ylim(0,4)+
  theme_minimal()+
  theme(legend.position = "none")

## add simulation
ggplot()+
  geom_mark_hull(data=lifes_set_alpha_complex_df[lifes_set_alpha_complex_df$complex %in% good_complex,],
                 aes(x=lifespan_rela, y=steepness_rela, 
                     color=complex, size=num_of_experi,group=complex,fill=complex),
                 alpha=0.2, colour=NA, expand=0.025,show.legend = FALSE) +
  geom_point(data=lifes_set_alpha_complex_df[lifes_set_alpha_complex_df$complex %in% good_complex,],
             aes(x=lifespan_rela, y=steepness_rela, 
                 color=complex, size=num_of_experi), alpha=1,show.legend=FALSE)+
  new_scale_color() +
  geom_line(data = simu_long, aes(x = lifespan_rela, y = steepness_rela, color = para, group = para), size = 0.7) +
  # geom_point(data = simu_long, aes(x = lifespan_rela, y = steepness_rela, color = para, group = para, shape = para)) +
  scale_color_brewer(palette = "Set1", name = "Parameters") +
  scale_shape_discrete(guide = "none")+
  xlim(0,2)+ylim(0,2)+
  theme(legend.position = "none")+
  labs(x="Lifespan Relative to WT", y="Steepness Relative to WT",
  title="Lifespan and Steepness of Genes in Complexes with Simulation")+
  theme_minimal()

## zoom into best complexes
selected_complexes=counts_complex$complex[order(-counts_complex$num_genes)][1:4]
lifes_set_alpha_complex_df$description <- sapply(lifes_set_alpha_complex_df$complex, function(x) {
  if (x %in% ComplexPortal_slim$Complex_ac) {
    return(ComplexPortal_slim$Recommended_name[ComplexPortal_slim$Complex_ac == x])
  } else {
    return(NA)
  }
})
ggplot()+
  geom_mark_hull(data=lifes_set_alpha_complex_df[lifes_set_alpha_complex_df$complex %in% selected_complexes,],
                 aes(x=lifespan_rela, y=steepness_rela, 
                     color=complex, size=num_of_experi,group=complex,fill=complex),
                 alpha=0.2, colour=NA, expand=0.025,show.legend = FALSE) +
  geom_point(data=lifes_set_alpha_complex_df[lifes_set_alpha_complex_df$complex %in% selected_complexes,],
             aes(x=lifespan_rela, y=steepness_rela, 
                 color=complex, size=num_of_experi), alpha=1)+
  new_scale_color() +
  geom_line(data = simu_long, aes(x = lifespan_rela, y = steepness_rela, color = para, group = para), size = 1,alpha=1) +
  # geom_point(data = simu_long, aes(x = lifespan_rela, y = steepness_rela, color = para, group = para, shape = para)) +
  scale_color_brewer(palette = "Set2", name = "Parameters") +
  scale_shape_discrete(guide = "none")+
  # xlim(0.2,2)+ylim(0.2,2)+
  theme(legend.position = "none")+
  labs(x="Lifespan Relative to WT", y="Steepness Relative to WT",
       title="Lifespan and Steepness of Genes in Complexes with Simulation (Zoom into Best Complexes)")+
  theme_minimal()

### prepare for enrichment
complex_df=ComplexPortal_slim[ComplexPortal_slim$Complex_ac %in% good_complex,]
complex_df=merge(complex_df,counts_complex,by.x="Complex_ac",by.y="complex",all=FALSE)
complex_df$ratio_genes_exist=as.numeric(complex_df$num_genes)/as.numeric(complex_df$num_all_genes)
################!!! FILTER CRETERIA II !!!################
complex_df=complex_df[complex_df$ratio_genes_exist>0.5,]

############### HARD clustering
distance_short_df=distance_df[,c("gene","eta_rela","epsilon_rela","xc_rela","scaling_rela")]
rownames(distance_short_df)=distance_short_df$gene; distance_short_df=subset(distance_short_df,select=-gene)
distance_short_df$hard_cluster=apply(distance_short_df, 1, function(x) names(distance_short_df)[which.min(x)])
distance_short_df$hard_cluster=gsub("_rela","",distance_short_df$hard_cluster)
distance_df$hard_cluster=distance_short_df$hard_cluster
# see percentage
arr=as.character(distance_short_df$hard_cluster);arr=prop.table(table(arr));label=paste(names(arr),round(100*arr,2),"%")
pie(arr,labels = label,main="percentage of each hard-cluster #genes assigned");rm(arr);rm(label)
# calculate p-value for each complexes > 1

N = length(distance_short_df$hard_cluster) #（总基因数）
K_eta = length(distance_short_df$hard_cluster[distance_short_df$hard_cluster=='eta']) #（属于eta类的基因）
K_epsilon = length(distance_short_df$hard_cluster[distance_short_df$hard_cluster=='epsilon']) 
K_xc = length(distance_short_df$hard_cluster[distance_short_df$hard_cluster=='xc']) 
K_scaling = length(distance_short_df$hard_cluster[distance_short_df$hard_cluster=='scaling']) 

for (i in 1:length(complex_df$Complex_ac)){
  n=complex_df$num_genes[i]
  genes=unlist(complex_df$genes[i])
  temp=na.omit(distance_short_df[genes,])
  k_eta=    length(temp[temp$hard_cluster=='eta',1])
  k_epsilon=length(temp[temp$hard_cluster=='epsilon',1])
  k_xc=     length(temp[temp$hard_cluster=='xc',1])
  k_scaling=length(temp[temp$hard_cluster=='scaling',1])
  complex_df$HARD_eta_p[i]=phyper(q = k_eta - 1, m = K_eta, n = N - K_eta, k = n, lower.tail = FALSE)
  complex_df$HARD_epsilon_p[i]=phyper(q = k_epsilon - 1, m = K_epsilon, n = N - K_epsilon, k = n, lower.tail = FALSE)
  complex_df$HARD_xc_p[i]=phyper(q = k_xc - 1, m = K_xc, n = N - K_xc, k = n, lower.tail = FALSE)
  complex_df$HARD_scaling_p[i]=phyper(q = k_xc - 1, m = K_xc, n = N - K_xc, k = n, lower.tail = FALSE)
}

############### SOFT clustering
distance_weight_df=subset(distance_short_df,select=-hard_cluster)
# distance_weight_df=1/distance_weight_df
# distance_weight_df=as.data.frame(t(apply(distance_weight_df, 1, function(x) x / sum(x, na.rm = TRUE))))
# see the propotion
temp=colSums(distance_weight_df);temp=temp/sum(temp);label=paste(names(temp),round(100*temp,2),"%")
pie(temp,labels = label,main="percentage of sum of weights of genes")
## calculate p-value for each complex
generate_null_dist <- function(n,N=1000,weight_df){
  gene_list=as.character(rownames(weight_df))
  result=matrix(NA,nrow = N,ncol=ncol(weight_df))
  for (i in 1:N){
    gene_null=sample(gene_list, n, replace = FALSE)
    temp=weight_df[gene_null,]
    scores=colSums(temp)
    result[i,]=scores
  }
  colnames(result)=colnames(weight_df)
  result=as.data.frame(result)
  return(result)
}
# temp=generate_null_dist(5,1000,distance_weight_df)
for (i in 1:length(complex_df$Complex_ac)){
  n=complex_df$num_genes[i]
  genes=unlist(complex_df$genes[i])
  temp=na.omit(distance_weight_df[genes,])
  scores=colSums(temp)
  N_sample=2000
  scores_null=generate_null_dist(n,N_sample,distance_weight_df)
  count_vector <- lapply(1:length(scores), function(j) {
    sum(scores_null[[j]] < scores[j])
  })
  count_vector <- unlist(count_vector);names(count_vector)=names(temp)
  complex_df$SOFT_eta_p[i]=(count_vector["eta_rela"]+1)/(N_sample+1)
  complex_df$SOFT_epsilon_p[i]=(count_vector["epsilon_rela"]+1)/(N_sample+1)
  complex_df$SOFT_xc_p[i]=(count_vector["xc_rela"]+1)/(N_sample+1)
  complex_df$SOFT_scaling_p[i]=(count_vector["scaling_rela"]+1)/(N_sample+1)
}

# remove proteins and genes columns
# complex_df=subset(complex_df,select=-c(proteins,genes))

# save into .xlsx
setwd("../data/")
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "Complexes")
writeData(wb, "Complexes", complex_df)
saveWorkbook(wb, "complex_df_with_p_values.xlsx", overwrite = TRUE)
# also save lifes_set_alpha_complex_df
wb2 <- createWorkbook()
addWorksheet(wb2, "lifes_set_alpha_complex")
writeData(wb2, "lifes_set_alpha_complex", lifes_set_alpha_complex_df)
saveWorkbook(wb2, "lifes_set_alpha_complex.xlsx", overwrite = TRUE)

## ITS GREAT TO SET p critical = 0.05  :>
################!!! FILTER CRETERIA III !!!################
p_critical=0.05
complex_list=list(
"eta"=subset(complex_df, SOFT_eta_p<p_critical),
"epsilon"=subset(complex_df, SOFT_epsilon_p<p_critical),
"xc"=subset(complex_df, SOFT_xc_p<p_critical),
"scaling"=subset(complex_df, SOFT_scaling_p<p_critical)
)


## plot
library(ggplot2)
library(ggforce)
library(ggnewscale)
plot_complex <-function(complex_selected,name,p_cri){
ggplot(data.frame(x = seq(0, 2, length.out = 1000)),aes(x=x))+
    geom_hline(yintercept=1, color = "black",linewidth = 0.3)+
    geom_line(aes(y = eta_x(x)),size = 0.3) +
    geom_line(aes(y = epsilon_x(x)),size = 0.3) +
    geom_line(aes(y = xc_x(x)),size = 0.3) +
    geom_point(data=simu_long,aes(x=lifespan_rela,y=steepness_rela,shape=para),size=1,alpha=0.4)+
    scale_shape_manual(values=c(1,2,3))+
    new_scale_color() +
    geom_mark_hull(data=lifes_set_alpha_complex_df[lifes_set_alpha_complex_df$complex %in% complex_selected,],
                   aes(x=lifespan_rela, y=steepness_rela, 
                       color=complex,group=complex,fill=complex),
                   alpha=0.2, colour=NA, expand=0.05,show.legend = FALSE) +
    geom_point(data=lifes_set_alpha_complex_df[lifes_set_alpha_complex_df$complex %in% complex_selected,],
               aes(x=lifespan_rela, y=steepness_rela, 
                   color=complex), alpha=1)+
    xlim(0,2)+ylim(0.2,3.5)+
    theme(legend.position = "none")+
    labs(x="Lifespan Relative to WT", y="Steepness Relative to WT",
        title=paste(name,": complexes with p-value < ",p_cri,sep = ""))+
    theme_minimal()
}

plot_complex(complex_list$eta$Complex_ac,"eta",p_critical)
plot_complex(complex_list$epsilon$Complex_ac,"epsilon",p_critical)
plot_complex(complex_list$xc$Complex_ac,"Xc",p_critical)
plot_complex(complex_list$scaling$Complex_ac,"scaling",p_critical)

## documented into .xlsx
# library(openxlsx)
# wb <- createWorkbook()
# addWorksheet(wb, paste("eta(p<",p_critical,")",sep = ""))
# writeData(wb, paste("eta(p<",p_critical,")",sep = ""), 
#           subset(ComplexPortal, Complex_ac %in% complex_eta$Complex_ac))
# addWorksheet(wb, paste("epsilon(p<",p_critical,")",sep = ""))
# writeData(wb, paste("epsilon(p<",p_critical,")",sep = ""), 
#           subset(ComplexPortal, Complex_ac %in% complex_epsilon$Complex_ac))
# addWorksheet(wb, paste("xc(p<",p_critical,")",sep = ""))
# writeData(wb, paste("xc(p<",p_critical,")",sep = ""), 
#           subset(ComplexPortal, Complex_ac %in% complex_xc$Complex_ac))
# addWorksheet(wb, paste("scaling(p<",p_critical,")",sep = ""))
# writeData(wb, paste("scaling(p<",p_critical,")",sep = ""), 
#           subset(ComplexPortal, Complex_ac %in% complex_scaling$Complex_ac))
# saveWorkbook(wb, "COMPLEXES_WITH_SIGNIFICANCE_PARAMETERS.xlsx", overwrite = TRUE)

################## paste complexes together by GO annotations
# complex_go_df=tibble(ComplexPortal[,c("Complex_ac","Go_Annotations")])
# complex_go_df$Go_Annotations=sapply(complex_go_df$Go_Annotations, function(x){
#   str=as.character(x)
#   vec=strsplit(str,"\\|")[[1]]
#   vec=gsub("\\(.*\\)","",vec)
#   return(list(vec))
# })
# 
# complex_go_df <- complex_go_df %>%
#   unnest(Go_Annotations)



#################### GO annotation
# GO_anno=read_tsv("./sgd.gaf/sgd_without_header.gaf",col_names = FALSE)
# colnames(GO_anno)=c("DB","DB_Object_ID","DB_Object_Symbol","Relation",
#                      "GO_ID","DB_Reference","Evidence_Code","With_or_From",
#                      "Aspect","DB_Object_Name","DB_Object_Synonym",
#                      "DB_Object_Type","Taxon","Date","Assigned_By",
#                      "Annotation_Extension","Gene_Product_Form_ID")
# length(GO_anno[[1]])
# # filter
# evi_codes=c("EXP","IDA","IPI","IMP","IGI","IEP")  # direct experimental evidence
# GO_anno=GO_anno[GO_anno$Evidence_Code %in% evi_codes,]
# length(GO_anno[[1]])
# GO_BP=GO_anno[GO_anno$Aspect=="P",];length(GO_BP[[1]])
# GO_CC=GO_anno[GO_anno$Aspect=="C",];length(GO_CC[[1]])
# GO_MF=GO_anno[GO_anno$Aspect=="F",];length(GO_MF[[1]])

# GO_terms_BP=unique(GO_BP$GO_ID)
# GO_terms_CC=unique(GO_CC$GO_ID)
# GO_terms_MF=unique(GO_MF$GO_ID)


# # BiocManager::install(c("GO.db", "org.Sc.sgd.db"))
# # library(GO.db)
# # library(org.Sc.sgd.db)

# # BiocManager::install("annotate")
# library(annotate)
# # install.packages("GOxploreR")
# library(GOxploreR)

# GOTermBPOnLevel(goterm = GO_terms_BP[102:200])

# library(GO.db)
# GOTERM[["GO:0141164"]]
# ancestors <- GOBPANCESTOR[["GO:0141164"]]
# print(ancestors)

# GO_terms_BP=GOTermBPOnLevel_fixed(goterm = GO_terms_BP)


# GOTermBPOnLevel_fixed <- function(goterm) {
#     if (is.numeric(goterm)) stop("The 'goterm' argument should be a GO-term(s)")
#     x <- goterm
#     ont <- lapply(x, function(y) Ontology(y))
#     isna <- which(is.na(ont))
#     nonretired <- which(ont != "BP")
#     if (length(isna) > 0 || length(nonretired) > 0) {
#         index <- unique(c(isna, nonretired))
#         warning(paste("Check that the terms ", paste(x[index], collapse = ", "), " are BP GO-terms and not obsolete"))
#         x <- x[-index]
#     }
#     if (length(x) > 0) {
#         dat <- data.frame(Term = character(), Level = numeric(), stringsAsFactors = FALSE)
#         go2h <- get("go2h", envir = getNamespace("GOxploreR"))
#         for (i in 1:length(x)) {
#             if (exists(x[i], envir = go2h)) {
#                 level <- go2h[[x[i]]][length(go2h[[x[i]]])] - 1
#                 dat[i, ] <- c(x[i], level)
#             } else {
#                 warning(paste("GO term", x[i], "not found in go2h, skipping"))
#             }
#         }
#         if (nrow(dat) > 0) {
#             colnames(dat) <- c("Term", "Level")
#             return(dat)
#         }
#     }
#     # warning("No valid BP terms processed")
#     return(data.frame(Term = character(), Level = numeric()))
# }

