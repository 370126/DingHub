## import data
setwd('C:/Users/20145/Desktop/westlake_summer/');
getwd()

library(openxlsx)
library(readxl)
yeast_life_clean <- read.xlsx("cleaner_yeastdata (YPD,30C,20+).xlsx")
yeast_life_clean$set_genotype <- gsub("\\*", "", yeast_life_clean$set_genotype)
# remove NA in set_mating_type
yeast_life_clean <- yeast_life_clean[!is.na(yeast_life_clean$set_mating_type), ]
yeast_life_clean$set_genotype=toupper(yeast_life_clean$set_genotype) # convert to uppercase


print(paste("Number of unique genotypes:", length(unique(yeast_life_clean$set_genotype))))

simu_result=list()
sheet_names <- excel_sheets("simu_result.xlsx")
for (sheet in sheet_names) {
  simu_result[[sheet]] <- read.xlsx("simu_result.xlsx", sheet = sheet, 
                                    colNames = TRUE, rowNames = TRUE)
}
simu_result$lifespan_rela=as.data.frame(lapply(simu_result$lifespan_absolute,function(x) 
  x/(simu_result$lifespan_absolute$`0`)))
simu_result$lifespan_rela$para=rownames(simu_result$lifespan_absolute)
simu_result$steepness_rela=as.data.frame(lapply(simu_result$steepness,function(x) 
  x/(simu_result$steepness$`0`)))
simu_result$steepness_rela$para=rownames(simu_result$steepness)
# NEW: SKEWNESS
simu_result$skewness_rela=as.data.frame(lapply(simu_result$skewness,function(x) 
  x/(simu_result$skewness$`0`)))
simu_result$skewness_rela$para=rownames(simu_result$skewness)
# combine into long dataframe for ggplot
library(tidyverse)
simu_X=pivot_longer(simu_result$lifespan_rela,cols = -para, names_to = "eff", values_to = "lifespan_rela")
simu_Y=pivot_longer(simu_result$steepness_rela,cols = -para, names_to = "eff", values_to = "steepness_rela")
simu_Z=pivot_longer(simu_result$skewness_rela,cols = -para, names_to = "eff", values_to = "skewness_rela")
simu_long <- simu_X %>% 
  left_join(simu_Y, by = c("para", "eff")) %>%
  left_join(simu_Z, by = c("para", "eff"))
rm(simu_X);rm(simu_Y);rm(simu_Z)
write.csv(simu_long, "./data/simu_Long.csv", row.names = FALSE)


library(tidyverse)
## check the data
allGenes <- unique(yeast_life_clean$set_genotype)
# Count the number of experiments and mating types for each gene
counts <- map_df(allGenes, function(gene) {
  total_exp <- sum(yeast_life_clean$set_genotype == gene,na.rm = TRUE)
  matalpha <- sum(yeast_life_clean$set_genotype == gene & yeast_life_clean$set_mating_type == "matalpha",na.rm = TRUE)
  mata <- sum(yeast_life_clean$set_genotype == gene & yeast_life_clean$set_mating_type == "MATa",na.rm = TRUE)
  tibble(gene = gene, exp = total_exp, matalpha = matalpha, MATa = mata)
})

## classify the genes
# "For each strain that showed a mean RLS increase of >30% over control, or 
# p < 0.05 for increased RLS, we measured RLS for 20 cells in the MATa strain 
# carrying the same gene deletion."
coupleGenes <- counts %>% filter(matalpha > 0) %>% pull(gene)  # ROI
deepGenes <- counts %>% filter(MATa == 0 & matalpha > 1) %>% pull(gene)
singleGenes <- counts %>% filter(MATa == 0 & matalpha == 1) %>% pull(gene)
coupleDeepGenes <- counts %>% filter(MATa > 0 & matalpha > 1) %>% pull(gene)
coupleShallowGenes <- counts %>% filter(MATa > 0 & matalpha == 1) %>% pull(gene)
catCounts <- c(length(deepGenes), length(singleGenes), length(coupleGenes), 
               length(coupleDeepGenes), length(coupleShallowGenes))
catLabels <- c("more than 1 MATalpha experiment but no MATa",
               "only 1 MATalpha experiment and no MATa",
               "both mating types",
               "both mating types and more than 1 MATalpha",
               "both mating types and only 1 MATalpha")
walk2(catCounts, catLabels, ~cat(.x, " genes with ", .y, "\n"))
gc()

## deal with replication
yeast_roi=yeast_life_clean[yeast_life_clean$set_genotype %in% coupleGenes,]
yeast_roi_alpha=yeast_roi[yeast_roi$set_mating_type=="matalpha",]
yeast_roi_a=yeast_roi[yeast_roi$set_mating_type=="MATa",]
yeast_NOT_roi=yeast_life_clean[!yeast_life_clean$set_genotype %in% coupleGenes,]

counts_alpha <- map_df(coupleGenes, function(gene) {
  exp <- sum(yeast_roi_alpha$set_genotype == gene)
  # calculate the mean RLS for each gene
  mean_rls <- mean(yeast_roi_alpha$set_lifespan_mean[
    yeast_roi_alpha$set_genotype == gene])
  mean_rls_rela <- mean(yeast_roi_alpha$set_lifespan_mean[
    yeast_roi_alpha$set_genotype == gene]/yeast_roi_alpha$ref_lifespan_mean[
      yeast_roi_alpha$set_genotype == gene])
  std_rls <- sd(yeast_roi_alpha$set_lifespan_mean[
    yeast_roi_alpha$set_genotype == gene])
  std_rls_rela <- sd(yeast_roi_alpha$set_lifespan_mean[
    yeast_roi_alpha$set_genotype == gene]/yeast_roi_alpha$ref_lifespan_mean[
      yeast_roi_alpha$set_genotype == gene])
  tibble(gene = gene, num_of_experi = exp, mean_rls=mean_rls, 
         mean_rls_rela = mean_rls_rela, std_mean_rls = std_rls,
         std_mean_rls_rela=std_rls_rela)
})
counts_alpha$cv_mean_rls_rela=counts_alpha$std_mean_rls_rela/counts_alpha$mean_rls_rela

counts_a <- map_df(coupleGenes, function(gene) {
  exp <- sum(yeast_roi_a$set_genotype == gene)
  # calculate the mean RLS for each gene
  mean_rls <- mean(yeast_roi_a$set_lifespan_mean[
    yeast_roi_a$set_genotype == gene])
  mean_rls_rela <- mean(yeast_roi_a$set_lifespan_mean[
    yeast_roi_a$set_genotype == gene]/yeast_roi_a$ref_lifespan_mean[
      yeast_roi_a$set_genotype == gene])
  std_rls <- sd(yeast_roi_a$set_lifespan_mean[
    yeast_roi_a$set_genotype == gene])
  std_rls_rela <- sd(yeast_roi_a$set_lifespan_mean[
    yeast_roi_a$set_genotype == gene]/yeast_roi_a$ref_lifespan_mean[
      yeast_roi_a$set_genotype == gene])
  tibble(gene = gene, num_of_experi = exp, mean_rls=mean_rls, 
         mean_rls_rela = mean_rls_rela, std_mean_rls = std_rls,
         std_mean_rls_rela=std_rls_rela)
})
counts_a$cv_mean_rls_rela=counts_a$std_mean_rls_rela/counts_a$mean_rls_rela
gc()

## Plot the number of experiments for each gene
library(ggplot2)
ggplot(counts_alpha, aes(x = num_of_experi)) +
  geom_bar(fill = "steelblue") +
  labs(title = "Number of Experiments for Each Gene in MATalpha",
       x = "Number of Experiments",
       y = "Count") +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5,size=1.5) +
  theme_minimal()

ggplot(counts_a, aes(x = num_of_experi)) +
  geom_bar(fill = "steelblue") +
  labs(title = "Number of Experiments for Each Gene in MATa",
       x = "Number of Experiments",
       y = "Count") +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5,size=1.5) +
  theme_minimal()

# filter ?
# counts_alpha=counts_alpha[counts_alpha$num_of_experi>1 & counts_alpha$cv_mean_rls_rela<0.3,]
# length(counts_alpha$gene) # genes with more than 1 experiment and cv < 0.3
# counts_a=counts_a[counts_a$num_of_experi>1 & counts_a$cv_mean_rls_rela<0.3,]
# length(counts_a$gene)


yeast_roi_alpha=yeast_roi_alpha[yeast_roi_alpha$set_genotype %in% counts_alpha$gene,] # nolint
yeast_roi_a=yeast_roi_a[yeast_roi_a$set_genotype %in% counts_a$gene,]

# calculate the steepness of the lifespans
cal_steep <- function(str_col) {
  str_vec= as.character(str_col)
  steepness_vec=c();
  for (i in 1:length(str_vec)){
    str=str_vec[i]
    vec=strsplit(str, ",")[[1]]
    vec=as.numeric(vec)
    vec=sort(vec,decreasing=TRUE)
    l=length(vec)
    vec_s=vec
    steepness_vec[i]=1/(sd(vec_s)/mean(vec_s))
  }
  return(steepness_vec)
}

yeast_roi_alpha$set_steepness=cal_steep(yeast_roi_alpha$set_lifespans)
yeast_roi_alpha$ref_steepness=cal_steep(yeast_roi_alpha$ref_lifespans)
yeast_roi_alpha$steepness_rela=yeast_roi_alpha$set_steepness/yeast_roi_alpha$ref_steepness # nolint
yeast_roi_alpha$lifespan_mean_rela=yeast_roi_alpha$set_lifespan_mean/yeast_roi_alpha$ref_lifespan_mean # nolint

yeast_roi_a$set_steepness=cal_steep(yeast_roi_a$set_lifespans)
yeast_roi_a$ref_steepness=cal_steep(yeast_roi_a$ref_lifespans)
yeast_roi_a$steepness_rela=yeast_roi_a$set_steepness/yeast_roi_a$ref_steepness # nolint
yeast_roi_a$lifespan_mean_rela=yeast_roi_a$set_lifespan_mean/yeast_roi_a$ref_lifespan_mean # nolint
gc()

## Plot the relative steepness of lifespans v.s. the relative lifespan mean
ggplot(yeast_roi_alpha, aes(x = lifespan_mean_rela, y = steepness_rela,color=set_genotype )) +
  geom_hline(yintercept=1, color = "grey",linewidth = 0.3)+
  geom_vline(xintercept=1, color = "grey",linewidth = 0.3)+
  geom_point(size=1,alpha=0.5) +
  labs(title = "Relative Steepness of Lifespans vs. Relative Lifespan Mean (MATalpha)",
       x = "Relative Lifespan Mean",
       y = "Relative Steepness") +
  theme_minimal()+
  ylim(0,4)+xlim(0,4)+
  # coord_fixed()+
  # remove the legend
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(size = 10),     # 调整主标题大小
    plot.subtitle = element_text(size = 8),  # 调整副标题大小
    axis.title = element_text(size = 8),     # 调整所有坐标轴标题大小
    axis.text = element_text(size = 5)       # 调整所有坐标轴标签大小
  )

############zoom into most focused genes
top20_genes=counts_alpha$gene[order(-counts_alpha$num_of_experi)][1:20]
top10_genes=counts_alpha$gene[order(-counts_alpha$num_of_experi)][1:10]
selected_genes=c("los1", "sir2", "hap4",'los1','fob1','rpn4') # zoom into these genes
selected_genes=toupper(selected_genes) # convert to uppercase
############

# zoom into top 20 genes
ggplot(yeast_roi_alpha[yeast_roi_alpha$set_genotype %in% top20_genes,], aes(x = lifespan_mean_rela, y = steepness_rela,color=set_genotype )) +
  geom_point(size=1,alpha=0.5) +
  labs(title = "Relative Steepness of Lifespans vs. Relative Lifespan Mean (MATalpha)",
       x = "Relative Lifespan Mean",
       y = "Relative Steepness") +
  theme_minimal()+
  theme(legend.position = "none")
ggplot(yeast_roi_alpha[yeast_roi_alpha$set_genotype %in% selected_genes,], aes(x = lifespan_mean_rela, y = steepness_rela,color=set_genotype )) +
  geom_point() +
  labs(title = "Relative Steepness of Lifespans vs. Relative Lifespan Mean (MATalpha)",
       x = "Relative Lifespan Mean",
       y = "Relative Steepness")+
  theme_minimal()
gc()


## plot MC
ggplot(data=simu_long,aes(x=lifespan_rela,y=steepness_rela,color=para,group=para,shape=para))+
  #geom_smooth(se=FALSE)+
  geom_line(size=0.5,alpha=0.75)+
  geom_point()+
  theme_minimal()+
  labs(x = "relative lifespan", y = "relative steepness",
       title='Simulation of parameter experiments')

# add experimental data
library(ggnewscale)
ggplot()+
  geom_point(data=yeast_roi_alpha,aes(x=lifespan_mean_rela, y=steepness_rela,color=set_genotype,
                                      shape = NULL, group = NULL),size=1,alpha=0.2,show.legend = FALSE)+
  new_scale_color() +
  geom_line(data=simu_long,aes(x=lifespan_rela,y=steepness_rela,group=para),size=0.5)+
  geom_point(data=simu_long,aes(x=lifespan_rela,y=steepness_rela,shape=para),size=0.5,alpha=0.5)+
  scale_shape_manual(values=c(1,2,3))+
  #geom_smooth(se=FALSE)+
  scale_color_brewer(palette = "Set1", name = "Parameters") +
  xlim(0.2, 2.0) +ylim(0.2, 4.0) +
  theme_minimal()+
  labs(x = "relative lifespan", y = "relative steepness",
       title='Simulation of parameter experiments & experimental data')
#theme(legend.position = 'none')

# zoom into top20 genes
ggplot()+
  geom_point(data=yeast_roi_alpha[yeast_roi_alpha$set_genotype %in% top20_genes,],aes(x=lifespan_mean_rela, y=steepness_rela,color=set_genotype,
                                                                                      shape = NULL, group = NULL),size=1,alpha=0.5,show.legend=FALSE)+
  new_scale_color() +
  geom_line(data=simu_long,aes(x=lifespan_rela,y=steepness_rela,color=para,group=para),size=0.5)+
  geom_point(data=simu_long,aes(x=lifespan_rela,y=steepness_rela,color=para,group=para,shape=para))+
  #geom_smooth(se=FALSE)+
  scale_color_brewer(palette = "Set1", name = "Parameters") +
  scale_shape_discrete(guide = "none")+
  xlim(0.2, 2.0) +
  theme_minimal()+
  labs(x = "relative lifespan", y = "relative steepness",
       title='Simulation of parameter experiments & experimental data (zoom into 20 genes)')
#theme(legend.position = 'none')

# zoom into 3 genes
selected_genes=c("SOD1","SOD2")
p1=ggplot()+
  geom_hline(yintercept=1, color = "black",linewidth = 0.7,alpha=0.5)+
  geom_point(data=yeast_roi_alpha[yeast_roi_alpha$set_genotype %in% selected_genes,],
             aes(x=lifespan_mean_rela, y=steepness_rela,color=set_genotype,
                 shape = NULL, group = NULL),size=2,alpha=0.25,stroke=0)+
  # geom_smooth(data=yeast_roi_alpha[yeast_roi_alpha$set_genotype %in% selected_genes,],
  #             aes(x=lifespan_mean_rela, y=steepness_rela,group=set_genotype, color=set_genotype),
  #             se=TRUE, size=0.5, alpha=0.1,method = "gam") +
  new_scale_color() +
  geom_line(data=simu_long,aes(x=lifespan_rela,y=steepness_rela,group=para),size=0.5)+
  geom_point(data=simu_long,aes(x=lifespan_rela,y=steepness_rela,color=para,group=para,shape=para))+
  scale_color_brewer(palette = "Set1", name = "Parameters") +
  scale_shape_discrete(guide = "none")+
  # xlim(0.2, 2.0) +
  theme_minimal()+
  labs(x = "relative lifespan", y = "relative steepness",
       title='Simulation of parameter experiments & experimental data (zoom into SOD)');p1



# integration -------------------------------------------------------------


#########################################  reorganize the data
merge_lifes_array <- function (df, col_name) {
  lifes_array=numeric() 
  col=df[[col_name]];
  l=length(col)
  for (i in 1:l) {
    seq=col[[i]]
    seq=as.numeric(unlist(strsplit(seq,",")))
    # add to the array
    lifes_array= c(lifes_array,seq)
  }
  return(lifes_array)
}
lifes_WT_alpha=merge_lifes_array(yeast_roi_alpha,"ref_lifespans")
lifes_WT_a=merge_lifes_array(yeast_roi_a,"ref_lifespans")
lifes_WT_not=merge_lifes_array(yeast_NOT_roi,"ref_lifespans")
length(lifes_WT_alpha);length(lifes_WT_a);length(lifes_WT_not)
# plot the distribution of WT lifespans of avove three on single graph
library(ggplot2)
ggplot(data = data.frame(lifespan = c(lifes_WT_alpha, lifes_WT_a, lifes_WT_not), 
                         type = rep(c("MATalpha", "MATa", "NOT ROI"), 
                                    times = c(length(lifes_WT_alpha), length(lifes_WT_a), length(lifes_WT_not)))), 
       aes(x = lifespan, fill = type)) +
  geom_histogram(binwidth = 2, position = "identity", alpha = 0.5, color = "black") +
  labs(title = "Distribution of WT Lifespan",
       x = "Lifespan",
       y = "Count") +
  theme_minimal() +
  scale_fill_manual(values = c("steelblue", "orange", "lightgreen")) +
  theme(legend.title = element_blank())


mean_lifespan_WT_alpha=mean(yeast_roi_alpha$ref_lifespan_mean)
vec_s=lifes_WT_alpha %>% 
  sort(decreasing = TRUE)
steepness_WT_alpha=1/(sd(vec_s)/mean(vec_s))
# NEW: calculate skewness
skew <- function(vec){
  n=length(vec)
  s=sd(vec)
  mean_vec=mean(vec)
  vec_3=((vec-mean_vec)/s)^3
  g1=(n/((n-1)*(n-2)))*sum(vec_3)
  return(g1)
}
skewness_WT_alpha=skew(lifes_WT_alpha)

genes_alpha=unique(yeast_roi_alpha$set_genotype)
lifes_set_alpha=list()
lifes_set_alpha_length=numeric()
for (gene in genes_alpha) {
  lifes_set_alpha[[gene]]=merge_lifes_array(yeast_roi_alpha[yeast_roi_alpha$set_genotype==gene,], "set_lifespans")
  lifes_set_alpha_length[gene]=length(lifes_set_alpha[[gene]])
}
# plot the distribution of number of lifes in for each gene
ggplot(data = data.frame(length = lifes_set_alpha_length), aes(x = length)) +
  geom_histogram(binwidth = 50, fill = "steelblue", color = "black")+
  labs(title = "Distribution of #Lifespan for Each Gene (MATalpha)",
       x = "Number of Lifespan Measurements",
       y = "Count") +
  theme_minimal()

lifes_set_alpha_df <- map_df(genes_alpha, function(gene) {
  vec <- lifes_set_alpha[[gene]]
  num=length(vec)
  vec <- sort(vec, decreasing = TRUE)
  vec_s <- vec
  tibble(
    gene = gene,
    steepness = 1/(sd(vec_s) / mean(vec_s)),
    mean = mean(vec_s),
    num=num
  )
})
lifes_set_alpha_df$steepness_rela=lifes_set_alpha_df$steepness/steepness_WT_alpha
lifes_set_alpha_df$mean_rela=lifes_set_alpha_df$mean/mean_lifespan_WT_alpha
# plot the steepness and mean of each gene
ggplot(lifes_set_alpha_df, aes(x = mean_rela, y = steepness_rela)) +
  geom_point(aes(color = gene), size = 1, alpha = 1) +
  labs(title = "Steepness vs. Mean Lifespan for Each Gene (MATalpha)",
       x = "Relative Mean Lifespan",
       y = "Relative Steepness") +
  theme_minimal() +
  theme(legend.position = "none")
# zoom into 20 genes
ggplot(lifes_set_alpha_df[lifes_set_alpha_df$gene %in% top20_genes, ], 
       aes(x = mean_rela, y = steepness_rela, color = gene)) +
  geom_point(size = 1, alpha = 1) +
  labs(title = "Steepness vs. Mean Lifespan for Selected Genes (MATalpha)",
       x = "Relative Mean Lifespan",
       y = "Relative Steepness") +
  theme_minimal() +
  theme(legend.position = "none")
# zoom into 3 genes
ggplot(lifes_set_alpha_df[lifes_set_alpha_df$gene %in% selected_genes, ], 
       aes(x = mean_rela, y = steepness_rela, color = gene)) +
  geom_point(size = 2.5, alpha = 1) +
  labs(title = "Steepness vs. Mean Lifespan for Selected Genes (MATalpha)",
       x = "Relative Mean Lifespan",
       y = "Relative Steepness") +
  theme_minimal() +
  theme(legend.position = "none")

# add the simulation result
ggplot() +
  geom_point(data = lifes_set_alpha_df, aes(x = mean_rela, y = steepness_rela, color = gene), 
             size = 1, alpha = 0.5,show.legend = FALSE) +
  new_scale_color() +
  geom_line(data = simu_long, aes(x = lifespan_rela, y = steepness_rela, color = para, group = para), size = 0.7) +
  geom_point(data = simu_long, aes(x = lifespan_rela, y = steepness_rela, color = para, group = para, shape = para)) +
  scale_color_brewer(palette = "Set1", name = "Parameters") +
  scale_shape_discrete(guide = "none") +
  xlim(0.2, 2.0) +
  theme_minimal() +
  labs(x = "Relative Mean Lifespan", y = "Relative Steepness",
       title='Simulation of parameter experiments & experimental data (MATalpha)')
# zoom into 20 genes
ggplot() +
  geom_point(data = lifes_set_alpha_df[lifes_set_alpha_df$gene %in% top20_genes, ], 
             aes(x = mean_rela, y = steepness_rela, color = gene), 
             size = 2.5, alpha = 1) +
  new_scale_color() +
  geom_line(data = simu_long, aes(x = lifespan_rela, y = steepness_rela, color = para, group = para), size = 0.7) +
  geom_point(data = simu_long, aes(x = lifespan_rela, y = steepness_rela, color = para, group = para, shape = para)) +
  scale_color_brewer(palette = "Set1", name = "Parameters") +
  scale_shape_discrete(guide = "none") +
  xlim(0.2, 2.0) +
  theme_minimal() +
  labs(x = "Relative Mean Lifespan", y = "Relative Steepness",
       title='Simulation of parameter experiments & experimental data (zoom into 3 genes)')

# zoom into 3 genes
p2=ggplot() +
  geom_hline(yintercept=1, color = "black",linewidth = 0.7,alpha=0.5)+
  geom_point(data = lifes_set_alpha_df[lifes_set_alpha_df$gene %in% selected_genes, ], 
             aes(x = mean_rela, y = steepness_rela, color = gene), 
             size = 3, alpha = 1) +
  new_scale_color() +
  geom_line(data = simu_long, aes(x = lifespan_rela, y = steepness_rela, color = para, group = para), size = 0.7) +
  geom_point(data = simu_long, aes(x = lifespan_rela, y = steepness_rela, color = para, group = para, shape = para)) +
  scale_color_brewer(palette = "Set1", name = "Parameters") +
  scale_shape_discrete(guide = "none") +
  # xlim(0.2, 2.0) +
  theme_minimal() +
  labs(x = "Relative Mean Lifespan", y = "Relative Steepness",
       title='Simulation of parameter experiments & experimental data (zoom into 3 genes)');p2


# bootstrap ---------------------------------------------------------------
## bootstrap the steepness and mean lifespan of each gene
# install.packages("MASS")
gc()
lifes_set_alpha_bootstrap_matrix <- list()
for (i in 1:length(lifes_set_alpha)){
  gene=names(lifes_set_alpha)[i]
  vec=lifes_set_alpha[[gene]]
  num=length(vec)
  bootstrap=5000
  bootstrap_matrix <- matrix(NA, nrow = bootstrap, ncol = num)
  for (j in 1:bootstrap) {
    bootstrap_matrix[j, ] <- sample(vec, num, replace = TRUE)
  }
  lifes_set_alpha_bootstrap_matrix[[gene]] <- bootstrap_matrix
}


lifes_set_alpha_bootstrap_df=tibble()
for (gene in names(lifes_set_alpha_bootstrap_matrix)) {
  bootstrap_matrix <- lifes_set_alpha_bootstrap_matrix[[gene]]
  steepness_rela_vec <- apply(bootstrap_matrix, 1, function(x) {
    vec_s <- sort(x, decreasing = TRUE)
    1/(sd(vec_s) / mean(vec_s))
  })
  steepness_rela_vec <- steepness_rela_vec / steepness_WT_alpha
  # NEW: SKEWNESS
  skewness_rela_vec <- apply(bootstrap_matrix,1,function(x){
    skew(x)
  })
  skewness_rela_vec <- skewness_rela_vec / skewness_WT_alpha
  mean_rela_vec <- rowMeans(bootstrap_matrix)
  mean_rela_vec <- mean_rela_vec / mean_lifespan_WT_alpha
  # Store the results in the data frame
  lifes_set_alpha_bootstrap_df <- rbind(lifes_set_alpha_bootstrap_df,
                                         tibble(gene = gene,
                                                steepness_rela_bootstrap = list(steepness_rela_vec),
                                                mean_rela_bootstrap = list(mean_rela_vec),
                                                skewness_rela_bootstrap=list(skewness_rela_vec)))
}
names(lifes_set_alpha_bootstrap_df)

lifes_set_alpha_bp_df=data.frame()
lifes_set_alpha_bp_df <- lifes_set_alpha_bootstrap_df %>%
  dplyr::select(gene, 
         steepness_rela_bootstrap, 
         mean_rela_bootstrap, 
         skewness_rela_bootstrap) %>%
  tidyr::unnest(c(steepness_rela_bootstrap, 
           mean_rela_bootstrap, 
           skewness_rela_bootstrap))
# # save into .csv file
# write.csv(lifes_set_alpha_bp_df, "lifes_set_alpha_bootstrap_df.csv", row.names = FALSE)
# 
# # load into lifes_set_alpha_bp_df
# lifes_set_alpha_bp_df <- read.csv("./lifes_set_alpha_bootstrap_df.csv")
head(lifes_set_alpha_bp_df)
# plot the bootstrap result using epllipse
# ggplot(lifes_set_alpha_bp_df, aes(x = mean_rela_bootstrap, y = steepness_rela_bootstrap,
#                                   colour = gene, group=gene)) +
#   geom_point(aes(color = gene), size = 0.1, alpha = 0.1) +
#   stat_ellipse(level = 0.95, alpha = 0.5,type="norm",linetype=1,
#                linewidth = 0.5) +
#   labs(title = "Bootstrap Steepness vs. Mean Lifespan for Each Gene (MATalpha)",
#        x = "Relative Mean Lifespan",
#        y = "Relative Steepness") +
#   theme_minimal() +
#   coord_fixed()+
#   theme(legend.position = "none")

# zoom into top 20 genes
ggplot(lifes_set_alpha_bp_df[lifes_set_alpha_bp_df$gene %in% top20_genes, ], 
       aes(x = mean_rela_bootstrap, y = steepness_rela_bootstrap, color = gene, group=gene)) +
  # geom_point(size = 0.01, alpha = 0.01, stroke=0, color="grey") +
  stat_ellipse(level = 0.95, alpha = 1,type="norm",
               linewidth = 0.7) +
  labs(title = "Bootstrap Steepness vs. Mean Lifespan for Selected Genes (MATalpha)",
       x = "Relative Mean Lifespan",
       y = "Relative Steepness") +
  theme_minimal() +
  coord_fixed()+
  theme(legend.position = "none")

# zoom into 3 genes
ggplot(lifes_set_alpha_bp_df[lifes_set_alpha_bp_df$gene %in% selected_genes, ], 
       aes(x = mean_rela_bootstrap, y = steepness_rela_bootstrap, color = gene, group=gene)) +
  geom_point(size = 0.1, alpha = 0.5,stroke=0) +
  stat_ellipse(level = 0.95, alpha = 1,type="norm",
               linewidth = 0.5) +
  labs(title = "Bootstrap Steepness vs. Mean Lifespan for Selected Genes (MATalpha)",
       x = "Relative Mean Lifespan",
       y = "Relative Steepness") +
  theme_minimal() +
  coord_fixed()

# plot ellipse with simulation result
# ggplot()+
#   # simulation data
#   geom_point(data=simu_long,aes(x=lifespan_rela,y=steepness_rela,color=para,group=para,shape=para),size=0.5)+
#   geom_line(data=simu_long,aes(x=lifespan_rela,y=steepness_rela,color=para,group=para),size=0.5)+
#   scale_color_brewer(palette = "Set1", name = "Parameters") +
#   new_scale_color() +
#   # bootstrap data
#   geom_point(data=lifes_set_alpha_bp_df, aes(x = mean_rela_bootstrap, y = steepness_rela_bootstrap, 
#                                              color = gene, group=gene), size = 0.1, alpha = 0.1,stroke=0) +
#   stat_ellipse(data=lifes_set_alpha_bp_df, aes(x = mean_rela_bootstrap, y = steepness_rela_bootstrap, 
#                                                color = gene, group=gene),
#                level = 0.95, alpha = 0.75,type="norm",linewidth=0.3) +
#   
#   scale_shape_discrete(guide = "none") +
#   coord_fixed() +
#   theme_minimal() +
#   labs(x = "Relative Mean Lifespan", y = "Relative Steepness",
#        title='Simulation of parameter experiments & bootstrap experimental data (MATalpha)') +
#   theme(legend.position = "none")

# zoom into top 20 genes
ggplot()+
  # simulation data
  geom_point(data=simu_long,aes(x=lifespan_rela,y=steepness_rela,color=para,group=para,shape=para),size=0.5)+
  geom_line(data=simu_long,aes(x=lifespan_rela,y=steepness_rela,color=para,group=para),size=0.5)+
  scale_color_brewer(palette = "Set1", name = "Parameters") +
  new_scale_color() +
  # bootstrap data
  geom_point(data=lifes_set_alpha_bp_df[lifes_set_alpha_bp_df$gene %in% top20_genes, ], 
             aes(x = mean_rela_bootstrap, y = steepness_rela_bootstrap, group=gene), 
             size = 0.01, alpha = 0.01, color="grey" ,stroke=0,show.legend=FALSE) +
  stat_ellipse(data=lifes_set_alpha_bp_df[lifes_set_alpha_bp_df$gene %in% top20_genes, ], 
               aes(x = mean_rela_bootstrap, y = steepness_rela_bootstrap, color = gene, group=gene),
               level = 0.95, alpha = 1,type="norm",
               linewidth = 0.5,show.legend=FALSE) +
  # scale_color_npg()+
  xlim(0.2, 2.0) + ylim(0.2, 2) +
  scale_shape_discrete(guide = "none") +
  coord_fixed() +
  theme_minimal() +
  labs(x = "Relative Mean Lifespan", y = "Relative Steepness",
       title='Simulation of parameter experiments & bootstrap experimental data (zoom into top20 genes)')
       
# zoom into 3 genes
p3=ggplot()+
  geom_hline(yintercept=1, color = "black",linewidth = 0.7,alpha=0.5)+
  # simulation data
  geom_point(data=simu_long,aes(x=lifespan_rela,y=steepness_rela,color=para,group=para,shape=para),size=0.5)+
  geom_line(data=simu_long,aes(x=lifespan_rela,y=steepness_rela,color=para,group=para),size=0.5)+
  scale_color_brewer(palette = "Set1", name = "Parameters") +
  new_scale_color() +
  # bootstrap data
  geom_point(data=lifes_set_alpha_bp_df[lifes_set_alpha_bp_df$gene %in% selected_genes, ], 
             aes(x = mean_rela_bootstrap, y = steepness_rela_bootstrap, color = gene, group=gene), 
             size = 0.5, alpha = 0.1,stroke=0) +
  stat_ellipse(data=lifes_set_alpha_bp_df[lifes_set_alpha_bp_df$gene %in% selected_genes, ], 
               aes(x = mean_rela_bootstrap, y = steepness_rela_bootstrap, color = gene, group=gene),
               level = 0.95, alpha = 1,type="norm",
               linewidth = 0.5) +
  scale_shape_discrete(guide = "none") +
  xlim(0.2, 2.0) + ylim(0.2, 2) +
  coord_fixed() +
  theme_minimal() +
  labs(x = "Relative Mean Lifespan", y = "Relative Steepness",
       title='Simulation of parameter experiments & bootstrap experimental data (zoom into 3 genes)');p3

