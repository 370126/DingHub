library(data.table)

parameter_all_df=fread("KVresult_Genotype_uM_EGi_A1_CVcut1_ZY.stat")


abundence_df=fread("TEM1_foldx_stability.stat")
colnames(abundence_df)=c("genotype", "ddG")

parameter_all_df$mutations[parameter_all_df$genotype=="WT"]=0

parameter_nostop_df=parameter_all_df[parameter_all_df$group=="nostop",]
# see mutatants with one mutation
paramter_1_df=parameter_nostop_df[parameter_nostop_df$mutations==1,]

## mutation site to mutation sequence
wt_seq="MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW"

# site2seq <- function(site_vec, wt_seq){
#   site2seq_df=data.frame(site=site_vec, seq=NA)
#   wt_vec=strsplit(wt_seq, split = "")[[1]]
#   for (site in site_vec){
#     sites=strsplit(site, split = ",")[[1]]
#     # each element like "A_123_D"
#     ref_seq=wt_vec
#     for (s in sites){
#       parts=strsplit(s, split = "_")[[1]]
#       wt=parts[1]
#       pos=as.integer(parts[2])
#       mut=parts[3]
#       if (wt!=wt_vec[pos]){
#         stop(paste0("Error: mutation site ", s, " does not match wt sequence!"))
#       }
#       ref_seq[pos]=mut
#     }
#     site2seq_df[site2seq_df$site==site, "seq"]=paste0(ref_seq, collapse = "")
#   }
#   return(site2seq_df)
# }


site2seq_fast <- function(site_vec, wt_seq, show_progress = TRUE) {
  wt_vec <- strsplit(wt_seq, "")[[1]]
  n <- length(site_vec)
  seq_out <- character(n)
  
  parsed_sites <- strsplit(site_vec, ",")

  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = n, style = 3)
  }
  
  for (i in seq_along(parsed_sites)) {
    ref_seq <- wt_vec
    sites <- parsed_sites[[i]]
    if (sites[1]=="WT") {
      seq_out[i] <- wt_seq
      next
    }
    parts <- do.call(rbind, strsplit(sites, "_"))
    wt <- parts[, 1]
    pos <- as.integer(parts[, 2])
    mut <- parts[, 3]
    
    if (any(wt_vec[pos] != wt)) {
      stop(paste0("Error: mutation site(s) ", paste(sites, collapse = ","), " do not match wt sequence!"))
    }
    
    ref_seq[pos] <- mut
    seq_out[i] <- paste0(ref_seq, collapse = "")
    
    if (show_progress && i %% 50 == 0) setTxtProgressBar(pb, i)
  }
  
  if (show_progress) {
    setTxtProgressBar(pb, n)
    close(pb)
  }
  
  data.frame(site = site_vec, seq = seq_out, stringsAsFactors = FALSE)
}


site2seq_df=site2seq_fast(site_vec=parameter_all_df$genotype, wt_seq=wt_seq)
parameter_all_df=merge(parameter_all_df, site2seq_df, by.x="genotype", by.y="site", all.x=TRUE)

site2seq_df=site2seq_fast(site_vec=abundence_df$genotype, wt_seq=wt_seq)
abundence_df=merge(abundence_df, site2seq_df, by.x="genotype", by.y="site", all.x=TRUE)


pymochi_catal_df=data.frame(aa_seq=parameter_all_df$seq, Nham_aa=parameter_all_df$mutations, WT=FALSE, fitness=NA, sigma=NA)

pymochi_catal_df$WT[pymochi_catal_df$aa_seq==wt_seq]=TRUE
V_wt=10^parameter_all_df$V[parameter_all_df$genotype=="WT"]
K_wt=10^parameter_all_df$K[parameter_all_df$genotype=="WT"]
pymochi_catal_df$fitness=log10((10^parameter_all_df$V/V_wt)/(10^parameter_all_df$K/K_wt))

pymochi_catal_df$K=log10(10^parameter_all_df$K/K_wt)
pymochi_catal_df$V=log10(10^parameter_all_df$V/V_wt)



pymochi_catal_df$sigma=1.00

pymochi_catal_nostop_df=pymochi_catal_df[parameter_all_df$group=="nostop",]


# count mutations manually
count_mutations <- function(genotype_vec){
  Nham_vec <- integer(length(genotype_vec))
  # split each element in genotype_vec, then count #
  for (i in seq_along(genotype_vec)){
    genotype <- genotype_vec[i]
    if (genotype == "WT"){
      Nham_vec[i] <- 0
    } else {
      mutations <- strsplit(genotype, ",")[[1]]
      Nham_vec[i] <- length(mutations)
    }
  }
  return(Nham_vec)
}
abundence_df$Nham_aa = count_mutations(abundence_df$genotype)
# add wild type
abundence_df=rbind(abundence_df, data.frame(genotype="WT", ddG=0, seq=wt_seq, Nham_aa=0))


pymochi_stab_df=data.frame(aa_seq=abundence_df$seq, Nham_aa=abundence_df$Nham_aa, WT=FALSE, ddG=abundence_df$ddG, weight=NA)
pymochi_stab_df$WT[pymochi_stab_df$aa_seq==wt_seq]=TRUE




### check the distribution of fitness
library(ggplot2)
p1=ggplot(pymochi_catal_nostop_df, aes(x=fitness))+
  geom_histogram(bins=50)+
  geom_vline(xintercept = 0, color="red", linetype="dashed")+
  theme_bw()+
  xlab("fitness (V/K)")+
  ggtitle("Distribution of fitness (V/K) for all mutants")
print(p1)

# get fitness distribution density function
library(MASS)
# make smoother density curve
fit_density <- density(pymochi_catal_nostop_df$fitness, bw="nrd0")
# visualize histogram and density
p2=ggplot(pymochi_catal_nostop_df, aes(x=fitness))+
  geom_histogram(aes(y=..density..), bins=50, fill="lightblue", alpha=0.7)+
  geom_line(data=data.frame(x=fit_density$x, y=fit_density$y), aes(x=x, y=y), color="blue", size=1)+
  geom_vline(xintercept = 0, color="red", linetype="dashed")+
  geom_hline(yintercept = 0.05, color="red", linetype="dashed")+
  theme_bw()+
  xlab("fitness (V/K)")+
  ggtitle("Fitness (V/K) Distribution with Density Curve")
print(p2)

# check saperately on K & V
fit_density_K <- density(pymochi_catal_nostop_df$K, bw="nrd0")
fit_density_V <- density(pymochi_catal_nostop_df$V, bw="nrd0")
# visualize histogram and density for K
p2_K=ggplot(pymochi_catal_nostop_df, aes(x=K))+
  geom_histogram(aes(y=..density..), bins=50, fill="lightblue", alpha=0.7)+
  geom_line(data=data.frame(x=fit_density_K$x, y=fit_density_K$y), aes(x=x, y=y), color="blue", size=0.5)+
  geom_vline(xintercept = 0, color="red", linetype="dashed")+
  theme_bw()+
  xlab("K (relative to WT)")+
  ggtitle("K (relative to WT) Distribution with Density Curve")
print(p2_K)
p2_V=ggplot(pymochi_catal_nostop_df, aes(x=V))+
  geom_histogram(aes(y=..density..), bins=50, fill="lightblue", alpha=0.7)+
  geom_line(data=data.frame(x=fit_density_V$x, y=fit_density_V$y), aes(x=x, y=y), color="blue", size=0.5)+
  geom_vline(xintercept = 0, color="red", linetype="dashed")+
  theme_bw()+
  xlab("V (relative to WT)")+
  ggtitle("V (relative to WT) Distribution with Density Curve")
print(p2_V)


## set weight of each mutant using density function
# interpolate density values for each fitness
density_values <- approx(fit_density$x, fit_density$y, xout = pymochi_catal_nostop_df$fitness)$y
# avoid zero density values
density_values[density_values == 0] <- min(density_values[density_values > 0])
# set density>0.1 all to 1
density_values[density_values > 0.05] <- 0.2
# set weight as modified density
pymochi_catal_nostop_df$weight <- 5*density_values


# fwrite(pymochi_catal_nostop_df, file="pymochi_catal_nostop_dataset_weighted.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)





### check the distribution of ddG
p3=ggplot(pymochi_stab_df, aes(x=ddG))+
  geom_histogram(bins=50)+
  geom_vline(xintercept = 0, color="red", linetype="dashed")+
  theme_bw()+
  xlab("ddG (kcal/mol)")+
  ggtitle("Distribution of ddG for all mutants")
print(p3)
# get ddG distribution density function
fit_density_ddG <- density(pymochi_stab_df$ddG, bw="nrd0")
# visualize histogram and density
p4=ggplot(pymochi_stab_df, aes(x=(ddG)))+
  geom_histogram(aes(y=..density..), bins=50, fill="lightblue", alpha=0.7)+
  geom_line(data=data.frame(x=fit_density_ddG$x, y=fit_density_ddG$y), aes(x=x, y=y), color="blue", size=1)+
  geom_vline(xintercept = 0, color="red", linetype="dashed")+
  geom_vline(xintercept = 40, color="red", linetype="dashed")+
  theme_bw()+
  xlab("ddG (kcal/mol)")+
  ggtitle("ddG Distribution with Density Curve")
print(p4)

pymochi_stab_df$weight = 1.0
pymochi_stab_df$weight[pymochi_stab_df$ddG>40]=5/pymochi_stab_df$ddG[pymochi_stab_df$ddG>40]
# fwrite(pymochi_stab_df, file="pymochi_stab_dataset_weighted.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


# COMBINE STABILITY AND CATALYTIC DATASETS
# merge by aa_seq
pymochi_combined_df=merge(pymochi_catal_nostop_df, pymochi_stab_df[, c("aa_seq", "ddG", "weight")], by="aa_seq", suffixes = c("_catal", "_stab"))
# fwrite(pymochi_combined_df, file="pymochi_combined_dataset_weighted.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

pymochi_combined_ddG_df=pymochi_combined_df[,c("aa_seq", "Nham_aa","WT","ddG","weight_stab")]
colnames(pymochi_combined_ddG_df)[5]="sigma"; colnames(pymochi_combined_ddG_df)[4]="fitness"
fwrite(pymochi_combined_ddG_df, file="data_ddG.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

pymochi_combined_K_df=pymochi_combined_df[,c("aa_seq", "Nham_aa","WT","K","weight_catal")]
colnames(pymochi_combined_K_df)[5]="sigma"; colnames(pymochi_combined_K_df)[4]="fitness"
# fwrite(pymochi_combined_K_df, file="data_K.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

########### check the correlation between ddG and V ###########
parameter_merged_df=merge(parameter_all_df, abundence_df[, .(genotype, ddG)], by="genotype", all.x=TRUE)
parameter_merged_df=parameter_merged_df[,-"V1"]


parameter_merged_less_df=parameter_merged_df[mutations<=9 & mutations>0, ]

parameter_merged_less_df$mutations=as.factor(parameter_merged_less_df$mutations)

p5=ggplot(parameter_merged_less_df, aes(x=ddG, y=V, color=mutations))+
  geom_point(alpha=0.1, size=0.5, stroke=0)+
  geom_smooth(method="loess", color="red", linewidth=0.5, se=TRUE)+
  theme_bw()+
  coord_cartesian(xlim=c(-10, 60), ylim=c(-2.5, 2.5))+
  facet_wrap(~mutations, nrow=3)+
  xlab("ddG (kcal/mol)")+
  ylab("log10(V relative to WT)")+
  # remove legend
  theme(legend.position = "none")+
  ggtitle("Correlation between ddG and V")
# print(p5)
# ggsave(p5, filename = "figures/ddG_vs_V_correlation.png", width = 6, height = 5, dpi = 300)
# ggsave(p5, filename = "figures/ddG_vs_V_correlation_saperate.png", width = 10, height = 8, dpi = 500)

# K v.s. ddG
p6=ggplot(parameter_merged_less_df, aes(x=ddG, y=K, color=mutations))+
  geom_point(alpha=0.1, size=0.5, stroke=0)+
  geom_smooth(method="lm", color="red", linewidth=0.5)+
  theme_bw()+
  coord_cartesian(xlim=c(-10, 60), ylim=c(-2.5, 3))+
  # facet_wrap(~mutations, nrow=4)+
  xlab("ddG (kcal/mol)")+
  ylab("log10(K relative to WT)")+
  # remove legend
  theme(legend.position = "none")+
  ggtitle("Correlation between ddG and K")
# print(p6)
# ggsave(p6, filename = "figures/ddG_vs_K_correlation.png", width = 10, height = 8, dpi = 500)


# K v.s V
p7=ggplot(parameter_merged_less_df, aes(x=V, y=K, color=mutations))+
  geom_point(alpha=0.1, size=0.5, stroke=0)+
  # geom_smooth(method="lm", color="red", linewidth=0.5, se=FALSE)+
  theme_bw()+
  # coord_cartesian(xlim=c(-2.5, 2.5), ylim=c(-2.5, 3))+
  # facet_wrap(~mutations, nrow=3)+
  xlab("log10(V relative to WT)")+
  ylab("log10(K relative to WT)")+
  # remove legend
  theme(legend.position = "none")+
  ggtitle("Correlation between V and K")
# print(p7)
# ggsave(p7, filename = "figures/V_vs_K_correlation_saperate.png", width = 10, height = 8, dpi = 500)
# ggsave(p7, filename = "figures/V_vs_K_correlation.png", width = 6, height = 5, dpi = 300)




### check distribution of mutations
p=ggplot(data=parameter_all_df, aes(x=mutations))+
  geom_bar(fill="lightblue", color="black")+
  geom_vline(xintercept = 10, color="red", linetype="dashed")+
  theme_bw()+
  xlab("Number of mutations")+
  ylab("Count")+
  ggtitle("Distribution of number of mutations in all genotypes")
# print(p)
# ggsave(p, filename = "figures/mutation_number_distribution.png", width=6, height=4, dpi=300)

## abundence
parameter_merged_less_df$abundence=1/(1+exp((parameter_merged_less_df$ddG-61)/(8.314*298)))
# distribution
p9=ggplot(parameter_merged_less_df, aes(x=abundence))+
  geom_histogram(bins=50, fill="lightblue", color="black")+
  geom_vline(xintercept = 0.5, color="red", linetype="dashed")+
  theme_bw()+
  xlab("Abundence (estimated from ddG)")+
  ggtitle("Distribution of abundence for all mutants")
# print(p9)


# V v.s. abundence
p8=ggplot(parameter_merged_less_df, aes(x=abundence, y=V, color=mutations))+
  geom_point(alpha=0.1, size=0.5, stroke=0)+
  geom_smooth(method="lm", color="red", linewidth=0.5)+
  theme_bw()+
  # coord_cartesian(xlim=c(0, 1), ylim=c(-2.5, 2.5))+
  facet_wrap(~mutations, nrow=4)+
  xlab("Abundence (estimated from ddG)")+
  ylab("log10(V relative to WT)")+
  # remove legend
  theme(legend.position = "none")+
  ggtitle("Correlation between abundence and V")
# print(p8)


##################### clutering mapping #####################
residue_cluster_df=read.csv("residue_features_clustering.csv")
res_0=residue_cluster_df$Residue_ID[residue_cluster_df$Cluster=="0"]

res_core=c(68,70,71,72,69,67)

# filter mutants with mutations only in res_core
is_core_mutation <- function(genotype, core_residues) {
  if (genotype == "WT") {
    return(FALSE)
  }
  mutations <- strsplit(genotype, ",")[[1]]
  for (mut in mutations) {
    parts <- strsplit(mut, "_")[[1]]
    pos <- as.integer(parts[2])
    if (!(pos %in% core_residues)) {
      return(FALSE)
    }
  }
  return(TRUE)
}
parameter_core_df=parameter_nostop_df[sapply(parameter_nostop_df$genotype, is_core_mutation, core_residues = res_core), ]





