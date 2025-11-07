## package loading
required_packages <- c("tidyverse", "ggplot2", "ggsci", "reshape2", "reticulate","e1071")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
np <- import("numpy")
rm(required_packages, pkg)

## import meta-information from python .npz files
path="D:/westlake_ZJ/GLM_downstream/seed_dna_minid_050_v3/seed_dna_minid_050_v3/"
# item_name=list.files(path=path,full.names = FALSE)
# item=data.frame(
#   name=item_name,
#   ID=gsub("\\_.*","",item_name)
# )
# item$CDS = mapply(gsub, item$ID, "", item$name) %>%
#   gsub("\\_\\(0\\,","",.) %>%
#   gsub("\\)","",.) %>%
#   as.numeric()
# item$CDS_token = ceiling(item$CDS/6)
# rownames(item)=item$name
# # save
# write.csv(item,"./data/item_CDS.csv",row.names = FALSE)
# # read
item=read.csv("./data/item_CDS.csv")


############ GLOBAL VARIABLES ############
rownames(item)=item$name
item_test=item$name[1000];print(item_test)
CDS_test=item[item_test,"CDS_token"];print(CDS_test)
# boundary check range
zoom_check_upper=c(CDS_test-50,CDS_test-5)
zoom_check_lower=c(CDS_test+3,CDS_test+25)
zoom_check=c(zoom_check_upper[1],zoom_check_lower[2])
# window size for f6
windows_size_upper=6
windows_size_lower=6
# suggested weights for the 6 functions
weight=c(0,0,1,0,1,1)  # f1,f2,f3,f4,f5,f6
##########################################



fetch_array <- function(item,layer_id) {
  file_name=paste0(item,"_attention_layer_",layer_id,".npz")
  path_full=paste(path,item,paste0("layer_",layer_id),file_name,sep="/")
  data <- np$load(path_full, allow_pickle = TRUE)
  attention_array <- py_to_r(data[["attention"]])
  return(attention_array)
}

plot_array <- function(item,layer_id,head_id) {
  attention_array=fetch_array(item = item,layer_id=layer_id)
  array <- attention_array[head_id, , ]
  array [upper.tri(array)] <- NA
  # array=t(array);# array[lower.tri(array)] <- NA
  df <- melt(array)
  # replace 0 with a small value to avoid issues with log scale
  # df$value[df$value == 0] <- 0
  colnames(df) <- c("Query", "Key", "Weight")
  ggplot(df, aes(x = Key, y = Query, fill = Weight)) +
    geom_tile(na.rm = FALSE) +
    scale_fill_viridis_c(trans = "log10", na.value = "transparent") +
    # add CDS line
    geom_vline(xintercept = CDS_test, color = "red",linewidth=0.8) +
    geom_hline(yintercept = CDS_test, color = "red",linewidth=0.8) +
    # add text label besides the cross of the lines
    annotate("text", x = CDS_test + 5, y = CDS_test - 5, label = "CDS end", color = "red", hjust = 0, vjust = 0, size=3.5) +
    coord_equal(expand = FALSE) +
    scale_x_continuous(expand = c(0,0), breaks = seq(0, ncol(array), by=20)) +
    scale_y_reverse(expand = c(0,0), breaks = seq(0, nrow(array), by=20)) +
    labs(title = paste(gsub("\\_.*","",item),": ","Layer",layer_id, "|", "Attention Head",head_id), 
         x = "Key position", y = "Query position", fill = "Attention") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = "black"), 
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(), 
      axis.text = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),
      plot.margin = margin(0,0,0,0),
      # remove xy axis labels
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
}
plot_array(item_test,layer_id=20,head_id=1)

# library(gridExtra)
# plots <- lapply(1:32, function(i) plot_array(item_test,20,i))
# 
# grid.arrange(grobs = plots[1:4], ncol = 2, nrow = 2)
# grid.arrange(grobs = plots[5:8], ncol = 2, nrow = 2)
# grid.arrange(grobs = plots[9:12], ncol = 2, nrow = 2)
# grid.arrange(grobs = plots[13:16], ncol = 2, nrow = 2)
# grid.arrange(grobs = plots[17:20], ncol = 2, nrow = 2)
# grid.arrange(grobs = plots[21:24], ncol = 2, nrow = 2)
# grid.arrange(grobs = plots[25:28], ncol = 2, nrow = 2)
# grid.arrange(grobs = plots[29:32], ncol = 2, nrow = 2)

# # save all head
# ggsave("./figures/test_attention_all_heads.png",
#        grid.arrange(grobs = plots, ncol = 4, nrow = 8),
#        width = 16, height = 32, dpi = 400)
## judging
# cal_entropy <- function(mat_h) {
#   mat[upper.tri(mat)] <- NA
#   vals <- mat_h[!is.na(mat_h)]
#   p <- vals / sum(vals)
#   H <- -sum(ifelse(p > 0, p * log(p), 0))
#   return(H)
# }

## functions to calculate scores
f1 <- function(df, n) {
  inter_df= subset(df,X <= n & Y > n)
  var=var(inter_df$value)
  score=(1+1/var)^(-1)  # for better visualization
  return(score)
}
f2 <- function(df, n) {
  protein_func_df= subset(df, X <= n & Y <= n)
  rna_2_cds_df = subset(df, X > n & Y > n)
  inter_df= subset(df,X <= n & Y > n)
  score=mean(c(protein_func_df$value, rna_2_cds_df$value))/(mean(inter_df$value)+1e-10)
  score=score  # for better visualization
  return(score)
}
f3 <- function(df, n) {
  protein_func_df= subset(df, X <= n & Y <= n)
  rna_2_cds_df = subset(df, X > n & Y > n)
  inter_df= subset(df,X <= n & Y > n)
  score=abs(mean(c(protein_func_df$value, rna_2_cds_df$value)) - mean(inter_df$value))/(mean(c(protein_func_df$value, rna_2_cds_df$value)) + mean(inter_df$value)+1e-10)
  return(score)
}
f4 <- function(df, n) {
  inter_df= subset(df,X <= n & Y > n)
  # top 3% inter value
  inter_top=quantile(inter_df$value,0.97)
  inter_top=mean(inter_df$value[inter_df$value>inter_top])
  score=1-mean(inter_top)/max(inter_df$value)
  return(score)
}
f5 <- function(df, n) {
  inter_df= subset(df,X <= n & Y > n)
  score=1-mean(inter_df$value)/max(inter_df$value)
  return(score)
}
f6 <- function(df, n) {
  if (n<=CDS_test) {
    window_size=windows_size_upper
  } else {
    window_size=windows_size_lower
  }
  n_left=max(1,n-window_size)
  n_right=min(max(df$X),n+window_size)
  window_left_df=subset(df,X > n_left & X <= n & Y > n_left & Y <= n)
  window_right_df=subset(df,X > n & X <= n_right & Y > n & Y <= n_right)
  gradient=abs(mean(window_right_df$value)-mean(window_left_df$value))/window_size
  return(gradient)
}
function_list=list(f1=f1,f2=f2,f3=f3,f4=f4,f5=f5,f6=f6)


cal_score <- function(item,layer_id,head_id){
  attention_array <- fetch_array(item = item_test, layer_id =  layer_id)
  mat=attention_array[head_id,,]
  mat[upper.tri(mat)] <- NA
  # log transform
  mat=log10(mat)
  df=melt(mat);colnames(df) <- c("Y", "X", "value");df=na.omit(df)
  df_upper=subset(df,X<=CDS_test & Y<=CDS_test)
  df_lower=subset(df,X>CDS_test & Y>CDS_test)
  n_upper=length(unique(df_upper$X))
  n_lower=length(unique(df_lower$X))
  function_names=names(function_list)
  # initialize scores data frame
  scores_upper=data.frame(n=1:n_upper)
  scores_upper[function_names]=NA
  scores_lower=data.frame(n=(n_upper+1):(n_upper+n_lower))
  scores_lower[function_names]=NA
  # calculate scores for each function
  for (n in 1:n_upper) {
    for (fname in function_names) {
      scores_upper[n,fname]=function_list[[fname]](df_upper,n)
    }
  }
  for (n in 1:n_lower) {
    for (fname in function_names) {
      scores_lower[n,fname]=function_list[[fname]](df_lower,n+n_upper)
    }
  }
  # plot scores together
  scores_together=rbind(scores_upper,scores_lower)
  return(scores_together)
}

norm_separate <- function(df,CDS) {
  # this function is only for visualization purpose
  # normalize to 0-1 for each column, separately for upper and lower
  df[1:CDS,2:ncol(df)]=apply(df[1:CDS,2:ncol(df)],2,function(x){
    (x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE))
  })
  df[(CDS+1):nrow(df),2:ncol(df)]=apply(df[(CDS+1):nrow(df),2:ncol(df)],2,function(x){
    (x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE))
  })
  return(df)
}


plot_score <- function(scores_together_raw,CDS,head_id=NA,layer_id=NA,zoom=c(NA,NA)){
  temp=scores_together_raw
  # normalize to 0-1 for each column, separately for upper and lower
  temp=norm_separate(temp,CDS)
  scores_melt=melt(temp,id.vars = "n")
  colnames(scores_melt)=c("n","function","score")
  p=ggplot(scores_melt,aes(x=n,y=score,color=`function`))+
    # shade CDS region
    annotate("rect", xmin = 0, xmax = CDS, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "grey50") +
    # shade check region
    annotate("rect", xmin = zoom_check_upper[1], xmax = zoom_check_upper[2], ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "blue") +
    annotate("rect", xmin = zoom_check_lower[1], xmax = zoom_check_lower[2], ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "green") +
    geom_vline(xintercept = CDS, color = "grey10",linewidth=0.8) +
    geom_line(size=0.8)+
    annotate("text", x = CDS - 25, y = 0.97*max(na.omit(scores_melt$score)), label = "CDS", color = "grey40", hjust = 0, vjust = 0, size=4) +
    labs(title = paste(gsub("\\_.*","",item_test),": ","Layer",layer_id, "|", "Attention Head",head_id), 
         x = "Position", y = "Score") +
    scale_x_continuous(expand = c(0,0), breaks = seq(0, (nrow(temp)), by=10)) +
    scale_color_bmj() +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = "black"), 
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(), 
      axis.text = element_text(color = "black"),
      axis.ticks = element_line(color = "black",size=0.8),
      plot.margin = margin(0,0,0,0)
    )
  if (!any(is.na(zoom))) {
    p=p+coord_cartesian(xlim = c(zoom[1]-5, zoom[2]+5))
  }
  return(p)
}

add_weight <- function(df,weight) {
  df_nor=norm_separate(df,CDS_test)
  weighted_score=df_nor %>%
    select(-n) %>%
    as.matrix(.) %*% weight
  weighted_score=as.vector(weighted_score)
  return(weighted_score)
}


find_peak <- function(x,min_height=0.1,neighbor=2){
  peaks=c()
  for (i in 3:(length(x)-2)) {
    # peak should be greater than its neighbors at least 0.2
    condition_a=
    if ((x[i]-x[i-neighbor])>min_height & (x[i]-x[i+neighbor])>min_height) {
      peaks=c(peaks,i)
    }
  }
  peaks=as.data.frame(cbind(peaks,x[peaks]))
  colnames(peaks)=c("position","score")
  return(peaks)
}

trimmed_cv <- function(x, trim = 0.1) {
  x=na.omit(x)
  x_sorted <- sort(x)
  n <- length(x_sorted)
  k <- floor(n * trim)
  x_trim <- x_sorted[(k+1):(n-k)]
  cv <- sd(x_trim) / mean(x_trim)
  return(cv)
}


# test
# head=1, layer=22 shows a clear pattern
# head=1, layer=25 is a bad example
# head=1, layer=26 is a weirdo
# layer_test=23
# head_test=20
# plot_array(item_test,layer_id=layer_test,head_id=head_test)
# temp=cal_score(item_test,layer_id=layer_test, head_id=head_test)
# temp$weighted_score=add_weight(temp,weight=weight)
# # find boundary
# upper_peaks=find_peak(temp$weighted_score[zoom_check_upper[1]:zoom_check_upper[2]])
# lower_peaks=find_peak(temp$weighted_score[zoom_check_lower[1]:zoom_check_lower[2]])
# # use max peak as boundary
# upper_boundary=upper_peaks$position[which.max(upper_peaks[,"score"])]+zoom_check_upper[1]-1
# lower_boundary=lower_peaks$position[which.max(lower_peaks[,"score"])]+zoom_check_lower[1]-1
# # upper_boundary=which.max(temp$weighted_score[zoom_check_upper[1]:zoom_check_upper[2]])+zoom_check_upper[1]-1
# # lower_boundary=which.max(temp$weighted_score[zoom_check_lower[1]:zoom_check_lower[2]])+zoom_check_lower[1]-1
# print(paste("Upper boundary is at",upper_boundary));print(paste("Lower boundary is at",lower_boundary))
# 
# 
# temp_window_upper=temp$weighted_score[zoom_check_upper[1]:zoom_check_upper[2]]
# temp_window_lower=temp$weighted_score[zoom_check_lower[1]:zoom_check_lower[2]]
# # TRIMMED CV
# # the bigger the trimmed CV, the better the pattern
# cv_upper=trimmed_cv(temp$weighted_score[1:CDS_test], trim = 0.1)
# cv_lower=trimmed_cv(temp$weighted_score[CDS_test:nrow(temp)], trim = 0.1)
# print(paste("Upper trimmed CV is",round(cv_upper,4)));print(paste("Lower trimmed CV is",round(cv_lower,4)))
# 
# p=plot_score(temp,CDS=CDS_test,head_id=head_test ,layer_id=layer_test)
# p=p+ geom_vline(xintercept = upper_boundary, color = "blue", linetype = "dashed",linewidth=0.8) +
#   geom_vline(xintercept = lower_boundary, color = "green", linetype = "dashed",linewidth=0.8)
# print(p)
# 
# ## according to the literature, upper boundary of reRNA are approximately 120bp(20 tokens) away from CDS end,
# ## and lower boundary of reRNA are approximately 30bp(5 tokens) away from CDS end.
# p=plot_score(temp,CDS=CDS_test,head_id=head_test,layer_id=layer_test,zoom=zoom_check)+
#   geom_vline(xintercept = upper_boundary, color = "blue", linetype = "dashed",linewidth=0.8) +
#   geom_vline(xintercept = lower_boundary, color = "green", linetype = "dashed",linewidth=0.8)
# print(p)
# # see array
# p=plot_array(item_test,layer_id=layer_test,head_id=head_test)
# p=p+ geom_vline(xintercept = upper_boundary, color = "blue", linetype = "dashed",linewidth=0.8) +
#   geom_hline(yintercept = upper_boundary, color = "blue", linetype = "dashed",linewidth=0.8) +
#   geom_vline(xintercept = lower_boundary, color = "green", linetype = "dashed",linewidth=0.8) +
#   geom_hline(yintercept = lower_boundary, color = "green", linetype = "dashed",linewidth=0.8) +
#   annotate("text", x = upper_boundary + 5, y = upper_boundary - 5, label = "Upper boundary", color = "blue", hjust = 0, vjust = 0, size=3.5) +
#   annotate("text", x = lower_boundary + 5, y = lower_boundary - 5, label = "Lower boundary", color = "green", hjust = 0, vjust = 0, size=3.5)
# print(p)


# simply check all all layers & heads
upper_index <- data.frame(
  layer = rep(13:26, each = 32),
  head = rep(1:32, times = 14),
  trimmed_cv_all = NA,
  max = NA,
  peak = NA,
  max_chu_mean = NA,
  peak_chu_mean = NA,
  max_times_trimmed_cv = NA,
  peak_times_trimmed_cv = NA,
  peak_times_trimmed_cv_chu_mean = NA
)
row_idx <- 1
for (layer_id in 13:26) {
  for (head_id in 1:32) {
    print(paste("Processing Layer", layer_id, "Head", head_id))
    temp = cal_score(item_test, layer_id = layer_id, head_id = head_id)
    temp$weighted_score = add_weight(temp, weight)
    temp_window_upper = temp$weighted_score[zoom_check_upper[1]:zoom_check_upper[2]]
    cv_upper = trimmed_cv(temp$weighted_score[1:CDS_test], trim = 0.1)
    max_upper = max(temp_window_upper, na.rm = TRUE)
    peaks_upper = find_peak(temp_window_upper)
    if (nrow(peaks_upper) > 0) {
      peak_upper = peaks_upper$position[which.max(peaks_upper[,"score"])] + zoom_check_upper[1] - 1
      peak_score = peaks_upper$score[which.max(peaks_upper[,"score"])]
    } else {
      peak_upper = NA
      peak_score = NA
    }
    mean_temp_window = mean(temp_window_upper, na.rm = TRUE)
    upper_index[row_idx, "trimmed_cv_all"] = cv_upper
    upper_index[row_idx, "max"] = max_upper
    upper_index[row_idx, "peak"] = peak_upper
    upper_index[row_idx, "max_chu_mean"] = max_upper / mean_temp_window
    upper_index[row_idx, "peak_chu_mean"] = peak_score / mean_temp_window
    upper_index[row_idx, "max_times_trimmed_cv"] = max_upper * cv_upper
    upper_index[row_idx, "peak_times_trimmed_cv"] = peak_score * cv_upper
    upper_index[row_idx, "peak_times_trimmed_cv_chu_mean"] = peak_score * cv_upper / mean_temp_window
    row_idx <- row_idx + 1
  }
}
# remove temporary variables
rm(temp,temp_window_upper,cv_upper,max_upper,peak_upper,
   max_chu_mean,peak_chu_mean,max_times_trimmed_cv,
   peak_times_trimmed_cv,peak_times_trimmed_cv_chu_mean,peaks_upper)

# save
write.csv(upper_index,"./data/upper_boundary_index_all_heads.csv",row.names = FALSE)
upper_index=read.csv("./data/upper_boundary_index_all_heads.csv")
upper_index$rank=rank(-upper_index$peak_times_trimmed_cv_chu_mean)

## see the best 6 together
best6=upper_index[order(upper_index$rank),][1:6,]
# first see array
plots_array_best6=list()
for (i in 1:nrow(best6)) {
  print(paste("Processing",i,"of",nrow(best6)))
  layer_id=best6$layer[i]
  head_id=best6$head[i]
  p=plot_array(item_test,layer_id=layer_id,head_id=head_id)
  temp=cal_score(item_test,layer_id=layer_id, head_id=head_id)
  temp$weighted_score=add_weight(temp,weight=weight)
  upper_peaks=find_peak(temp$weighted_score[zoom_check_upper[1]:zoom_check_upper[2]])
  upper_boundary=upper_peaks$position[which.max(upper_peaks[,"score"])]+zoom_check_upper[1]-1
  p=p+ geom_vline(xintercept = upper_boundary, color = "blue", linetype = "dashed",linewidth=0.8) +
    geom_hline(yintercept = upper_boundary, color = "blue", linetype = "dashed",linewidth=0.8) +
    annotate("text", x = upper_boundary + 5, y = upper_boundary - 5, label = "Upper boundary", color = "blue", hjust = 0, vjust = 0, size=3.5)
  plots_array_best6[[paste0("L",layer_id,"H",head_id)]]=p
}
library(gridExtra)
p1=grid.arrange(grobs = plots_array_best6, ncol = 2, nrow = 3)
print(p1)
# then see score
plots_score_best6=list()
for (i in 1:nrow(best6)) {
  layer_id=best6$layer[i]
  head_id=best6$head[i]
  temp=cal_score(item_test,layer_id=layer_id, head_id=head_id)
  temp$weighted_score=add_weight(temp,weight=weight)
  upper_peaks=find_peak(temp$weighted_score[zoom_check_upper[1]:zoom_check_upper[2]])
  upper_boundary=upper_peaks$position[which.max(upper_peaks[,"score"])]+zoom_check_upper[1]-1
  p=plot_score(temp,CDS=CDS_test,head_id=head_id,layer_id=layer_id,zoom=zoom_check)+
    geom_vline(xintercept = upper_boundary, color = "blue", linetype = "dashed",linewidth=0.8) +
    annotate("text", x = upper_boundary + 5, y = 0.97*max(na.omit(temp$weighted_score)), label = "Upper boundary", color = "blue", hjust = 0, vjust = 0, size=4)
  plots_score_best6[[paste0("L",layer_id,"H",head_id)]]=p
}
p2=grid.arrange(grobs = plots_score_best6, ncol = 2, nrow = 3)
print(p2)
# save both
ggsave("./figures/best6_upper_boundary_array.png",
       p1,
       width = 20, height = 18, dpi = 400)
ggsave("./figures/best6_upper_boundary_score.png",
       p2,
       width = 20, height = 18, dpi = 400)



# # check all heads
# all_scores=list()
# plots_list=list()
# for (layer_id in 13:26) {
#   for (head_id in 1:32) {
#     print(paste("Processing Layer",layer_id,"Head",head_id))
#     temp=cal_score(item_test,layer_id=layer_id, head_id=head_id)
#     temp=add_weight(temp,weight)
#     p=plot_score(temp,CDS=CDS_test,head_id=head_id,layer_id=layer_id)
#     plots_list[[paste0("L",layer_id,"H",head_id)]]=p
#     temp$layer=layer_id
#     temp$head=head_id
#     all_scores[[paste0("L",layer_id,"H",head_id)]]=temp
#   }
# }
# # store
# all_scores_df=do.call(rbind,all_scores)
# write.csv(all_scores_df,"./data/attention_scores_all_heads_test.csv",row.names = FALSE)
# # save all plots
# library(gridExtra)
# for (layer_id in 13:26) {
#   plots_sub=plots_list[grep(paste0("L",layer_id,"H"),names(plots_list))]
#   ggsave(paste0("./figures/test_attention_scores_layer_",layer_id,"_all_heads.png"),
#          grid.arrange(grobs = plots_sub, ncol = 4, nrow = 8),
#          width = 32, height = 32, dpi = 400)
# }




# all_scores_df=read.csv("./data/attention_scores_all_heads_test.csv")
# # shrink layer and head into 1 column
# all_scores_df$layer_head=paste0("L",all_scores_df$layer,"H",all_scores_df$head)
# all_scores_df=all_scores_df %>%
#   select(-layer,-head)
# # calculate index for each graph
# graph_index_upper_df=data.frame(
#   graph=unique(all_scores_df$layer_head)
# )
# graph_index_upper_df[names(function_list)]=NA
# graph_index_lower_df=graph_index_upper_df

# # calculate graph index from score sequence
# cal_graph_index <- function(df) {
#   if (max(df$n)<=CDS_test) {
#     boundary=zoom_check[1]
#     df_ROI=df[df$n>=boundary,]
#   } else {
#     boundary=zoom_check[2]
#     df_ROI=df[df$n<=boundary,]
#   }
#   temp=df_ROI[,names(function_list)]
  
#   index = c()
#   for (col in names(temp)) {
#     vals = na.omit(temp[[col]])
#     vals_10_percentile = quantile(vals, 0.95)
#     top_vals = vals[vals >= vals_10_percentile]
#     left_vals = vals[vals < vals_10_percentile]
#     idx= (mean(top_vals)-mean(left_vals))/(var(left_vals)/mean(left_vals))
#     index[col] = 1/idx
#   }
#   # index=1/index
#   # dim(index)
#   # [1] 1 6
#   return(index)
# }

# for(i in unique(all_scores_df$layer_head)){
#   df=all_scores_df[all_scores_df$layer_head==i,]
#   df_upper=df[df$n<=CDS_test,]
#   df_lower=df[df$n>CDS_test,]
#   # upper
#   graph_index_upper_df[graph_index_upper_df$graph==i,2:ncol(graph_index_upper_df)]=cal_graph_index(df_upper)
#   # lower
#   graph_index_lower_df[graph_index_lower_df$graph==i,2:ncol(graph_index_lower_df)]=cal_graph_index(df_lower)
# }

# # document the rank of each function in each graph
# rank_graph_index <- function(df) {
#   df_rank=df
#   df_rank[2:ncol(df_rank)]=apply(df[2:ncol(df)],2,rank)
#   return(df_rank)
# }
# graph_rank_upper_df=rank_graph_index(graph_index_upper_df)
# graph_rank_lower_df=rank_graph_index(graph_index_lower_df)

# # weighting again
# # weight=c(0,0,0,0,0,1)
# graph_rank_upper_df$weighted_rank <- graph_rank_upper_df[,names(function_list)] %>%
#   as.matrix(.) %*% weight
# graph_rank_lower_df$weighted_rank <- graph_rank_lower_df[,names(function_list)] %>%
#   as.matrix(.) %*% weight
# graph_rank_upper_df$weighted_rank_rank=rank(graph_rank_upper_df$weighted_rank)
# graph_rank_lower_df$weighted_rank_rank=rank(graph_rank_lower_df$weighted_rank)

# # fetch best layer and head ids
# best_upper=graph_rank_upper_df[which.min(graph_rank_upper_df$weighted_rank_rank),"graph"]
# layer_best_upper=as.numeric(gsub("L(.*)H.*","\\1",best_upper))
# head_best_upper=as.numeric(gsub("L.*H(.*)","\\1",best_upper))
# print(paste("Best upper graph is",best_upper))

# # see the plot
# best_lower=graph_rank_lower_df[which.min(graph_rank_lower_df$weighted_rank_rank),"graph"]
# layer_best_lower=as.numeric(gsub("L(.*)H.*","\\1",best_lower))
# head_best_lower=as.numeric(gsub("L.*H(.*)","\\1",best_lower))
# print(paste("Best lower graph is",best_lower))

# # see the plot
# p=plot_array(item_test,layer_best_upper,head_best_upper)
# temp=all_scores_df[all_scores_df$layer_head==best_upper,]
# upper_boundary=which.max(temp$weighted_score[zoom_check_upper[1]:zoom_check_upper[2]])+zoom_check_upper[1]-1
# lower_boundary=which.max(temp$weighted_score[zoom_check_lower[1]:zoom_check_lower[2]])+zoom_check_lower[1]-1
# p=p+ geom_vline(xintercept = upper_boundary, color = "blue", linetype = "dashed",linewidth=0.8) +
#   geom_hline(yintercept = upper_boundary, color = "blue", linetype = "dashed",linewidth=0.8) +
#   geom_vline(xintercept = lower_boundary, color = "green", linetype = "dashed",linewidth=0.8) +
#   geom_hline(yintercept = lower_boundary, color = "green", linetype = "dashed",linewidth=0.8) +
#   annotate("text", x = upper_boundary + 5, y = upper_boundary - 5, label = "Upper boundary", color = "blue", hjust = 0, vjust = 0, size=3.5) +
#   annotate("text", x = lower_boundary + 5, y = lower_boundary - 5, label = "Lower boundary", color = "green", hjust = 0, vjust = 0, size=3.5)
# print(p)
# temp=temp %>%
#   select(-layer_head)

# p=plot_score(temp,CDS=CDS_test,head_id=head_best_upper,layer_id=layer_best_upper)+
#   geom_vline(xintercept = upper_boundary, color = "blue", linetype = "dashed",linewidth=0.8) +
#   geom_vline(xintercept = lower_boundary, color = "green", linetype = "dashed",linewidth=0.8)
# print(p)

