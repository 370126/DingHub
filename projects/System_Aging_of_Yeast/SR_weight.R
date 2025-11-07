# Load necessary libraries
library(data.table)
distance_all_bp_df <- fread("./data/distance_all_bp_df.csv", sep="," ,header=TRUE, stringsAsFactors = FALSE)
distance_all_bp_df=as.data.frame(distance_all_bp_df)
head(distance_all_bp_df)
distance_df = read.csv("./data/distance_df.csv", header = TRUE, stringsAsFactors = FALSE)
head(distance_df)


# calculate SD of each gene, each parameter across bootstrap samples
library(dplyr)
distance_sd_df <- distance_all_bp_df %>%
  group_by(gene) %>%
  summarise(
    eta_sd = sd(eta, na.rm = TRUE),
    epsilon_sd = sd(epsilon, na.rm = TRUE),
    xc_sd = sd(xc, na.rm = TRUE),
    scaling_sd = sd(scaling, na.rm = TRUE)
  )
head(distance_sd_df)
distance_df = merge(distance_df, distance_sd_df, by="gene")
head(distance_df)

# for each gene, calculate their minimal 2 parameter distance
distance_df$d2_minus_d1 <- NA
distance_df$sqrt_sd1_2_plus_sd2_2 <- NA
distance_df$weight_new <- NA

for (i in seq_len(nrow(distance_df))) {
  # collect distances and SDs for this row
  distances <- c(distance_df$eta[i], distance_df$epsilon[i], distance_df$xc[i], distance_df$scaling[i])
  sd_values <- c(distance_df$eta_sd[i], distance_df$epsilon_sd[i], distance_df$xc_sd[i], distance_df$scaling_sd[i])
  # order distances, dropping NAs
  ord <- order(distances, na.last = NA)
  if (length(ord) < 2) {
    # fewer than 2 valid distances -> cannot compute d2-d1 or combined SD
    distance_df$d2_minus_d1[i] <- NA_real_
    distance_df$sqrt_sd1_2_plus_sd2_2[i] <- NA_real_
    next
  }
  idx1 <- ord[1]
  idx2 <- ord[2]
  d1 <- distances[idx1]
  d2 <- distances[idx2]
  distance_df$d2_minus_d1[i] <- d2 - d1
  # get corresponding SDs by index (robust to duplicate values and unnamed vectors)
  sd1 <- sd_values[idx1]
  sd2 <- sd_values[idx2]
  if (is.na(sd1) || is.na(sd2)) {
    distance_df$sqrt_sd1_2_plus_sd2_2[i] <- NA_real_
  } else {
    distance_df$sqrt_sd1_2_plus_sd2_2[i] <- sqrt(sd1^2 + sd2^2)
  }
  distance_df$weight_new[i] <- distance_df$d2_minus_d1[i] / distance_df$sqrt_sd1_2_plus_sd2_2[i]
}
head(distance_df)
distance_df$weight_new_normalized <- distance_df$weight_new/sum(distance_df$weight_new, na.rm=TRUE)

write.csv(distance_df, "./data/distance_df_new_weight.csv", row.names = FALSE)

# check distributions

distance_df <- read.csv("./data/distance_df_new_weight.csv", header = TRUE, stringsAsFactors = FALSE)
# head(distance_df)
library(ggplot2)
p2=ggplot(distance_df, aes(x=d2_minus_d1)) +
  geom_histogram(binwidth=0.1, fill="green", color="black", alpha=0.7) +
  theme_minimal() +
  labs(title="Distribution of d2 - d1", x="d2 - d1", y="Count")
p1=ggplot(distance_df, aes(x=weight_new)) +
  geom_histogram(binwidth=0.1, fill="blue", color="black", alpha=0.7) +
  theme_minimal() +
  labs(title="Distribution of New Weights", x="New Weight", y="Count")
ggsave("./figure/d2_minus_d1_distribution.png", plot = p2, width = 6, height = 4)
ggsave("./figure/weight_new_distribution.png", plot = p1, width = 6, height = 4)
