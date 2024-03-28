library(tidyr)
library(dplyr)
library(pheatmap)

data <- read.csv("merge.diff.txt",  sep="\t",check.names = FALSE, header = TRUE)

#caculate FC, all0 are delete :3808
data <- data %>%
  mutate(
    meanCD_over_meanHC = mean_CD / mean_HC,
    medianCD_over_medianHC = median_CD / median_HC,
    meanCD_minus_meanHC = mean_CD - mean_HC,
    medianCD_minus_medianHC = median_CD - median_HC
  )

data <- data %>%
  mutate(
    FC = case_when(
      medianCD_minus_medianHC > 0 ~ 1,
      medianCD_minus_medianHC < 0 ~ -1,
      medianCD_minus_medianHC == 0 & meanCD_minus_meanHC > 0 ~ 1,
      medianCD_minus_medianHC == 0 & meanCD_minus_meanHC < 0 ~ -1,
      TRUE ~ NA_real_  # Assign NA for conditions that do not match any of the above
    )
  )

data <- data %>%filter(!is.na(FC))

#sig score
data <- data %>%
  mutate(significant_count = rowSums(select(., linear_mixed_all_cofactor_p.x, linear_mixed_all_nocofactor_p.y, pvalue_lm, pvalue_kw) < 0.05, na.rm = TRUE))  %>%
  mutate(min_0.05 = as.integer(pmin(linear_mixed_all_cofactor_p.x, linear_mixed_all_nocofactor_p.y, pvalue_lm, pvalue_kw, na.rm = TRUE) < 0.05))

#rank score--队列内rank-desending, 未有p值的不排序，四种方法相加/方法数量
data <- data %>%
  group_by(cohort) %>%
  mutate(
    rank_linear_mixed_all_cofactor_p.x = rank(desc(linear_mixed_all_cofactor_p.x), ties.method = "average", na.last = "keep"),
    rank_linear_mixed_all_nocofactor_p.y = rank(desc(linear_mixed_all_nocofactor_p.y), ties.method = "average", na.last = "keep"),
    rank_pvalue_lm = rank(desc(pvalue_lm), ties.method = "average", na.last = "keep"),
    rank_pvalue_kw = rank(desc(pvalue_kw), ties.method = "average", na.last = "keep")
  ) %>%
  ungroup()

data <- data %>%
  rowwise() %>% 
  mutate(
    average_rank = (
      coalesce(rank_linear_mixed_all_cofactor_p.x, 0) + 
        coalesce(rank_linear_mixed_all_nocofactor_p.y, 0) + 
        coalesce(rank_pvalue_lm, 0) + 
        coalesce(rank_pvalue_kw, 0)
    ) / sum(!is.na(c(rank_linear_mixed_all_cofactor_p.x, 
                     rank_linear_mixed_all_nocofactor_p.y, 
                     rank_pvalue_lm, 
                     rank_pvalue_kw)))
  ) %>%
  ungroup()

data <- data %>%
  mutate(
    rank_FC = average_rank*FC,
    sig_FC = significant_count*FC,
    min_FC = min_0.05*FC,
  )

write.csv(data,"spceis.score.csv")


#根据3个关键Score长转宽rank_FC，sig_FC，min_FC
library(tidyr)

#替换宽变量
data_wide_rank_FC <- data %>%
  pivot_wider(
    names_from = cohort,
    values_from = rank_FC,  #rank_FC，sig_FC，min_FC
    names_prefix = "rank_FC_", # rank_FC，sig_FC，min_FC
    id_cols = gene_name2.x
  )

#only for rank_FC, normalize column to each column sums to 1
#nor-sum1 (not correlation)
#data_wide_rank_FC_normalized <- data_wide_rank_FC %>%
  mutate(across(contains("rank"), ~ .x / sum(.x, na.rm = TRUE)))
#nor-min-max (no-positive and negative)
#data_wide_rank_FC_normalized <- data_wide_rank_FC %>%
  mutate(across(contains("rank"), ~ (.x - min(.x, na.rm = TRUE)) / (max(.x, na.rm = TRUE) - min(.x, na.rm = TRUE))))
# Percentile Rank 
#data_wide_rank_FC_normalized <- data_wide_rank_FC %>%
  mutate(across(contains("rank"), function(x) {
    ecdf_x <- ecdf(x)  
    ecdf_x(x)  
  }))
# Z-Score Normalization(choose)
data_wide_rank_FC_normalized <- data_wide_rank_FC %>%
  mutate(across(contains("rank"), ~ ( .x - mean(.x, na.rm = TRUE)) / sd(.x, na.rm = TRUE)))

data_wide_rank_FC_normalized <- data_wide_rank_FC_normalized %>%
  rename_with(~paste0(., "_normalized"), -gene_name2.x)
merged_data <- merge(data_wide_rank_FC, data_wide_rank_FC_normalized, by = "gene_name2.x")
numeric_columns <- select(merged_data, -gene_name2.x)
cor_matrix <- cor(numeric_columns, use = "complete.obs")
print(cor_matrix)

#查看队列相关性，以便排除差的队列
numeric_data_wide_rank_FC <- data_wide_rank_FC %>% select(-gene_name2.x)
cor_matrix <- cor(numeric_data_wide_rank_FC, use = "complete.obs", method = "pearson")
print(cor_matrix)
png("normal-zscore_rank_cor.png", width = 800, height = 800)
pheatmap(cor_matrix,
         display_numbers = TRUE, # Optional: Display correlation coefficients
         cluster_cols=FALSE,
         cluster_rows=FALSE,
         #clustering_distance_rows = "euclidean",
         #clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("blue", "white", "red"))(50), # Custom color scale
         title = "Correlation Matrix Heatmap")
dev.off()

p_value_matrix <- matrix(nrow = ncol(numeric_data_wide_rank_FC ), ncol = ncol(numeric_data_wide_rank_FC ), dimnames = list(colnames(numeric_data_wide_rank_FC ), colnames(numeric_data_wide_rank_FC )))

# Loop through each pair of columns and perform cor.test
for (i in 1:(ncol(numeric_data_wide_rank_FC )-1)) {
  for (j in (i+1):ncol(numeric_data_wide_rank_FC )) {
    test_result <- cor.test(numeric_data_wide_rank_FC [[i]], numeric_data_wide_rank_FC [[j]], method = "pearson")
    
    # Store the p-value in the matrix
    p_value_matrix[i, j] <- test_result$p.value
    p_value_matrix[j, i] <- test_result$p.value # Matrix is symmetric
  }
}

write.csv(p_value_matrix,"p_value_matrix.csv")
write.csv(cor_matrix,"cor_matrix.csv")
write.csv(data_wide_rank_FC,"data_wide_min_FC.csv")  #name
write.csv(data_wide_rank_FC_normalized,"data_wide_rank_FC_normalized.csv")

          