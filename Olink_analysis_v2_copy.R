#https://cran.r-project.org/package=OlinkAnalyze
#https://github.com/Olink-Proteomics/OlinkRPackage
#https://cran.r-project.org/web/packages/OlinkAnalyze/vignettes/Vignett.html

#install.packages("OlinkAnalyze")
#install.packages("arrow")
rm(list = ls())

setwd("XX/Analysis2")
# Load OlinkAnalyze
library(OlinkAnalyze)

# Load other libraries used in Vignette
library(dplyr)
library(ggplot2)
library(stringr)
library(arrow)
library(readxl)
library(tidyr)
library(ggfortify)
library(umap)

#######
## Olink Analysis ##

df1 <- read_excel("De-identified.xlsx", sheet = "export")

#NPX: Normalized Protein Expressions

data <- read_NPX("ALL_NPX.parquet")
#View(data)

merged_df <- left_join(data, df1 %>% select(SpecimenBarCode, `Strategy A`), by = c("SampleID" = "SpecimenBarCode"))
filtered_df <- merged_df %>% filter(`Strategy A` %in% c(1, 2)) %>% mutate(`Strategy A` = ifelse(`Strategy A` == 1, "B", ifelse(`Strategy A` == 2, "A", `Strategy A`)))
filtered_df$`Strategy A` <- as.factor(filtered_df$`Strategy A`)
write.csv(filtered_df, "ALL_NPX_Strategy_A_cleaned.csv", row.names = FALSE)
results_ttest <- olink_ttest(df = filtered_df, variable = 'Strategy A')
#View(results_ttest)
geneDataA <- results_ttest %>% select(c(Assay, OlinkID))
write.csv(results_ttest, "ALL_NPX_Strategy_A_results.csv", row.names = FALSE)

df1 <- read_excel("De-identified.xlsx", sheet = "export")
merged_df <- left_join(data, df1 %>% select(SpecimenBarCode, `Strategy B`), by = c("SampleID" = "SpecimenBarCode"))
filtered_df <- merged_df %>% filter(`Strategy B` %in% c(3, 4)) %>% mutate(`Strategy B` = ifelse(`Strategy B` == 4, "B", ifelse(`Strategy B` == 3, "A", `Strategy B`)))
filtered_df$`Strategy B` <- as.factor(filtered_df$`Strategy B`)
write.csv(filtered_df, "ALL_NPX_Strategy_B_cleaned.csv", row.names = FALSE)
results_ttest <- olink_ttest(df = filtered_df, variable = 'Strategy B')
geneDataB <- results_ttest %>% select(c(Assay, OlinkID))
#View(results_ttest)
# Debugging source: https://rdrr.io/cran/OlinkAnalyze/src/R/Olink_ttest.R
#var_levels <- levels(filtered_df[['Strategy B']])
write.csv(results_ttest, "ALL_NPX_Strategy_B_results.csv", row.names = FALSE)


#######
### Deregulated Proteins ###
results_A <- read.csv("ALL_NPX_Strategy_A_results.csv")
results_A_pvalue_cutoff <- results_A %>% select(Assay, estimate, A, B, p.value, Adjusted_pval) %>% filter(p.value<0.05) %>% rename("log2FoldChange" = estimate, "Log2 Average Group 1" = B, "Log2 Average Group 2" = A)
write.csv(results_A_pvalue_cutoff, "results_A_pvalue_cutoff.csv", row.names=FALSE)

results_B <- read.csv("ALL_NPX_Strategy_B_results.csv")
results_B_pvalue_cutoff <- results_B %>% select(Assay, estimate, A, B, p.value, Adjusted_pval) %>% filter(p.value<0.05)  %>% rename("log2FoldChange" = estimate, "Log2 Average Group 3" = A, "Log2 Average Group 4" = B)
write.csv(results_B_pvalue_cutoff, "results_B_pvalue_cutoff.csv", row.names=FALSE)

#######
### Protein Expression Data

# All data cleaned Strategy A
Strategy_A <- read.csv("ALL_NPX_Strategy_A_cleaned.csv")

# Summarize duplicates by taking the mean
Strategy_A <- Strategy_A %>% group_by(SampleID, Assay) %>% summarise(NPX = mean(NPX, na.rm = TRUE), .groups = 'drop')
# Unnormalize by 2^n which gives actual expression values
result_A_temp <- Strategy_A %>% select(Assay, SampleID, NPX) %>% pivot_wider(names_from = Assay, values_from = NPX) %>% mutate(across(-SampleID, ~ 2^(.))) %>% replace(is.na(.), 0) %>% mutate(SampleID = as.numeric(trimws(as.character(SampleID))))
df1 <- read_excel("De-identified.xlsx", sheet = "export")
# Clean data: Remove leading/trailing spaces and ensure SpecimenBarCode is numeric
df1 <- df1 %>% mutate(SpecimenBarCode = as.numeric(trimws(as.character(SpecimenBarCode))), `Strategy A` = as.numeric(trimws(as.character(`Strategy A`)))) %>% filter(`Strategy A` %in% c(1, 2)) %>% select(SpecimenBarCode, `Strategy A`)
result_A_temp <- left_join(result_A_temp, df1, by = c("SampleID" = "SpecimenBarCode"))
#View(result_A_temp)
result_A <- result_A_temp %>% rename(Patient.ID = SampleID, Class = `Strategy A`) %>% select(Patient.ID, Class, everything())
#View(result_A)
write.csv(result_A, "Strategy_A_Protein_Expression.csv", row.names=FALSE)


# Cleaned and deregulated Strategy A
filtered_genelist_A <- results_A_pvalue_cutoff %>% select(Assay)
# Unnormalize by 2^n which gives actual expression values
result_A_temp <- Strategy_A %>% 
  semi_join(filtered_genelist_A, by = "Assay") %>%
  select(Assay, SampleID, NPX) %>%
  pivot_wider(names_from = Assay, values_from = NPX) %>%
  mutate(across(-SampleID, ~ 2^(.))) %>%
  replace(is.na(.), 0) %>%
  mutate(SampleID = as.numeric(trimws(as.character(SampleID))))
# Summarize duplicates by taking the mean
df1 <- read_excel("De-identified.xlsx", sheet = "export")
# Clean data: Remove leading/trailing spaces and ensure SpecimenBarCode is numeric
df1 <- df1 %>% mutate(SpecimenBarCode = as.numeric(trimws(as.character(SpecimenBarCode))), `Strategy A` = as.numeric(trimws(as.character(`Strategy A`)))) %>% filter(`Strategy A` %in% c(1, 2)) %>% select(SpecimenBarCode, `Strategy A`)
result_A_temp <- left_join(result_A_temp, df1, by = c("SampleID" = "SpecimenBarCode"))
#View(result_A_temp)
result_A <- result_A_temp %>% rename(Patient.ID = SampleID, Class = `Strategy A`) %>% select(Patient.ID, Class, everything())
#View(result_A)
write.csv(result_A, "filtered_Strategy_A_Protein_Expression.csv", row.names=FALSE)





# All data cleaned Strategy B
Strategy_B <- read.csv("ALL_NPX_Strategy_B_cleaned.csv")

# Summarize duplicates by taking the mean
Strategy_B <- Strategy_B %>% group_by(SampleID, Assay) %>% summarise(NPX = mean(NPX, na.rm = TRUE), .groups = 'drop')
# Unnormalize by 2^n which gives actual expression values
result_B_temp <- Strategy_B %>% select(Assay, SampleID, NPX) %>% pivot_wider(names_from = Assay, values_from = NPX) %>% mutate(across(-SampleID, ~ 2^(.))) %>% replace(is.na(.), 0) %>% mutate(SampleID = as.numeric(trimws(as.character(SampleID))))
df1 <- read_excel("De-identified.xlsx", sheet = "export")
# Clean data: Remove leading/trailing spaces and ensure SpecimenBarCode is numeric
df1 <- df1 %>% mutate(SpecimenBarCode = as.numeric(trimws(as.character(SpecimenBarCode))), `Strategy B` = as.numeric(trimws(as.character(`Strategy B`)))) %>% filter(`Strategy B` %in% c(3, 4)) %>% select(SpecimenBarCode, `Strategy B`)
result_B_temp <- left_join(result_B_temp, df1, by = c("SampleID" = "SpecimenBarCode"))
#View(result_A_temp)
result_B <- result_B_temp %>% rename(Patient.ID = SampleID, Class = `Strategy B`) %>% select(Patient.ID, Class, everything())
#View(result_A)
write.csv(result_B, "Strategy_B_Protein_Expression.csv", row.names=FALSE)


# Cleaned and deregulated Strategy B
filtered_genelist_B <- results_B_pvalue_cutoff %>% select(Assay)
# Unnormalize by 2^n which gives actual expression values
result_B_temp <- Strategy_B %>% 
  semi_join(filtered_genelist_B, by = "Assay") %>%
  select(Assay, SampleID, NPX) %>%
  pivot_wider(names_from = Assay, values_from = NPX) %>%
  mutate(across(-SampleID, ~ 2^(.))) %>%
  replace(is.na(.), 0) %>%
  mutate(SampleID = as.numeric(trimws(as.character(SampleID))))
# Summarize duplicates by taking the mean
df1 <- read_excel("De-identified.xlsx", sheet = "export")
# Clean data: Remove leading/trailing spaces and ensure SpecimenBarCode is numeric
df1 <- df1 %>% mutate(SpecimenBarCode = as.numeric(trimws(as.character(SpecimenBarCode))), `Strategy B` = as.numeric(trimws(as.character(`Strategy B`)))) %>% filter(`Strategy B` %in% c(3, 4)) %>% select(SpecimenBarCode, `Strategy B`)
result_B_temp <- left_join(result_B_temp, df1, by = c("SampleID" = "SpecimenBarCode"))
#View(result_A_temp)
result_B <- result_B_temp %>% rename(Patient.ID = SampleID, Class = `Strategy B`) %>% select(Patient.ID, Class, everything())
#View(result_A)
write.csv(result_B, "filtered_Strategy_B_Protein_Expression.csv", row.names=FALSE)





#######
### PCA ###

# PCA Entire Protein Expression
result_A <- read.csv("Strategy_A_Protein_Expression.csv")
PCA_data_A <- result_A %>% select(-Class, -Patient.ID)
result_A$`Class` <- as.factor(result_A$`Class`)
PCA_data_A <- PCA_data_A %>% select_if(~ var(.) != 0)
pca_result <- prcomp(PCA_data_A, scale. = TRUE)
pca_plot <- autoplot(pca_result, label=FALSE, data=result_A, colour='Class')
ggsave("PCA_plot_Strategy_A.png", pca_plot)

# PCA deregulated filtered Protein Expression
result_A <- read.csv("filtered_Strategy_A_Protein_Expression.csv")
PCA_data_A <- result_A %>% select(-Class, -Patient.ID)
result_A$`Class` <- as.factor(result_A$`Class`)
PCA_data_A <- PCA_data_A %>% select_if(~ var(.) != 0)
pca_result <- prcomp(PCA_data_A, scale. = TRUE)
pca_plot <- autoplot(pca_result, label=FALSE, data=result_A, colour='Class')
ggsave("filtered_PCA_plot_Strategy_A.png", pca_plot)



# PCA Entire Protein Expression
result_B <- read.csv("Strategy_B_Protein_Expression.csv")
PCA_data_B <- result_B %>% select(-Class, -Patient.ID)
result_B$`Class` <- as.factor(result_B$`Class`)
PCA_data_B <- PCA_data_B %>% select_if(~ var(.) != 0)
pca_result <- prcomp(PCA_data_B, scale. = TRUE)
pca_plot <- autoplot(pca_result, label=FALSE, data=result_B, colour='Class')
ggsave("PCA_plot_Strategy_B.png", pca_plot)

# PCA deregulated filtered Protein Expression
result_B <- read.csv("filtered_Strategy_B_Protein_Expression.csv")
PCA_data_B <- result_B %>% select(-Class, -Patient.ID)
result_B$`Class` <- as.factor(result_B$`Class`)
PCA_data_B <- PCA_data_B %>% select_if(~ var(.) != 0)
pca_result <- prcomp(PCA_data_B, scale. = TRUE)
pca_plot <- autoplot(pca_result, label=FALSE, data=result_B, colour='Class')
ggsave("filtered_PCA_plot_Strategy_B.png", pca_plot)




#######
### UMAP

# Currently PCA is not capturing properly.

result_A <- read.csv("filtered_Strategy_A_Protein_Expression.csv")
umap_data_A <- result_A %>% select(-Class, -Patient.ID)
umap_data_A <- umap_data_A %>% select_if(~ var(.) != 0)
umap_config <- umap.defaults
umap_config$n_neighbors <- 10
umap_config$min_dist <- 0.05
umap_config$metric <- "cosine" 

umap_result <- umap(umap_data_A, config = umap_config)
umap_df <- data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")

umap_df$Class <- result_A$Class

umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Class)) +
  geom_point() +
  theme_minimal() +
  ggtitle("UMAP Projection")

print(umap_plot)
ggsave("filtered_UMAP_Strategy_A.png", umap_plot)




result_A <- read.csv("filtered_Strategy_B_Protein_Expression.csv")
umap_data_A <- result_A %>% select(-Class, -Patient.ID)
umap_data_A <- umap_data_A %>% select_if(~ var(.) != 0)
umap_config <- umap.defaults
umap_config$n_neighbors <- 10
umap_config$min_dist <- 0.05
umap_config$metric <- "cosine" 

umap_result <- umap(umap_data_A, config = umap_config)
umap_df <- data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")

umap_df$Class <- result_A$Class

umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Class)) +
  geom_point() +
  theme_minimal() +
  ggtitle("UMAP Projection")

print(umap_plot)
ggsave("filtered_UMAP_Strategy_B.png", umap_plot)


#######
### Extract immune-realted genes of interest


results_A <- read.csv("U54_ALL_NPX_2024-05-31_Strategy_A_results.csv")
# TAU -> MAPT
# NFL -> NEFL
genes <- c("GFAP", "MAPT", "NEFL", "SNCA", "NPTX2", "IL6", "VEGFA", "TNF", "IL1A", "IL1B", "IL2", "IL4", "IL10", "FOXP3", "GATA3", "TBX21", "STAT1", "STAT3", "STAT6", "FOXP1", "FOXP2", "NFKB1")
results_A_genes_of_interest <- results_A %>% filter(Assay %in% genes) %>% select(Assay, estimate, A, B, p.value, Adjusted_pval) %>% rename("log2FoldChange" = estimate, "Log2 Average Group 1" = B, "Log2 Average Group 2" = A)
print(results_A_genes_of_interest)
write.csv(results_A_genes_of_interest, "Strategy_A_immune_related_genes_log2FC_pvalue.csv", row.names = FALSE)

results_B <- read.csv("U54_ALL_NPX_2024-05-31_Strategy_B_results.csv")
genes <- c("GFAP", "MAPT", "NEFL", "SNCA", "NPTX2", "IL6", "VEGFA", "TNF", "IL1A", "IL1B", "IL2", "IL4", "IL10", "FOXP3", "GATA3", "TBX21", "STAT1", "STAT3", "STAT6", "FOXP1", "FOXP2", "NFKB1")
results_B_genes_of_interest <- results_B %>% filter(Assay %in% genes) %>% select(Assay, estimate, A, B, p.value, Adjusted_pval)  %>% rename("log2FoldChange" = estimate, "Log2 Average Group 3" = A, "Log2 Average Group 4" = B)
print(results_B_genes_of_interest)
write.csv(results_B_genes_of_interest, "Strategy_B_immune_related_genes_log2FC_pvalue.csv", row.names = FALSE)



#############################################################################
#############################################################################