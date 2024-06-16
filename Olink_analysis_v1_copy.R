#https://cran.r-project.org/package=OlinkAnalyze
#https://github.com/Olink-Proteomics/OlinkRPackage
#https://cran.r-project.org/web/packages/OlinkAnalyze/vignettes/Vignett.html

#install.packages("OlinkAnalyze")
#install.packages("arrow")

setwd("XXX/Data_2024_06_03")
# Load OlinkAnalyze
library(OlinkAnalyze)

# Load other libraries used in Vignette
library(dplyr)
library(ggplot2)
library(stringr)
library(arrow)
library(readxl)
library(tidyr)
library("ggfortify")

#######
## Olink Analysis ##

df1 <- read_excel("De-identified4.xlsx", sheet = "export")

#NPX: Normalized Protein Expressions

#Step1: Preprocessing (NPX=normalized protein expression)
data <- read_NPX("ALL_NPX.parquet")
#View(data)

merged_df <- left_join(data, df1 %>% select(SpecimenBarCode, `Strategy A`), by = c("SampleID" = "SpecimenBarCode"))
filtered_df <- merged_df %>% filter(`Strategy A` %in% c(1, 2))
filtered_df$`Strategy A` <- as.factor(filtered_df$`Strategy A`)
write.csv(filtered_df, "ALL_NPX_Strategy_A_cleaned.csv", row.names = FALSE)
#Step2: Differential protein expression testing
results_ttest <- olink_ttest(df = filtered_df, variable = 'Strategy A')
#View(results_ttest)

write.csv(results_ttest, "ALL_NPX_Strategy_A_results.csv", row.names = FALSE)


merged_df <- left_join(data, df1 %>% select(SpecimenBarCode, `Strategy B`), by = c("SampleID" = "SpecimenBarCode"))
filtered_df <- merged_df %>% filter(`Strategy B` %in% c(3, 4))
filtered_df$`Strategy B` <- as.factor(filtered_df$`Strategy B`)
write.csv(filtered_df, "ALL_NPX_Strategy_B_cleaned.csv", row.names = FALSE)
#Step2: Differential protein expression testing
results_ttest <- olink_ttest(df = filtered_df, variable = 'Strategy B')
#View(results_ttest)

write.csv(results_ttest, "ALL_NPX_Strategy_B_results.csv", row.names = FALSE)



#######
### PCA ###

Strategy_A <- read.csv("ALL_NPX_Strategy_A_cleaned.csv")

# Unnormalize by 2^n which gives actual expression values
result_A_temp <- Strategy_A %>%
  select(OlinkID, SampleID, NPX) %>%
  pivot_wider(names_from = OlinkID, values_from = NPX) %>%
  mutate(across(-SampleID, ~ 2^(.))) %>%
  replace(is.na(.), 0) %>%
  mutate(SampleID = as.numeric(trimws(as.character(SampleID))))

# Load the Excel file
df1 <- read_excel("De-identified4.xlsx", sheet = "export")

# Clean data: Remove leading/trailing spaces and ensure SpecimenBarCode is numeric
df1 <- df1 %>%
  mutate(SpecimenBarCode = as.numeric(trimws(as.character(SpecimenBarCode))),
         `Strategy A` = as.numeric(trimws(as.character(`Strategy A`)))) %>%
  filter(`Strategy A` %in% c(1, 2)) %>% select(SpecimenBarCode, `Strategy A`)

# Debugging: Check for matching values
matching_values <- intersect(result_A_temp$SampleID, df1$SpecimenBarCode)
print(matching_values)

# Perform the left join
result_A_temp <- left_join(result_A_temp, df1, by = c("SampleID" = "SpecimenBarCode"))
View(result_A_temp)

result_A <- result_A_temp %>%
  rename(Patient.ID = SampleID, Class = `Strategy A`) %>%
  select(Patient.ID, Class, everything())

View(result_A)

write.csv(result_A, "Strategy_A_Protein_Expression.csv", row.names=FALSE)


#result_A <- read.csv("Strategy_A_Protein_Expression.csv")
PCA_data_A <- result_A %>% select(-Class, -Patient.ID)
result_A$`Class` <- as.factor(result_A$`Class`)
PCA_data_A <- PCA_data_A %>% select_if(~ var(.) != 0)
pca_result <- prcomp(PCA_data_A, scale. = TRUE)
autoplot(pca_result, label=FALSE, data=result_A, colour='Class')





Strategy_B <- read.csv("ALL_NPX_Strategy_B_cleaned.csv")

# Unnormalize by 2^n which gives actual expression values
result_B_temp <- Strategy_B %>%
  select(OlinkID, SampleID, NPX) %>%
  pivot_wider(names_from = OlinkID, values_from = NPX) %>%
  mutate(across(-SampleID, ~ 2^(.))) %>%
  replace(is.na(.), 0) %>%
  mutate(SampleID = as.numeric(trimws(as.character(SampleID))))

# Load the Excel file
df1 <- read_excel("De-identified4.xlsx", sheet = "export")

# Clean data: Remove leading/trailing spaces and ensure SpecimenBarCode is numeric
df1 <- df1 %>%
  mutate(SpecimenBarCode = as.numeric(trimws(as.character(SpecimenBarCode))),
         `Strategy B` = as.numeric(trimws(as.character(`Strategy B`)))) %>%
  filter(`Strategy B` %in% c(3, 4)) %>% select(SpecimenBarCode, `Strategy B`)

# Debugging: Check for matching values
matching_values <- intersect(result_B_temp$SampleID, df1$SpecimenBarCode)
print(matching_values)

# Perform the left join
result_B_temp <- left_join(result_B_temp, df1, by = c("SampleID" = "SpecimenBarCode"))
View(result_B_temp)

result_B <- result_B_temp %>%
  rename(Patient.ID = SampleID, Class = `Strategy B`) %>%
  select(Patient.ID, Class, everything())

View(result_B)

write.csv(result_B, "Strategy_B_Protein_Expression.csv", row.names=FALSE)

PCA_data_B <- result_B %>% select(-Class, -Patient.ID)
result_B$`Class` <- as.factor(result_B$`Class`)
PCA_data_B <- PCA_data_B %>% select_if(~ var(.) != 0)
pca_result <- prcomp(PCA_data_B, scale. = TRUE)
autoplot(pca_result, label=FALSE, data=result_B, colour='Class')



#######
### Deregulated Proteins ###
results_A <- read.csv("ALL_NPX_Strategy_A_results.csv")
results_A_pvalue_cutoff <- results_A %>% select(Assay, estimate, X1, X2, p.value, Adjusted_pval) %>% filter(p.value<0.05) %>% rename("log2FoldChange" = estimate, "Log2 Average Group 1" = X1, "Log2 Average Group 2" = X2)
write.csv(results_A_pvalue_cutoff, "results_A_pvalue_cutoff.csv", row.names=FALSE)

results_B <- read.csv("ALL_NPX_Strategy_B_results.csv")
results_B_pvalue_cutoff <- results_B %>% select(Assay, estimate, X3, X4, p.value, Adjusted_pval) %>% filter(p.value<0.05)  %>% rename("log2FoldChange" = estimate, "Log2 Average Group 3" = X3, "Log2 Average Group 4" = X4)
write.csv(results_B_pvalue_cutoff, "results_B_pvalue_cutoff.csv", row.names=FALSE)

#results_wilcox <- olink_wilcox(df = npx_data1, variable = 'Treatment')
#View(results_wilcox)

#############################################################################
#############################################################################