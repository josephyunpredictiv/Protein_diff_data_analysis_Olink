# Protein_diff_data_analysis_Olink
Protein differential expression data analysis from Olink. Bascially we're using olink software... 

The file also contains code for PCA and Umap. This is created as a memo to self. However, hopefully this helps others trying Olink software. 

#Important note, make sure you carefuly check Log2FoldChange!!!!!! 
In the Olink source code: https://rdrr.io/cran/OlinkAnalyze/src/R/Olink_ttest.R 

The line where it says:
var_levels <- levels(df[[variable]])

Orders the numbers numerically or letteres alphabetically which means if your groups are E for experimental and C for control, the Log2FoldChange is actually backwards. You can account for this by manually changing around control and experiemental in the code or editing the results afterwards. Like set C to A and E to B so the t_test orders correctly. All other expression values should be the same no matter what order they are put into. And yes, as of June 15th, the Olink still does not have a way to specify t_test control group. Please add it please!!!!
