#Read the significant k-mer data filtered out in the previous steps and put it into a model. Use linear mixed effects to filter
#After filtering, the data is nested, i.e., K+1 contains k
library(readr)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(MASS)
library(caret)
merged_data_2mer <- read.csv("merged_data_2mer.csv")
merged_data_2mer <- data.frame(merged_data_2mer)
X_pca_2mer <- read.csv("/X_pca_2mer.csv")
X_pca_2mer <- data.frame(X_pca_2mer)
merged_pca_2mer <- cbind(merged_data_2mer[,c(1,2,3)],X_pca_2mer)
write.csv(merged_pca_2mer,"merged_pca_2mer.csv",row.names = FALSE)

merged_data_3mer <- read.csv("merged_data_3mer.csv")
merged_data_3mer <- data.frame(merged_data_3mer)
X_pca_3mer <- read.csv("X_pca_3mer.csv")
X_pca_3mer <- data.frame(X_pca_3mer)
merged_pca_3mer <- cbind(merged_data_3mer[,c(1,2,3)],X_pca_3mer)
write.csv(merged_pca_3mer,"merged_pca_3mer.csv",row.names = FALSE)


merged_data_4mer <- read.csv("merged_data_4mer.csv")
merged_data_4mer <- data.frame(merged_data_4mer)
X_pca_4mer <- read.csv("X_pca_significant_4mer.csv")
X_pca_4mer <- data.frame(X_pca_4mer)
merged_pca_4mer <- cbind(merged_data_4mer[,c(1,2,3)],X_pca_4mer)
write.csv(merged_pca_4mer,"merged_pca_4mer.csv",row.names = FALSE)

merged_data_5mer <- read.csv("merged_data_5mer.csv")
merged_data_5mer <- data.frame(merged_data_5mer)
X_pca_5mer <- read.csv("X_pca_significant_5mer.csv")
X_pca_5mer <- data.frame(X_pca_5mer)
merged_pca_5mer <- cbind(merged_data_5mer[,c(1,2,3)],X_pca_5mer)
write.csv(merged_pca_5mer,"merged_pca_5mer.csv",row.names = FALSE)

# Aggregate all data together
merged_pca_total01 <- merge(merged_pca_2mer,merged_pca_3mer,by = c("Gene","SampleID","Expression","Initial_value"))
merged_pca_total02 <- merge(merged_pca_total01,merged_pca_4mer,by = c("Gene","SampleID","Expression","Initial_value"))
merged_pca_total03 <- merge(merged_pca_total02,merged_pca_5mer,by = c("Gene","SampleID","Expression","Initial_value"))
write.csv(merged_pca_total03,"merged_pca_total03.csv",row.names = FALSE)

#Use linear mixed effects to filter nested data, select the final significant k-mer data, and perform modelling
response_kmer <- merged_pca_total03$Expression
predictors_kmer <- merged_pca_total03[, c(-1,-2,-3)]  # Remove gene ranking, etc.

anova_results_kmer <- list()
i <- 1
for (kmer in colnames(predictors_2mer)) {
  print(paste("Fitting model for k-mer:", kmer, " - ", i))  
  
  formula <- as.formula(paste("Expression ~ ", kmer, "+ Initial_value + (1 | Gene) + (1 | SampleID)"))
  model <- lmer(formula, data = merged_pca_total03)
  anova_results_kmer[[kmer]] <- model
  
  i <- i + 1
}

# Define a function to determine the analysis results and output the conclusion
analyze_results <- function(anova_results, p_value_threshold = 0.05, output_prefix) {
  # Initialise the result list and the list of significant 2-mer names
  results_list <- list()
  significant_k_mer <- data.frame(kmer = character())
  
  # Results of traversing each k-mer
  for (kmer in names(anova_results)) {
    model_summary <- summary(anova_results[[kmer]])
    coef_summary <- coef(model_summary)
    
    if (nrow(coef_summary) > 1) {
      p_value <- coef_summary[2, "Pr(>|t|)"]
      t_value <- coef_summary[2, "t value"]
      
     
      print(paste("Checking k-mer:", kmer))
      print(paste("P-value:", p_value))
      print(paste("T-value:", t_value))
      
      if (!is.na(p_value) && p_value < p_value_threshold) {
        # Extract significant results
        estimate <- coef_summary[2, "Estimate"]
        std_error <- coef_summary[2, "Std. Error"]
        
        # Build result data frame
        result_text <- paste("The k-mer", kmer, "has an estimate of", estimate, 
                             "with a standard error of", std_error, 
                             "a t value of", t_value, 
                             "and a p value of", p_value, 
                             ". This suggests that the effect of", kmer, 
                             "on gene expression is significant.")
        
        results_list[[kmer]] <- data.frame(
          kmer = kmer,
          Estimate = estimate,
          StdError = std_error,
          TValue = t_value,
          PValue = p_value,
          ResultText = result_text
        )
        
        # Store the significant 2-mer names in the data frame
        significant_k_mer <- rbind(significant_k_mer, data.frame(kmer = kmer))
      }
    }
  }
  
  # Convert the result list to a data frame
  results_df <- bind_rows(results_list)
  
  # Multiple test correction
  if (nrow(results_df) > 0) {
    results_df$AdjustedPValue <- p.adjust(results_df$PValue, method = "BH")
    significant_results <- results_df[results_df$AdjustedPValue < 0.05, ]
    
    # Build result text
    significant_results$ResultText <- paste("The k-mer", significant_results$kmer, 
                                            "has an estimate of", significant_results$Estimate, 
                                            "with a standard error of", significant_results$StdError, 
                                            "a t value of", significant_results$TValue, 
                                            "and a p value of", significant_results$PValue, 
                                            ". This suggests that the effect of", significant_results$kmer, 
                                            "on gene expression is significant after adjustment.")
    
    write.csv(significant_results, paste0(output_prefix, "_motifs.csv"), row.names = FALSE)
    write.csv(significant_k_mer, paste0(output_prefix, ".csv"), row.names = FALSE)
    return(significant_results)
  } else {
    return("No significant motifs found.")
  }
}



significant_motifs_df_kmer <- analyze_results(anova_results_kmer,,"significant_k_mer")







