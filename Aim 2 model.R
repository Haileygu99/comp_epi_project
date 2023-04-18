# Pre-processing code: 
rm(list=ls()) 
library(tidyverse) 
library(openxlsx) 
library(omics)
library(cluster)
library(gower)
library(factoextra)
library(Rtsne)
library(clustMixType)
library(glmnet)
library(nnet)
library(grpreg)
library(clusterCons)
library(clustMixType)
library(glue)
library(caret)
#Loading Data
df <- read.csv('./imputed_baseline.csv') 
df$case<- ifelse(df$cohort.x %in% c("cohort_a", "cohort_b","cohort_c"), "1", "0")
df<-subset(df,case==1)
table(df$case)

#Data pre-processing-----------------------
cols_subtypes<-c('eosinophil_sputum_percentage','neutrophil_sputum_percentage','eosinophil_serum_count',
                 'neutrophil_serum_count','nasal_polyps','age_of_onset',
                 'atopy')
df_subtypes<-df[,c(cols_subtypes)]
categorical<-c('nasal_polyps', 'atopy')
df_subtypes <- df_subtypes %>% mutate(across(categorical, as.factor))

#Re-scaling-----------------------
# Re-scaling
# Standardize only continuous columns
continuous_columns <- setdiff(colnames(df_subtypes), categorical)
df_subtypes[, continuous_columns] <- scale(df_subtypes[, continuous_columns])

#Hierarchical Clustering--------------------------
gower_dist <- daisy(df_subtypes, metric = "gower")
hc <- hclust(gower_dist, method = "ward.D2")
fviz_nbclust(df_subtypes, FUNcluster = hcut, method = "wss", diss = gower_dist) + theme_minimal() # elbow plot
clusters_hier <- factor(cutree(hc, k = 3))
plot(hc, cex = 0.6, hang = -1)
rect.hclust(hc, k = 3, border = "red")

#K-prototypes --------------------------
 
# Set the number of repetitions for each k value
n_reps <- 50
# Initialize a vector to store average silhouette scores
avg_silhouette_scores <- numeric()
# Perform k-prototypes clustering for different k values
for (k in 2:7) {
  # Initialize a vector to store silhouette scores for current k
  silhouette_scores <- numeric()
  for (rep in 1:n_reps) {
    # Perform k-prototypes clustering with a random seed
    set.seed(rep)
    kproto_result <- kproto(df_subtypes, k, verbose=FALSE)
    # Calculate silhouette scores
    silhouette_result <- silhouette(kproto_result$cluster, gower_dist)
    avg_silhouette <- mean(silhouette_result[, "sil_width"])
    # Append silhouette score to the vector
    silhouette_scores <- append(silhouette_scores, avg_silhouette)
  }
  # Calculate the average silhouette score for the current k and append it to the vector
  avg_silhouette_scores <- append(avg_silhouette_scores, mean(silhouette_scores))
}
# Plot average silhouette scores
plot(2:7, avg_silhouette_scores, type = "b", xlab = "Number of clusters (k)", ylab = "Average silhouette score", main = glue("Average silhouette scores over {n_reps} iterations"))

#Use best K
set.seed(123)
k <- 2 # Number of clusters
kproto_result <- kproto(df_subtypes, k, verbose = FALSE)
clusters_kproto <- factor(kproto_result$cluster)

#T_SNE -------------------
# Loop through perplexities from 10 to 90, in intervals of 10
for (perplexity in seq(10, 90, by = 10)) {
  # Apply the t-SNE algorithm on the Gower distance matrix with the current perplexity
  set.seed(123) # Set a seed for reproducibility
  tsne_result <- Rtsne(df_subtypes, distance_matrix = TRUE, dims = 2, perplexity = perplexity)
  
  # Create a data frame with t-SNE results
  tsne_data <- as.data.frame(tsne_result$Y)
  colnames(tsne_data) <- c("TSNE1", "TSNE2")
  
  # Add cluster labels or categories if you have them
  tsne_data$label <- clusters_kproto
  
  # Plot the t-SNE result
  plot_title <- paste0("t-SNE Visualization with Gower Distance (Perplexity = ", perplexity, ")")
  p <- ggplot(tsne_data, aes(x = TSNE1, y = TSNE2, color = label)) +
    geom_point() +
    theme_minimal() +
    labs(title = plot_title)
  
  print(p) # Print the plot in the loop
}

#Describe clusters ---------------
df_subtypes$clusters<-clusters_kproto
df$clusters<-clusters_kproto
cluster_summary <- df_subtypes %>%
  group_by(clusters) %>%
  summarise(across(.cols = where(is.numeric), .fns = list(mean = mean), .names = "{col}_{.fn}"), na.rm = TRUE)
print(cluster_summary)

# Denoising---------------------------
proteins <- c("CCL18", "IL18", "CCL5", "i.TAC", "TARC", "IL.10", "IL.16", "MIP.1a", "PAPP.A", "IL17", "TSLP", "SP.A", "CTACK", "KL.6", "MPO", "EDN", "IL.6", "TNF.A", "MIG", "CCL20", "IP.10", "IL.4", "IL.5")
X <- df[, proteins]
model = mlmer(X ~ (1|PLATE.ID), data=df, save.residuals = TRUE, save.ranks = TRUE)
saveRDS(model$residuals, "proteins_denoised.rds")
X <- readRDS("proteins_denoised.rds")
# Covariates
covars <- c("Education", "Occupation", "Sex", "Age", "Race", "Body.Mass.Index..kg.m2.")
X_covars <- df[, covars]
X_covars$Education<-as.factor(X_covars$Education)
X_covars$Occupation<-as.factor(X_covars$Occupation)
X_covars$Sex<-as.factor(X_covars$Sex)
X_covars$Race<-as.factor(X_covars$Race)
# Create dummy variables (one-hot encoding) for categorical covariates
dummies <- dummyVars(~ ., data = X_covars, fullRank = TRUE)
X_covars_onehot <- data.frame(predict(dummies, newdata = X_covars))
# Combine the denoised proteins (X) and the one-hot encoded covariates (X_covars_onehot) into a single data frame
X_combined <- cbind(X, X_covars_onehot)
# Create a penalty factor vector
penalty_factors <- rep(1, ncol(X_combined))
# Set the penalty factors for the covariates to 0
penalty_factors[(ncol(X) + 1):ncol(X_combined)] <- 0
# Group Lasso-----------------------------
pr_groups <- c()
pr_groups[c("CCL18", "CCL5", "i.TAC", "TARC", "MIP.1a", "MIG", "CCL20", "IP.10", "CTACK")] <- 1
pr_groups[c("IL.10", "IL.16", "IL17", "TSLP", "IL.4", "IL.5")] <- 2
pr_groups[c("PAPP.A", "KL.6", "MPO", "EDN")] <- 3
pr_groups[c("IL18", "SP.A", "IL.6", "TNF.A")] <- 4
pr_groups <- pr_groups[match(proteins, names(pr_groups))]
pr_groups_covars <- rep("0", length(X_covars_onehot)) # Assign group 5 to all covariates
pr_groups <- c(pr_groups, pr_groups_covars)
# Perform group Lasso with cross-validation
cv_fit <- cv.grpreg(X_combined, clusters_kproto, groups=pr_groups, penalty="grLasso", family="binomial", penalty.factor=penalty_factors)
plot(cv_fit)
# Extract the coefficients at lambda.1se
coefficients_lambda_1se <- coef(cv_fit, s = cv_fit$lambda.1se)
print(coefficients_lambda_1se)
# Remove the intercept
coefficients_lambda_1se <- coefficients_lambda_1se[-1]
# Create a data frame with feature names, group names, and coefficients
features <- colnames(X_combined)
coef_df <- data.frame(
  feature = features,
  group = pr_groups,
  coefficient = coefficients_lambda_1se
)
# Order the data frame by group and feature
coef_df <- coef_df[order(coef_df$group, coef_df$feature), ]
# Plot the coefficients in a bar plot
barplot(
  coef_df$coefficient,
  names.arg = coef_df$feature,
  xlab = "Features",
  ylab = "Coefficients",
  main = "Coefficients for each feature and group at lambda.1se",
  col = coef_df$group,
  las = 2,
  cex.names = 0.5
)
# Add a legend
#legend("right", legend = unique(coef_df$group), fill = unique(coef_df$group), title = "Group")
