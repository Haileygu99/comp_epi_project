library(lme4)
library(sgPLS)
library(utils)
library(stringr)
library(pheatmap)
library(mixOmics)
library(pls)
library(RColorBrewer)
source("/Users/wenjiazhang/Documents/MSc_HDA/Computational Epi/4/04-practical/BCL_analyses/Scripts/pls_functions.R")
MyPal = brewer.pal("Paired", n = 12)

### loading data
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd("/Users/wenjiazhang/Desktop")
data <- read.csv("imputed_baseline.csv")

### making a case/control col
data$case<- ifelse(data $cohort.x %in% c("cohort_a", "cohort_b","cohort_c"), "1", "0")
table(data$case)

# re-code categorical variables back 
data$Education <- ifelse(data$Education == 1,"Secondary_or_lower",
                                 ifelse(data$Education ==2,"University", 
                                        ifelse(data$Education ==3,"Higher",NA)))
data$Occupation <- ifelse(data$Occupation == 1,"employed",
                                  ifelse(data$Occupation ==2, "housekeeping",
                                         ifelse(data$Occupation ==3, "other", 
                                                ifelse(data$Occupation == 4, "retired", 
                                                       ifelse(data$Occupation == 5, "student",
                                                              ifelse(data$Occupation == 6, "unemployed",NA))))))

data$Ethnic_father <- ifelse(data$Ethnic_father == 1,"white_caucasian",
                                     ifelse(data$Ethnic_father == 2, "asian",
                                            ifelse(data$Ethnic_father == 3, "black_african",
                                                   ifelse(data$Ethnic_father == 4, "other",NA))))

data$Ethnic_mother <- ifelse(data$Ethnic_mother == 1, "white_caucasian",
                                     ifelse(data$Ethnic_mother == 2, "asian",
                                            ifelse(data$Ethnic_mother == 3, "black_african", 
                                                   ifelse(data$Ethnic_mother == 4, "other",NA))))

data$Race <- ifelse(data$Race == 1,"white_caucasian",
                            ifelse(data$Race ==2, "asian",
                                   ifelse(data$Race == 3, "black_african",
                                          ifelse(data$Race == 4,"other",NA))))
data$ocs <- ceiling(data$ocs)
data$ocs <- ifelse(data$ocs == 3, "at_least_once_per_day",
                           ifelse(data$ocs == 2, "less_than_once_per_day",
                                  ifelse(data$ocs == 0, "never",
                                         ifelse(data$ocs ==1, "Previous", NA))))

names(data)[names(data) == "Body.Mass.Index..kg.m2."] <- "BMI"

# Preparing data - linear mixed models 
denoised = NULL
Beta_pooled = NULL
pvalue_pooled = NULL

Groups=c("CCL18","CCL5","i.TAC","TARC","MIP.1a","MIG","CCL20","IP.10",
         "IL.10","IL.16","IL17","TSLP","IL.4","IL.5",
         "PAPP.A","KL.6","MPO","EDN",
         "IL18","SP.A","IL.6","TNF.A","CTACK")
df_protein = data[,Groups]

f0='df_protein[,k] ~ Age + Sex + Race + BMI + Smoking.Status + Occupation + Education + icu + (1 | PLATE.ID)'
f1=paste(f0, '+ case')
 
for (k in seq(1:ncol(df_protein))){
  print(k)
  model=lmer(as.formula(f1), data=data, REML=FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-04)))
  model0=lmer(as.formula(f0), data=data, REML=FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-04)))
  pvalue_pooled=c(pvalue_pooled, anova(model, model0)$'Pr(>Chisq)'[2])
  beta=fixef(model)[c('(Intercept)', 'case1')]
  Beta_pooled=c(Beta_pooled, fixef(model)['case1'])
  X=cbind(rep(1, length(data$case)), as.numeric(data$case))
  denoised=cbind(denoised, (X%*%beta + resid(model)))
  }
colnames(denoised)=colnames(df_protein)
rownames(denoised)=rownames(df_protein)
saveRDS(denoised, "Proteins_denoised.rds")

Table_pooled=cbind(Beta_pooled, pvalue_pooled)
rownames(Table_pooled)=colnames(df_protein)

# Preparing data 
Groups=c("CCL18","CCL5","i.TAC","TARC","MIP.1a","MIG","CCL20","IP.10",
         "IL.10","IL.16","IL17","TSLP","IL.4","IL.5",
         "PAPP.A","KL.6","MPO","EDN",
         "IL18","SP.A","IL.6","TNF.A","CTACK")
X_pooled = data[,Groups]

# non-penalised PLS-DA model
MyPLSDA_pooled <- plsda(X_pooled,data$case,ncomp = 1)
MyPLSDA_pooled$loadings$X
MyPLSDA_pooled$loadings$Y
MyPLSDA_pooled$explained_variance

# running an uncalibrated sparse PLS-DA model
MysPLSDA_pooled <- splsda(X_pooled,data$case, ncomp = 1,
                          keepX = 5)
MysPLSDA_pooled$loadings$X
MysPLSDA_pooled$loadings$X[MysPLSDA_pooled$loadings$X != 0, ]

# setting group membership 
Xgroups = c(8,14,18)

MygPLSDA_pooled <- gPLSda(X_pooled, data$case, ncomp = 1,
                          ind.block.x = Xgroups, keepX = 1)

MygPLSDA_pooled$loadings$X

print(MygPLSDA_pooled)

# we are adding a second layer of sparsity by selecting the variables within the group
MysgPLSDA_pooled <- sgPLSda(X_pooled,data$case, ncomp = 1,
                            ind.block.x = Xgroups, keepX = 1, alpha.x = 0.1)
MysgPLSDA_pooled$loadings$X

#calibration - opti num of variables = 6
set.seed(1)
res_splsda = CalibratesPLSDA(dataX = X_pooled, dataY = data$case,
                             ncomp = 1, Nrepeat = 5)
PlotCalib(res = res_splsda)

# sppls calibration - 2 groups with a sparsity parameter of 0.6 leading to smallest MSEP
set.seed(1)
res_sgplsda = CalibratesgPLSDA(dataX = X_pooled, dataY = data$case,
                               ncomp = 1, Nrepeat = 5, Xgroups = Xgroups)
pdf("my_sgPLSDA_caliplot.pdf")
PlotCalib(res = res_sgplsda, type = "sgPLSDA")
dev.off()

# running the pls-da -
Loadings = cbind(MysPLSDA_pooled$loadings$X, MysgPLSDA_pooled$loadings$X,
                 rep(NA, 23), rep(NA, 23))
Loadings = as.vector(t(Loadings))
Loadings = Loadings[-c(length(Loadings) - 1, length(Loadings))]
par(mar = c(10, 5, 3, 3))
pdf("pls_loading.pdf")
plot(Loadings, col = c(MyPal[6], MyPal[10], NA, NA),
     xaxt = "n", ylab = "Loadings Coefficients", type = "h",
     lwd = 3, xlab = "")
axis(1, at = seq(1.5, 23 * 4, by = 4), labels = colnames(MyPLSDA_pooled$X),
     las = 2)
axis(1, at = c(0, Xgroups, 23) * 4, line = 6, labels = NA)
axis(1, at = c(0, 2, 4, 6) * 4, labels = c("Chemotaxis","Th1/2 cell diff", "Tissue repair","Inflammation"), line = 6, tick = FALSE)
abline(v = c(0, Xgroups, 23) * 4, lty = 3, col = "black")
abline(h = 0, lty = 2)
legend("bottomright", legend = c("sPLS-DA", "sgPLS-DA"),
       lty = 1, lwd = 3, col = c(MyPal[6], MyPal[10]),
       cex = 0.75)
dev.off()

# Stability analyses 
set.seed(1)
Stability_results = StabilityPlot(X = X_pooled, Y = data$case,
                                  NIter = 100)
pheatmap(Stability_results, cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, filename = "Comp_PLS_stability.pdf",
         height = 5, width = 10)

# ------ severity score calculation
# data preparation 
var = c("eosinophil_sputum_percentage","neutrophil_sputum_percentage","fev1_percentage",
        "exacerbation_per_year","severe_exacerbation_per_year","acq5","ocs","age_of_onset",
        "icu","icu_last_year","atopy","nasal_polyps","Smoking.Status","Second.Hand.Smoke","BMI","case","id")
data_score <-  data[,var]

categorical_cols <- c("ocs", "nasal_polyps","Smoking.Status","Second.Hand.Smoke")
mydata_categorical <- data_score[categorical_cols]
encoded_data <- model.matrix(~ . - 1, data = mydata_categorical)
data_score <- cbind(data_score, encoded_data)
data_score <- data_score[,-c(7,12:14)]
names(data_score)[names(data_score) == "Second.Hand.Smokeyes"] <- "Second_hand_smoke"

# Create a response variable (e.g., a severity score)
# non-penalised PLS-DA model
X_var <- data_score[, -which(colnames(data_score) == "case")]
X_var <- data_score[,-which(colnames(data_score) == "id")].   #keepi id for further use
PLSDA_score <- plsda(X_var,data_score$case,ncomp = 1)
PLSDA_score$loadings$X
PLSDA_score$loadings$Y
PLSDA_score$explained_variance

# calculating severity score
# Extract PLS coefficients
coef <- as.vector(PLSDA_score$loadings$X)
data_score$severity_score <- NA
severity_score <- as.matrix(data_score[c("eosinophil_sputum_percentage", "neutrophil_sputum_percentage", "fev1_percentage", "exacerbation_per_year", "severe_exacerbation_per_year", 
                                         "acq5", "ocsat_least_once_per_day","ocsless_than_once_per_day", "ocsnever","ocsPrevious","age_of_onset", "icu", "icu_last_year", 
                                         "atopy", "nasal_polypsuncertain", "nasal_polypsyes","Smoking.Statusnon_smoker","Second.Hand.Smokeuncertain",
                                         "Smoking.Statusex_smoker", "Second_hand_smoke", "BMI")]) %*% coef
# Assign severity score to corresponding row in new column
data_score$severity_score <- severity_score

pdf("his_severityscore.pdf")
plot <- hist(data_score$severity_score, breaks = 10, col = "blue", xlab = "Severity score", ylab = "Frequency",main = "Histogram of severity scores")
dev.off()

# if we want to make the severity score to all positive 
pdf("his_severityscore_shifted.pdf")
min_score <- min(severity_score)
shifted_scores <- severity_score - min_score
dev.off()

# ----- Not necessary
# running an uncalibrated sparse PLS-DA model
sPLSDA_score <- splsda(X_var,data_score$case, ncomp = 1,
                          keepX = 5)
sPLSDA_score$loadings$X
sPLSDA_score$loadings$X[sPLSDA_score$loadings$X != 0, ]
