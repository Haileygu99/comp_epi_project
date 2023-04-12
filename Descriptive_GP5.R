# Pre-processing code: 
  rm(list=ls()) 
library(tidyverse) 
library(openxlsx) 

setwd('/rds/general/project/hda-22-23/live/Comp_Epi/General/Group5') 
baseline_clinical <- readRDS('./Data/baseline_data_clinical.rds') 
baseline_covariates <- readRDS('./Data/baseline_data_covariates.rds') 
longitudinal_clinical <- readRDS('./Data/longitudinal_data_clinical.rds') 

#data_clinical <- read_tsv(file = './data/data_clinical.tsv') 

length(unique(baseline_clinical$id)) == nrow(baseline_clinical) 
length(unique(baseline_covariates$id)) == nrow(baseline_covariates) 
length(unique(longitudinal_clinical$id)) == nrow(longitudinal_clinical) 

baseline_all <- merge(x=baseline_clinical,y=baseline_covariates,  
                      by="id", all.x=TRUE) 

baseline_all_cohort_A <- subset(baseline_all, cohort == 'cohort_a') 
baseline_all_cohort_B <- subset(baseline_all, cohort == 'cohort_b') 
baseline_all_cohort_C <- subset(baseline_all, cohort == 'cohort_c') 
baseline_all_cohort_D <- subset(baseline_all, cohort == 'cohort_d') 

num_cohorts <- c(nrow(baseline_all_cohort_A), nrow(baseline_all_cohort_B), nrow(baseline_all_cohort_C), nrow(baseline_all_cohort_D)) 

barplot(num_cohorts, names=c('A', 'B', 'C', 'D'), xlab = 'Cohort', ylab = 'Number of Participants') 
# A = severe nonsmoking asthma (SAn) 
# B = smokers and ex-smokers with severe asthma (SAs/ex) 
# C = mild/moderate nonsmoking asthmatics (MMA) 
# D = healthy nonsmoking controls (HC) 

# Proteins ---- 
getwd() 
proteins = read.csv("Data/UB1st_submit20220729.csv") 
#baseline = readRDS("../data/baseline_data_clinical.rds") 
sample_key = read.xlsx("Data/Sample_key_16Aug2022.xlsx") 

dim(sample_key) 
head(proteins) 


sample_key = sample_key[11:18] 
mydata = merge(sample_key[c("KIT.ID", "TranSMART.Allocation.ID","Comment","PLATE.ID")], proteins, by="KIT.ID")
table(table(mydata$TranSMART.Allocation.ID)) 

newdata = mydata %>% mutate(id = TranSMART.Allocation.ID) %>% select(-TranSMART.Allocation.ID) %>% 
  merge(baseline_all[c("id","cohort")], by = "id") 
table(table(newdata$id)) # 299 participants with repeating measurements 
table(newdata$cohort[!duplicated(newdata$id)]) 

proteins_baseline <- subset(newdata, Comment=='Adult baseline') 
baseline_all_no_proteins <- baseline_all 
baseline_all <- merge(baseline_all_no_proteins, proteins_baseline, by='id') 

# Creating case and control column - 0 control; 1 case
baseline_all$case<- ifelse(baseline_all$cohort.x %in% c("cohort_a", "cohort_b","cohort_c"), "1", "0")
table(baseline_all$case)

# ------------- Recoding 
# Education #secondary_or_lower=1, university=2, higher=3 
baseline_all$Education <- ifelse(baseline_all$`Highest Level Education` %in% c("COMPLETE SECONDARY SCHOOL","High school","PLUMBING AND HEATING ENGINEER QUALIFICATIONS","course in bakery","INDUSTRY QUALIFICATION","PROFESSION QUALIFICATIONS","not_complete_secondary", "no_formal_education"), "Secondary_or_lower",
                                 ifelse(baseline_all$`Highest Level Education` %in% c("higher_or_professional", "POST GRADUATE"), "Higher",
                                        ifelse(baseline_all$`Highest Level Education` %in% c("undergraduate_universuty", "post_secondary"), "University",NA)))
table(baseline_all$Education)

baseline_all$Education <- ifelse(baseline_all$Education == "Secondary_or_lower", 1,
                           ifelse(baseline_all$Education =="University", 2, 
                            ifelse(baseline_all$Education =="Higher", 3, NA)))

#plot for descriptive
ggplot(df_edu_1, aes(x = Education, fill = case)) +
  geom_bar(alpha = 0.5) +
  labs(x = "Education", y = "Density") +
  ggtitle("Distribution Plot of Education by Case/Control Status") +
  scale_fill_manual(values = c("blue", "red")) +
  guides(fill = guide_legend(title = "Status"))

# Occupation - group volunteer to other (volunteer = 8 ppl + other =44)
baseline_all$Occupation <- ifelse(baseline_all$`occupational employed` %in% c("1"), "employed",
                                  ifelse(baseline_all$`occupational unemployed` %in% c("1"), "unemployed",
                                         ifelse(baseline_all$`occupational retired` %in% c("1"), "retired",
                                                ifelse(baseline_all$`occupational keeping house` %in% c("1"), "housekeeping",
                                                       ifelse(baseline_all$`occupational student`%in% c("1"), "student",
                                                              ifelse(baseline_all$`occupational volunteer` %in% c("1"), "other",
                                                                     ifelse(baseline_all$`occupational other` %in% c("1"), "other",NA)))))))

# Recoding occcupation to integers 1= employed, housekeeping=2, other=3, retired=4, student=5, unemployed=6
baseline_all$Occupation <- ifelse(baseline_all$Occupation == "employed", 1,
                                 ifelse(baseline_all$Occupation =="housekeeping", 2, 
                                        ifelse(baseline_all$Occupation =="other", 3, 
                                               ifelse(baseline_all$Occupation =="retired", 4, 
                                                      ifelse(baseline_all$Occupation =="student", 5, 
                                                             ifelse(baseline_all$Occupation =="unemployed",6, NA))))))
table(baseline_all$Occupation)

# Ethnic_father - White_caucasian=1, Asian (2)= central_asian +east_asian+south_asian+south_east_asian, black_african (3), other (4)= arabic_north_heritage (6 ppl) + other (6 ppl) + uncertain (2 ppl)
baseline_all$Ethnic_father <- ifelse(baseline_all$`Ethnic origin Father` %in% c("white_caucasian"), "white_caucasian",
                                 ifelse(baseline_all$`Ethnic origin Father` %in% c("central_asian", "east_asian","south_asian","south_east_asian"), "asian",
                                  ifelse(baseline_all$`Ethnic origin Father`%in% c("black_african"), "black_african",
                                          ifelse(baseline_all$`Ethnic origin Father`%in% c("arabic_north_heritage","other","uncertain"), "other", NA))))
baseline_all$Ethnic_father <- ifelse(baseline_all$Ethnic_father == "white_caucasian",1,
                                 ifelse(baseline_all$Ethnic_father =="asian",2, 
                                        ifelse(baseline_all$Ethnic_father =="black_african",3, 
                                         ifelse(baseline_all$Ethnic_father =="other",4,NA))))
table(baseline_all$Ethnic_father)

# Ethnic_mother - White_caucasian=1, Asian (2)= central_asian + east_asian + south_asian + south_east_asian, black_african (3), other (4)= arabic_north_heritage (8 ppl) + other (5 ppl) + uncertain (2 ppl) + multiple_races (1 person)
baseline_all$Ethnic_mother <- ifelse(baseline_all$`Ethnic origin Mother` %in% c("white_caucasian"), "white_caucasian",
                                     ifelse(baseline_all$`Ethnic origin Mother` %in% c("central_asian", "east_asian","south_asian","south_east_asian"), "asian",
                                            ifelse(baseline_all$`Ethnic origin Mother`%in% c("black_african"), "black_african",
                                                   ifelse(baseline_all$`Ethnic origin Mother`%in% c("arabic_north_heritage","other","uncertain","multiple_races"), "other", NA))))
baseline_all$Ethnic_mother <- ifelse(baseline_all$Ethnic_mother == "white_caucasian",1,
                                     ifelse(baseline_all$Ethnic_mother =="asian",2, 
                                            ifelse(baseline_all$Ethnic_mother =="black_african",3, 
                                                   ifelse(baseline_all$Ethnic_mother =="other",4,NA))))
table(baseline_all$Ethnic_mother)

# Race - white_caucasia(1), asian(2) = central_asian + east_asian + south_asian+ south_east_asian, black_african (3), other (4) = arabic_north_heritage (5 ppl) + other (7 ppl) + multiple_races ( 7 ppl)
baseline_all$Race <- ifelse(baseline_all$Race %in% c("white_caucasian"), "white_caucasian",
                                     ifelse(baseline_all$Race %in% c("central_asian", "east_asian","south_asian","south_east_asian"), "asian",
                                            ifelse(baseline_all$Race %in% c("black_african"), "black_african",
                                                   ifelse(baseline_all$Race %in% c("arabic_north_heritage","other","multiple_races"), "other", NA))))
baseline_all$Race <- ifelse(baseline_all$Race == "white_caucasian",1,
                                     ifelse(baseline_all$Race =="asian",2, 
                                            ifelse(baseline_all$Race =="black_african",3, 
                                                   ifelse(baseline_all$Race =="other",4,NA))))

##################################### 
# icu_times: 
sum(is.na(baseline_all$icu_times)) # 504 NA's: not sure what to do with these 

# icu - no/uncertainty=0, yes=1: 
sum(is.na(baseline_all$icu)) # 100 NA's 
baseline_all$icu = as.character(baseline_all$icu) 
baseline_all$icu = ifelse(baseline_all$icu == 'yes', 'yes', 'no/uncertain') 
baseline_all$icu = as.factor(baseline_all$icu) 
baseline_all$icu <- ifelse(baseline_all$icu == "no/uncertain",0,
                             ifelse(baseline_all$icu =="yes",1, NA))
table(baseline_all$icu) 

#ocs: 
# Combine to : at least once per day=3, less than once per day=2, previous=1, never=0
baseline_all$ocs = as.character(baseline_all$ocs) 
baseline_all$ocs = ifelse(baseline_all$ocs == "daily", "at_least_once_per_day", baseline_all$ocs) 
baseline_all$ocs = ifelse(baseline_all$ocs == "at_least_twice_a_day", "at_least_once_per_day", baseline_all$ocs) 
baseline_all$ocs = ifelse(baseline_all$ocs == "weekly", "less_than_once_per_day", baseline_all$ocs) 
baseline_all$ocs = ifelse(baseline_all$ocs == "more_than_twice_week", "less_than_once_per_day", baseline_all$ocs) 
baseline_all$ocs = ifelse(baseline_all$ocs == "more_than_once_a_month", "less_than_once_per_day", baseline_all$ocs) 
baseline_all$ocs = ifelse(baseline_all$ocs == "once_a_month", "less_than_once_per_day", baseline_all$ocs) 
baseline_all = subset(baseline_all, is.na(baseline_all$ocs) | baseline_all$ocs!= 'Free_text') 
baseline_all$ocs = as.factor(baseline_all$ocs) 
baseline_all$ocs <- ifelse(baseline_all$ocs == "at_least_once_per_day",3,
                     ifelse(baseline_all$ocs =="less_than_once_per_day",2,
                            ifelse(baseline_all$ocs =="never",0,
                                   ifelse(baseline_all$ocs =="Previous",1,NA))))
table(baseline_all$ocs) 

# atopy: remove uncertain - negative=0, positive=1
sum(is.na(baseline_all$atopy)) # No NA's 
baseline_all$atopy = as.character(baseline_all$atopy) 
baseline_all = subset(baseline_all, atopy != 'uncertain') 
baseline_all$atopy = as.factor(baseline_all$atopy) 
table(baseline_all$atopy) 
baseline_all$atopy <- ifelse(baseline_all$atopy == "negative",0,
                                     ifelse(baseline_all$atopy =="positive",1, NA))


# icu_last_year: remove uncertain - no = 0, yes = 1
sum(is.na(baseline_all$icu_last_year)) # 100 NA's 
baseline_all$icu_last_year = as.character(baseline_all$icu_last_year) 
baseline_all = subset(baseline_all, is.na(icu_last_year) | icu_last_year != "uncertain" ) 
table(baseline_all$icu_last_year) 
baseline_all$icu_last_year <- ifelse(baseline_all$icu_last_year == "no",0,
                            ifelse(baseline_all$icu_last_year =="yes",1, NA))

# removing useful columns 
baseline_all <- baseline_all[, -c(28:29,31:39)]

# Protein Descriptions ---------------------------------------------------- 
library(tableone) 
library(corrplot) 
proteins<-c("CCL18","IL18","CCL5","i.TAC","TARC","IL.10","IL.16", "MIP.1a","PAPP.A","IL17","TSLP","SP.A", "CTACK","KL.6","MPO","EDN","IL.6","TNF.A", "MIG","CCL20","IP.10","IL.4","IL.5") 
binary<-c("occupational employed","occupational unemployed","occupational retired","occupational keeping house","occupational student","occupational volunteer", 
          "occupational other") 
# Define a vector of column names to apply the outlier detection and removal to 
numeric_column_names <- colnames(select_if(baseline_all, is.numeric)) 
numeric_column_names<-numeric_column_names[!(numeric_column_names %in% binary)] 
numeric_column_names<-numeric_column_names[!(numeric_column_names %in% c('KIT.ID'))] 
non_protein_cols<-numeric_column_names[!(numeric_column_names %in% proteins)] 
characterize_numeric<- function(df, variables){  
  #Table 1 
  df$'case'<-ifelse(df$cohort.x=='cohort_d', "control", "case") 
  summary(df) 
  table1<-CreateTableOne(data = df[-1],strata = 'case',vars=variables) 
  print(table1) 
  # Create a long-format data frame with the protein data 
  baseline_long <- tidyr::pivot_longer(df, cols = variables, names_to = "protein") 
  # Create a histogram of protein expression, faceted by protein name 
  histogram<-ggplot(baseline_long, aes(x = value,fill=case)) + 
    geom_density(alpha=0.5) + 
    facet_wrap(~ protein, scales = "free") + 
    xlab("Protein Expression") + 
    ylab("Count") + 
    ggtitle("Histograms of Protein Expression") 
  # Load the corrplot package 
  baseline_proteins<-df[,variables] 
  sum(complete.cases(df[, c(variables)])) 
  cor_matrix <- cor(df[,variables], use = "pairwise.complete.obs") 
  correlation_plot<-corrplot(cor_matrix, method="color") 
  return(list(histogram,correlation_plot)) 
} 

characterize_numeric(baseline_all, proteins) #Proteins before log transformation 

#Log transform proteins 
# Define a custom function to add half the minimum positive value to each column 
add_half_min_positive <- function(column) { 
  min_positive <- min(column[column > 0], na.rm = TRUE) 
  return(column + min_positive / 2) 
} 
baseline_all<-mutate(baseline_all,across(proteins, ~log(add_half_min_positive(.)))) 

characterize_numeric(baseline_all, proteins) #Proteins after log transformation 
#Other numeric variables 
characterize_numeric(baseline_all, non_protein_cols)  
cols_to_log<-c('acq5','Body Mass Index (kg/m2)','eosinophil_sputum_percentage','eosinophil_serum_count','exacerbation_per_year','neutrophil_serum_count','Pack Years','severe_exacerbation_per_year') 
baseline_all<-mutate(baseline_all,across(cols_to_log, ~log(add_half_min_positive(.)))) 
characterize_numeric(baseline_all, non_protein_cols)  

# Outlier Removal ---------------------------------------------------- 
# Initialize a vector to store row indices of outliers 
outlier_rows <- c() 

# Remove proteins with outliers, but not for proteins where the outliers are mostly equal to some specific minimum value 
outlier_vars<-setdiff(proteins, "TARC") 
outlier_vars<-setdiff(outlier_vars, "i.TAC") 
outlier_vars<-setdiff(outlier_vars, "IL.4") 
for (col in outlier_vars) { 
  # Identify outliers in the current column, excluding missing values 
  col_no_na <- baseline_all[[col]][!is.na(baseline_all[[col]])] 
  mean_col <- mean(col_no_na) 
  sd_col <- sd(col_no_na) 
  outliers <- col_no_na[col_no_na < mean_col - 3*sd_col | col_no_na > mean_col + 3*sd_col] 
  print(col) 
  print(outliers) 
  # Get the row indices of the outliers in the current column 
  outlier_indices <- which(!is.na(baseline_all[[col]]) & baseline_all[[col]] %in% outliers) 
  # Add the outlier indices to the outlier_rows vector 
  outlier_rows <- unique(c(outlier_rows, outlier_indices)) 
} 
# Remove rows with outliers from the baseline_all dataframe 
baseline_all <- baseline_all[-outlier_rows, ] 

# saving clean data
write.csv(baseline_all, "baseline_all_clean.csv")

