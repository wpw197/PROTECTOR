## ------------------------------------------------------------------------
## 
## Script name: Design Option 2 Sample Size Calculation
##
## Author: Emily Lane
##
## Date Created: 18/10/23
##
## ------------------------------------------------------------------------
##
## Notes: This file contains a function that simulates the number of ovarian
##   cancer diagnoses in a cohort of BRCA1 and BRCA2 patients (conditional on
##   age) and calculates the power of a study with a specified sample size.
##
## ------------------------------------------------------------------------
## load up the packages we will need: 

packages <- c("tidyverse", "stringr", "readxl")
install.packages(setdiff(packages, rownames(installed.packages())))
lapply(packages, library, character.only = TRUE, quietly = TRUE)
rm(packages)

## ------------------------------------------------------------------------
BRCA1_Risks <- read_excel('Chen Adjusted BOADICEA Risks.xlsx',1)
BRCA2_Risks <- read_excel('Chen Adjusted BOADICEA Risks.xlsx',2)
BRCA1_Ages <- read.csv('BRCA1_Ages.csv')
BRCA2_Ages <- read.csv('BRCA2_Ages.csv')

################################################################################
## Inputs:
## meno_age: age at menopause/censoring
## Sample_Size: total number of patients undergoing surgical intervention
## follow_up_years: number of years after surgery patients are followed up
## Risk_reduction: the expected reduction in ovarian cancer risk from undergoing surgery
## Margin: relative non-inferiority (observed vs. expected) used in the power calculation
## nSims: number of times the study is simulated

OC_Simulation <- function(meno_age = 51, Sample_Size, follow_up_years, Risk_reduction1, Risk_reduction2, Margin = 0.8, nSims = 1000){
  
  ## Checking all inputs 
  if(Sample_Size < 2){
    message("Minimum sample size of 2")
    stop()
  }else if(follow_up_years %% 1 != 0){
    follow_up_years <- ceiling(follow_up_years)
    warning(paste("Follow-up years rounded to", follow_up_years, sep = " "))
  }else if(!(Risk_reduction1 < 1 & Risk_reduction1 > 0)){
    message("Risk reduction must be between 0 and 1")
    stop()
  }else if(!(Risk_reduction2 < 1 & Risk_reduction2 > 0)){
    message("Risk reduction must be between 0 and 1")
    stop()
  }
  
  n <- nSims
  ## Adjusting for drop out - 5%
  Sample_Size <- round(Sample_Size*0.95)
  
  ## Assuming 1:1 BRCA1:BRCA2 gene mutation
  N_BRCA1 <- round(Sample_Size/2)
  N_BRCA2 <- round(Sample_Size/2)
  BRCA1_Ages <- BRCA1_Ages %>% 
    filter(Age <= meno_age - 1.5)
  BRCA2_Ages <- BRCA2_Ages %>% 
    filter(Age <= meno_age - 1.5)
  
  ## output data frame
  My_Data <- data.frame(BRCA1_N = rep(NA, n) , BRCA2_N = rep(NA, n), RRSO_BRCA1 = rep(NA, n), RRESDO_BRCA1 = rep(NA, n), 
                        RRSO_BRCA2 = rep(NA, n), RRESDO_BRCA2 = rep(NA, n), Mean_Age_BRCA1 = rep(NA, n), 
                        Mean_Age_BRCA2 = rep(NA, n), Power = rep(NA, n))
  
  
  
  
}

Data_Generator <- function(meno_age, n_sims, n_sample, follow_up, RRES_reduction){
  
  RRES_reduction <- 1-RRES_reduction
  
  sample_size <- data.frame(sample_size = n_sample*0.95)
  sample_size <- sample_size %>% 
    mutate(n_BRCA1 = round(sample_size*BRCA1_prop), n_BRCA2 = round(sample_size*BRCA2_prop))
  sample_size <- sample_size %>% 
    mutate(n_BRCA1_Control = round(n_BRCA1*BRCA_ARM_Nums %>% 
                                     filter(BRCA == 'BRCA1') %>% 
                                     filter(ARM_NAME == 'Control') %>% 
                                     select(prop) %>% 
                                     pull()), 
           n_BRCA1_RRES = round(n_BRCA1*BRCA_ARM_Nums %>% 
                                  filter(BRCA == 'BRCA1') %>% 
                                  filter(ARM_NAME == 'Risk-reducing early salpingectomy delayed oophorectomy (RRESDO)') %>% 
                                  select(prop) %>% 
                                  pull()), 
           n_BRCA1_RROS = round(n_BRCA1*BRCA_ARM_Nums %>% 
                                  filter(BRCA == 'BRCA1') %>% 
                                  filter(ARM_NAME == 'Risk-reducing salpingo-oophorectomy (RRSO)') %>% 
                                  select(prop) %>% 
                                  pull()), 
           n_BRCA2_Control = round(n_BRCA1*BRCA_ARM_Nums %>% 
                                     filter(BRCA == 'BRCA2') %>% 
                                     filter(ARM_NAME == 'Control') %>% 
                                     select(prop) %>% 
                                     pull()), 
           n_BRCA2_RRES = round(n_BRCA1*BRCA_ARM_Nums %>% 
                                  filter(BRCA == 'BRCA2') %>% 
                                  filter(ARM_NAME == 'Risk-reducing early salpingectomy delayed oophorectomy (RRESDO)') %>% 
                                  select(prop) %>% 
                                  pull()), 
           n_BRCA2_RROS = round(n_BRCA1*BRCA_ARM_Nums %>% 
                                  filter(BRCA == 'BRCA2') %>% 
                                  filter(ARM_NAME == 'Risk-reducing salpingo-oophorectomy (RRSO)') %>% 
                                  select(prop) %>% 
                                  pull()))
  
  BRCA1_Ages <- BRCA1_Ages %>% 
    filter(Age <= meno_age)
  
  BRCA2_Ages <- BRCA2_Ages %>% 
    filter(Age <= meno_age)
  
  BRCA1_RRSO_Data <- data.frame(id = c(1:sample_size[1,6]), Age = sample(BRCA1_Ages$Age, sample_size[1,6], replace = TRUE), 
                                Menopause_Age = c(rep(meno_age, sample_size[1,6])))
  BRCA1_RRES1_Data <- data.frame(id = c(1:sample_size[1,5]), Age = sample(BRCA1_Ages$Age, sample_size[1,5], replace = TRUE),
                                 Menopause_Age = c(rep(meno_age, sample_size[1,5])))
  
  BRCA2_RRSO_Data <- data.frame(id = c(1:sample_size[1,9]), Age = sample(BRCA2_Ages$Age, sample_size[1,9], replace = TRUE),
                                Menopause_Age = c(rep(meno_age, sample_size[1,9])))
  BRCA2_RRES1_Data <- data.frame(id = c(1:sample_size[1,8]), Age = sample(BRCA2_Ages$Age, sample_size[1,8], replace = TRUE),
                                 Menopause_Age = c(rep(meno_age, sample_size[1,8])))
  
  ## Risk reduction selection
  
  BRCA1_RRSO_Data$Age_at_Exit <- ifelse(BRCA1_RRSO_Data$Age + follow_up < meno_age, BRCA1_RRSO_Data$Age + follow_up, meno_age)
  BRCA1_RRES1_Data$Age_at_Exit <- ifelse(BRCA1_RRES1_Data$Age + follow_up < meno_age, BRCA1_RRES1_Data$Age + follow_up, meno_age)
  BRCA2_RRES1_Data$Age_at_Exit <- ifelse(BRCA2_RRES1_Data$Age + follow_up < meno_age, BRCA1_RRES1_Data$Age + follow_up, meno_age)
  BRCA2_RRSO_Data$Age_at_Exit <- ifelse(BRCA2_RRSO_Data$Age + follow_up < meno_age, BRCA2_RRSO_Data$Age + follow_up, meno_age)
  
  ## Adding cumulative Hazard at study entry and exit:
  
  BRCA1_RRSO_Data <- BRCA1_RRSO_Data %>% 
    mutate(Age_At_Entry = Age, Age = Age - 1) %>% 
    left_join(BRCA1_Risks %>% 
                select(Age, Cumulative_Hazard), by = 'Age') %>% 
    left_join(BRCA1_Risks %>%
                select(Age, Cumulative_Hazard) %>% 
                rename(Age_at_Exit = Age), by = 'Age_at_Exit') %>% 
    rename(Cum_Hazard_Entry = Cumulative_Hazard.x, Cum_Hazard_Exit = Cumulative_Hazard.y) %>% 
    mutate(Study_Hazard = (Cum_Hazard_Exit*0.04)-(Cum_Hazard_Entry*0.04)) %>% 
    mutate(Expected_Hazard = Cum_Hazard_Exit-Cum_Hazard_Entry) 
  
  ## BRCA2 
  
  BRCA2_RRSO_Data <- BRCA2_RRSO_Data %>% 
    mutate(Age_At_Entry = Age, Age = Age - 1) %>% 
    left_join(BRCA2_Risks %>% 
                select(Age, Cumulative_Hazard), by = 'Age') %>% 
    left_join(BRCA2_Risks %>%
                select(Age, Cumulative_Hazard) %>% 
                rename(Age_at_Exit = Age), by = 'Age_at_Exit') %>% 
    rename(Cum_Hazard_Entry = Cumulative_Hazard.x, Cum_Hazard_Exit = Cumulative_Hazard.y) %>% 
    mutate(Study_Hazard = (Cum_Hazard_Exit*0.04)-(Cum_Hazard_Entry*0.04)) %>% 
    mutate(Expected_Hazard = Cum_Hazard_Exit-Cum_Hazard_Entry)  
  
  ## RRESDO BRCA1
  
  BRCA1_RRES1_Data <- BRCA1_RRES1_Data %>% 
    mutate(Age_At_Entry = Age, Age = Age - 1) %>% 
    left_join(BRCA1_Risks %>% 
                select(Age, Cumulative_Hazard), by = 'Age') %>% 
    left_join(BRCA1_Risks %>%
                select(Age, Cumulative_Hazard) %>% 
                rename(Age_at_Exit = Age), by = 'Age_at_Exit') %>% 
    rename(Cum_Hazard_Entry = Cumulative_Hazard.x, Cum_Hazard_Exit = Cumulative_Hazard.y) %>% 
    mutate(Study_Hazard = (Cum_Hazard_Exit*RRES_reduction)-(Cum_Hazard_Entry*RRES_reduction)) %>% 
    mutate(Expected_Hazard = Cum_Hazard_Exit-Cum_Hazard_Entry) 
  
  ## BRCA 2
  
  BRCA2_RRES1_Data <- BRCA2_RRES1_Data %>% 
    mutate(Age_At_Entry = Age, Age = Age - 1) %>% 
    left_join(BRCA2_Risks %>% 
                select(Age, Cumulative_Hazard), by = 'Age') %>% 
    left_join(BRCA2_Risks %>%
                select(Age, Cumulative_Hazard) %>% 
                rename(Age_at_Exit = Age), by = 'Age_at_Exit') %>% 
    rename(Cum_Hazard_Entry = Cumulative_Hazard.x, Cum_Hazard_Exit = Cumulative_Hazard.y) %>% 
    mutate(Study_Hazard = (Cum_Hazard_Exit*RRES_reduction)-(Cum_Hazard_Entry*RRES_reduction)) %>% 
    mutate(Expected_Hazard = Cum_Hazard_Exit-Cum_Hazard_Entry) 
  
  #############
  
  ## BRCA 1 
  
  dat_frame <- data.frame(matrix(NA, ncol = n_sims, nrow = dim(BRCA1_RRES1_Data)[1]))
  colnames(dat_frame) <- paste0('Binom_Trial_', 1:n_sims, by = '')
  
  BRCA1_RRES1_Data <- BRCA1_RRES1_Data %>% 
    bind_cols(dat_frame)
  
  
  dat_frame <- data.frame(matrix(NA, ncol = n_sims, nrow = dim(BRCA1_RRSO_Data)[1]))
  colnames(dat_frame) <- paste0('Binom_Trial_', 1:n_sims, by = '')
  
  BRCA1_RRSO_Data <- BRCA1_RRSO_Data %>% 
    bind_cols(dat_frame)
  
  ## BRCA 2 
  
  dat_frame <- data.frame(matrix(NA, ncol = n_sims, nrow = dim(BRCA2_RRES1_Data)[1]))
  colnames(dat_frame) <- paste0('Binom_Trial_', 1:n_sims, by = '')
  
  BRCA2_RRES1_Data <- BRCA2_RRES1_Data %>% 
    bind_cols(dat_frame)
  
  
  dat_frame <- data.frame(matrix(NA, ncol = n_sims, nrow = dim(BRCA2_RRSO_Data)[1]))
  colnames(dat_frame) <- paste0('Binom_Trial_', 1:n_sims, by = '')
  
  BRCA2_RRSO_Data <- BRCA2_RRSO_Data %>% 
    bind_cols(dat_frame)
  
  ############
  
  Binomial <- function(data){
    
    for(i in 1:dim(data)[1]){
      
      data[i, 10:(9+n_sims)] <- rbinom(n_sims, 1, prob = data$Study_Hazard[i])
      
      
    }
    return(data)
    
  }
  
  BRCA1_RRSO_Data <- Binomial(BRCA1_RRSO_Data)
  BRCA1_RRES1_Data <- Binomial(BRCA1_RRES1_Data)
  BRCA2_RRSO_Data <- Binomial(BRCA2_RRSO_Data)
  BRCA2_RRES1_Data <- Binomial(BRCA2_RRES1_Data)
  
  My_Data <- data.frame(Simulation = 1:n_sims, RRES_BRCA1_N = NA , RRES_BRCA2_N = NA, RRES_Expected_BRCA1 = NA, RRES_Observed_BRCA1 = NA, 
                        RRES_Expected_BRCA2 = NA, RRES_Observed_BRCA2 = NA, 
                        RRSO_BRCA1_N = NA, RRSO_BRCA2_N = NA,
                        RRSO_Expected_BRCA1 = NA, RRSO_Observed_BRCA1 = NA,
                        RRSO_Expected_BRCA2 = NA, RRSO_Observed_BRCA2 = NA, Lower_CI = NA, 
                        Upper_CI = NA, Estimate = NA)
  
  for(i in 1:n_sims){
    
    My_Data[i, 2] <- dim(BRCA1_RRES1_Data)[1]
    My_Data[i, 3] <- dim(BRCA2_RRES1_Data)[1]
    My_Data[i, 4] <- sum(BRCA1_RRES1_Data$Expected_Hazard)
    My_Data[i, 5] <- sum(BRCA1_RRES1_Data[, i + 9])
    My_Data[i, 6] <- sum(BRCA2_RRES1_Data$Expected_Hazard)
    My_Data[i, 7] <- sum(BRCA2_RRES1_Data[, i + 9])
    My_Data[i, 8] <- dim(BRCA1_RRSO_Data)[1]
    My_Data[i, 9] <- dim(BRCA2_RRSO_Data)[1]
    My_Data[i, 10] <- sum(BRCA1_RRSO_Data$Expected_Hazard)
    My_Data[i, 11] <- sum(BRCA1_RRSO_Data[, i + 9])
    My_Data[i, 12] <- sum(BRCA2_RRSO_Data$Expected_Hazard)
    My_Data[i, 13] <- sum(BRCA2_RRSO_Data[, i + 9])
    
  }
  
  for(i in 1:n_sims){
    
    E2_E1 <-(My_Data$RRES_Expected_BRCA1[i] + My_Data$RRES_Expected_BRCA2[i])/(My_Data$RRSO_Expected_BRCA1[i] + My_Data$RRSO_Expected_BRCA2[i])
    myci <- binom.wilson(My_Data$RRSO_Observed_BRCA1[i] + My_Data$RRSO_Observed_BRCA2[i], My_Data$RRSO_Observed_BRCA1[i] + My_Data$RRSO_Observed_BRCA2[i] +
                           My_Data$RRES_Observed_BRCA1[i] + My_Data$RRES_Observed_BRCA2[i])
    
    My_Data$Lower_CI[i] <- (myci$lower/(1 - myci$lower)) * E2_E1
    My_Data$Upper_CI[i] <- (myci$upper/(1 - myci$upper)) *E2_E1
    My_Data$Estimate[i] <- (myci$proportion/(1 - myci$proportion)) * E2_E1
    
  }
  
  return(My_Data)
  
}


