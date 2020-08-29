#library(foreign)
#library(mvtnorm)
#library(MASS)
#library(sandwich)
#library(lmtest)
#library(ggplot2)
#library(grf)
#library(bindata)
#library(fitdistrplus)
#library(Rfast)
#library(grf)

#setwd("/Users/timohaller/Desktop/CausalForests/Analyse/R")

#rm(list=ls())
#
## Load and select data
#data <- read.csv("data.csv")
#
## Select relevant variables 
#
#depvar <- "expec_stockret_perc_t01"
#covariates <- c("high","female", "atleast_bachelor", "z_conf_prior","prior_t01","prior_squ_t01")
#selection_variables <- c("housing")
#
#data = subset(data,select=c(depvar,covariates,selection_variables))
#
## Rename the variables
#
#names(data) = c("stockret", "high", "female","edu","conf","prior","prior_squ","housing")
#
## Select relevant observations and drop missing values
#
#data = data[data$housing == 0,]
#
## Drop missing values
#
#data =  data[complete.cases(data), ]
#data$housing = NULL
#
#
#n_noise = 20
#N = 1000
#
#set.seed(123)
#gamma_true_val = runif(n=2000,-1,1)
#Sigma_noise_val = diag(2000)

crf_sim_noise <- function(n_noise_vec=c(20,50,100),gamma_true=gamma_true_val, Sigma_noise = Sigma_noise_val, n_obs = 1000, m_runs=10, dgp_complex=T) {
  
  n_sim = length(n_noise_vec)
  
  ### Set up results data_frame
  
  final_table = c()
  
  ### Run simulation study for different n 
  
  for (k in seq(1,n_sim,1)) {
    
    n_noise = as.integer(n_noise_vec[k])
    
    #print(n_noise)
    message(n_noise)
    col = c()
    col = cbind(col,n_noise)
    
    Sigma_noise_sel = Sigma_noise[1:n_noise, 1:n_noise]
    gamma_true_sel = gamma_true[1:n_noise]
    
    beta_true = matrix(data=c(3.5,1.7,0.3,-0.9,-0.5, 0.8,-0.5,0.3,-0.5,0.8,0.15,0.001,0.3,-0.013,-0.000037,0.0001392,gamma_true_sel), nrow=16+length(gamma_true_sel))
    
    set.seed(123)
    N = n_obs
    M = m_runs
    
    counter = 123
    
    
    for (i in seq(1:M)) {
      
      set.seed(counter)
      
      #print(i)
      message(i)
      
      # Prepare dataframes
      
      if (i == 1) {
        
        results_tau_crf = data.frame(x=rep(NA,N))
        
        results_true_Y0 = data.frame(x=rep(NA,N))
        results_true_Y1 = data.frame(x=rep(NA,N))
        results_ATE_crf = data.frame(x=rep(NA,1))
        
      }
      
      # Draw datasample (training and test)
      
      for (j in seq(1:2)) {
        
        set.seed(counter+j)
        
        if (dgp_complex == TRUE) {
          
          constant = rep(1,N)
          epsilon = rnorm(N,sd=4.5)
          
          high = rbinom(N,1,0.5)
          prior = rexp(n=N, rate=0.1176)
          prior_squ = prior*prior
          
          ### Draw noise variables from multivariate distribution
          
          noise_matrix = as.matrix(rmvnorm(n=N,mu = rep(0, n_noise),sigma = Sigma_noise_sel))
          
          ##### Draw female, education and conf from a multivariate distribution ##### 
          
          ### First draw female and education (multivariate binomial distribution, correlation structure taken from data)
          
          commprob = matrix(c(sum(data$edu)/dim(data)[1],dim(data[data$female == 1 & data$edu == 1,])[1]/dim(data)[1]
                              ,dim(data[data$female == 1 & data$edu == 1,])[1]/dim(data)[1],sum(data$female)/dim(data)[1]),nrow=2,ncol=2)
          mprob = c(sum(data$female)/dim(data)[1],sum(data$edu)/dim(data)[1])
          
          sample_edu_fem = data.frame(rmvbin(n=N,margprob=as.matrix(mprob),commonprob=commprob))
          names(sample_edu_fem) = c("female","edu")
          
          ### Let conf depend on gender (correlation structure taken from data)
          
          cat_vector = as.numeric(names(table(data[data$female == 0,]$conf)))
          
          prob_male = c(data.frame(table(data[data$female == 0,]$conf)/length(data[data$female == 0,]$conf))$Freq)
          prob_female = c(data.frame(table(data[data$female == 1,]$conf)/length(data[data$female == 1,]$conf))$Freq)
          
          sample_edu_fem_conf = sample_edu_fem
          sample_edu_fem_conf$conf = NA
          
          for (i in seq(1,N,1)) {
            
            if (sample_edu_fem_conf[i,]$female == 0) {
              
              sample_edu_fem_conf[i,]$conf = sum(rmultinom(n=1,size=1,prob=prob_male)*cat_vector)
            }
            else {
              
              sample_edu_fem_conf[i,]$conf = sum(rmultinom(n=1,size=1,prob=prob_female)*cat_vector)
            }
          }
          
          female = sample_edu_fem_conf$female
          edu = sample_edu_fem_conf$edu
          conf = sample_edu_fem_conf$conf
          
          stockret = 3.5*constant + 1.7*high + 0.3*edu - 0.9*high*edu - 0.5*female + 0.8*high*female - 0.5*conf + 0.3*high*conf -0.5*edu*conf+0.8*high*edu*conf+0.15*prior+0.001*prior_squ+0.3*high*prior-0.013*high*prior_squ-0.000037*prior*prior_squ+0.0001392*high*prior*prior_squ + noise_matrix %*% gamma_true_sel + epsilon  
          
        }
        
        else {
          
          constant = rep(1,N)
          epsilon = rnorm(N,sd=4.5)
          high = rbinom(N,1,0.5)
          
          ### Draw noise variables from multivariate distribution
          
          noise_matrix = as.matrix(rmvnorm(n=N,mu = rep(0, n_noise),sigma = Sigma_noise[1:n_noise,1:n_noise]))
          
          female = rbinom(N,1,0.79)
          edu = rbinom(N,1,0.47)
          prior = rexp(n=N, rate=0.1176)
          prior_squ = prior*prior
          conf = rnorm(N,mean=mean(data$conf),sd=1)
          
          
          
          stockret = 3.5*constant + 1.7*high + 0.3*edu - 0.9*high*edu - 0.5*female + 0.8*high*female - 0.5*conf + 0.3*high*conf -0.5*edu*conf+0.8*high*edu*conf+0.15*prior+0.001*prior_squ+0.3*high*prior-0.013*high*prior_squ-0.000037*prior*prior_squ+0.0001392*high*prior*prior_squ + noise_matrix %*% gamma_true_sel + epsilon  
        }
        
        if (j == 1) {
          
          # Set up training sample for crf
          
          covariates_train = cbind(data.frame(edu,female,conf,prior,prior_squ),data.frame(noise_matrix))
          high_train = high
          stockret_train = stockret
          
          # Fit CRF with training data
          
          cf_fitting <- causal_forest(
            
            X <- as.matrix(covariates_train),
            Y <- as.matrix(stockret_train),
            W <- as.matrix(high_train),
            
            
            seed <- 123,
            #tune.parameters = "all",
            num.trees = 2000
            
          )
        }
        
        else {
          
          # Set up test sample for crf 
          
          covariates_test = cbind(data.frame(edu,female,conf,prior,prior_squ),data.frame(noise_matrix))
          high_test = high
          stockret_test = stockret   
          
          # Generate dataframe for true tau
          
          for (high_val in range(0,1)) {
            
            high = rep(high_val,N)  
            
            high_edu = high*edu
            high_female = high*female
            high_conf = high*conf
            edu_conf = edu*conf
            high_edu_conf = high*edu*conf
            high_prior = high*prior
            high_priorsqu = high*prior_squ
            prior_priorsqu = prior*prior_squ
            high_prior_priorsqu = high*prior*prior_squ
            
            dmatrix_true_test = cbind(data.frame(constant,high,edu,high_edu,female,high_female,conf,high_conf,edu_conf,high_edu_conf,prior,prior_squ,high_prior,high_priorsqu,prior_priorsqu,high_prior_priorsqu),data.frame(noise_matrix))
            
            if (high_val == 0) {
              
              Y0_true = as.matrix(dmatrix_true_test) %*% beta_true
              
            }
            else {
              
              Y1_true = as.matrix(dmatrix_true_test) %*% beta_true
              
            }
          }
        }
        
        
      }      
      
      # CRF predictions with test data and append to results 
      
      tau_crf = predict(cf_fitting,newdata=covariates_test)
      results_tau_crf = cbind(results_tau_crf,tau_crf)
      
      # CRF predictions of ATE
      
      ATE_crf = average_treatment_effect(cf_fitting)[1]
      
      # Append Y0_true and Y1_true to results dataframe  
      
      results_true_Y0 = cbind(results_true_Y0, Y0_true) 
      results_true_Y1 = cbind(results_true_Y1, Y1_true) 
      
      counter = counter + 2
      
      results_ATE_crf = cbind(results_ATE_crf, ATE_crf)
    }
    
    results_tau_true = results_true_Y1 - results_true_Y0
    
    ### Calculate CATE EMSE 
      
    MSE_CATE_crf = colMeans((results_tau_crf - results_tau_true)^2)
    EMSE_CATE_crf = mean(MSE_CATE_crf,na.rm = T)
    MSE_CATE_crf_sd = sd(MSE_CATE_crf,na.rm = T)
    
    ### Calculate ATE EMSE 
    
    ATE_true = colMeans(results_true_Y1) - colMeans(results_true_Y0)  
    MSE_ATE_crf = (results_ATE_crf - ATE_true)^2
    EMSE_ATE_crf =  rowMeans(MSE_ATE_crf,na.rm=TRUE)
    MSE_ATE_crf_sd = sd(MSE_ATE_crf, na.rm=TRUE)  
    
    col = rbind(col,EMSE_CATE_crf,MSE_CATE_crf_sd,EMSE_ATE_crf,MSE_ATE_crf_sd)
    
    final_table = cbind(final_table,col)
    
  }
  
  row.names(final_table)[1] = "Noise"
  return(final_table)
  
}

ols_sim_noise <- function(n_noise_vec=c(20,50,100),gamma_true=gamma_true_val, Sigma_noise = Sigma_noise_val, n_obs = 4000, m_runs=10, dgp_complex=T) {
  
  n_sim = length(n_noise_vec)
  
  ### Set up results data_frame
  
  final_table = c()
  
  ### Run simulation study for different n 
  
  for (k in seq(1,n_sim,1)) {
    
    n_noise = as.integer(n_noise_vec[k])
    
    #print(n_noise)
    message(n_noise)
    col = c()
    col = cbind(col,n_noise)
    
    Sigma_noise_sel = Sigma_noise[1:n_noise, 1:n_noise]
    gamma_true_sel = gamma_true[1:n_noise]
    
    beta_true = matrix(data=c(3.5,1.7,0.3,-0.9,-0.5, 0.8,-0.5,0.3,-0.5,0.8,0.15,0.001,0.3,-0.013,-0.000037,0.0001392,gamma_true_sel), nrow=16+length(gamma_true_sel))
    
    set.seed(123)
    N = n_obs
    M = m_runs
    
    counter = 123
    
    
    for (i in seq(1:M)) {
      
      set.seed(counter)
      
      #print(i)
      message(i)
      
      # Prepare dataframes
      
      if (i == 1) {
        
        results_simple_Y0 = data.frame(x=rep(NA,N))
        results_simple_Y1 = data.frame(x=rep(NA,N))
        
        results_complex_Y0 = data.frame(x=rep(NA,N))
        results_complex_Y1 = data.frame(x=rep(NA,N))
        
        results_true_Y0 = data.frame(x=rep(NA,N))
        results_true_Y1 = data.frame(x=rep(NA,N))
        
      }
      
      # Draw datasample (training and test)
      
      for (j in seq(1:2)) {
        
        set.seed(counter+j)
        
        if (dgp_complex == TRUE) {
          
          constant = rep(1,N)
          epsilon = rnorm(N,sd=4.5)
          
          high = rbinom(N,1,0.5)
          prior = rexp(n=N, rate=0.1176)
          prior_squ = prior*prior
          
          ### Draw noise variables from multivariate distribution
          
          noise_matrix = as.matrix(rmvnorm(n=N,mu = rep(0, n_noise),sigma = Sigma_noise_sel))
          
          ##### Draw female, education and conf from a multivariate distribution ##### 
          
          ### First draw female and education (multivariate binomial distribution, correlation structure taken from data)
          
          commprob = matrix(c(sum(data$edu)/dim(data)[1],dim(data[data$female == 1 & data$edu == 1,])[1]/dim(data)[1]
                              ,dim(data[data$female == 1 & data$edu == 1,])[1]/dim(data)[1],sum(data$female)/dim(data)[1]),nrow=2,ncol=2)
          mprob = c(sum(data$female)/dim(data)[1],sum(data$edu)/dim(data)[1])
          
          sample_edu_fem = data.frame(rmvbin(n=N,margprob=as.matrix(mprob),commonprob=commprob))
          names(sample_edu_fem) = c("female","edu")
          
          ### Let conf depend on gender (correlation structure taken from data)
          
          cat_vector = as.numeric(names(table(data[data$female == 0,]$conf)))
          
          prob_male = c(data.frame(table(data[data$female == 0,]$conf)/length(data[data$female == 0,]$conf))$Freq)
          prob_female = c(data.frame(table(data[data$female == 1,]$conf)/length(data[data$female == 1,]$conf))$Freq)
          
          sample_edu_fem_conf = sample_edu_fem
          sample_edu_fem_conf$conf = NA
          
          for (i in seq(1,N,1)) {
            
            if (sample_edu_fem_conf[i,]$female == 0) {
              
              sample_edu_fem_conf[i,]$conf = sum(rmultinom(n=1,size=1,prob=prob_male)*cat_vector)
            }
            else {
              
              sample_edu_fem_conf[i,]$conf = sum(rmultinom(n=1,size=1,prob=prob_female)*cat_vector)
            }
          }
          
          female = sample_edu_fem_conf$female
          edu = sample_edu_fem_conf$edu
          conf = sample_edu_fem_conf$conf
          
          stockret = 3.5*constant + 1.7*high + 0.3*edu - 0.9*high*edu - 0.5*female + 0.8*high*female - 0.5*conf + 0.3*high*conf -0.5*edu*conf+0.8*high*edu*conf+0.15*prior+0.001*prior_squ+0.3*high*prior-0.013*high*prior_squ-0.000037*prior*prior_squ+0.0001392*high*prior*prior_squ + noise_matrix %*% gamma_true_sel + epsilon  
          
        }
        
        else {
          
          constant = rep(1,N)
          epsilon = rnorm(N,sd=4.5)
          high = rbinom(N,1,0.5)
          
          ### Draw noise variables from multivariate distribution
          
          noise_matrix = as.matrix(rmvnorm(n=N,mu = rep(0, n_noise),sigma = Sigma_noise[1:n_noise,1:n_noise]))
          
          female = rbinom(N,1,0.79)
          edu = rbinom(N,1,0.47)
          prior = rexp(n=N, rate=0.1176)
          prior_squ = prior*prior
          conf = rnorm(N,mean=mean(data$conf),sd=1)
          
          
          stockret = 3.5*constant + 1.7*high + 0.3*edu - 0.9*high*edu - 0.5*female + 0.8*high*female - 0.5*conf + 0.3*high*conf -0.5*edu*conf+0.8*high*edu*conf+0.15*prior+0.001*prior_squ+0.3*high*prior-0.013*high*prior_squ-0.000037*prior*prior_squ+0.0001392*high*prior*prior_squ + noise_matrix %*% gamma_true_sel + epsilon  
        }
        
        if (j == 1) {
          
          high_edu = high*edu
          high_female = high*female
          high_conf = high*conf
          edu_conf = edu*conf
          high_edu_conf = high*edu*conf
          high_prior = high*prior
          high_priorsqu = high*prior_squ
          prior_priorsqu = prior*prior_squ
          high_prior_priorsqu = high*prior*prior_squ
          
          dmatrix_simple_train = cbind(data.frame(constant,high,edu,high_edu,female,high_female,conf,high_conf,prior,prior_squ,high_prior),data.frame(noise_matrix))
          dmatrix_complex_train = cbind(data.frame(constant,high,edu,high_edu,female,high_female,conf,high_conf,edu_conf,high_edu_conf,prior,prior_squ,high_prior,high_priorsqu,prior_priorsqu,high_prior_priorsqu),data.frame(noise_matrix))
          stockret_train = stockret
          
          beta_simple = beta_hat(y=stockret_train,design_matrix = dmatrix_simple_train)
          beta_complex = beta_hat(y=stockret_train,design_matrix = dmatrix_complex_train)
          
          
        }
        
        else {
          
          # Generate Y0 and Y1 (real and counterfactual) 
          
          for (high_val in range(0,1)) {
            
            high = rep(high_val,N)  
            
            high_edu = high*edu
            high_female = high*female
            high_conf = high*conf
            edu_conf = edu*conf
            high_edu_conf = high*edu*conf
            high_prior = high*prior
            high_priorsqu = high*prior_squ
            prior_priorsqu = prior*prior_squ
            high_prior_priorsqu = high*prior*prior_squ
            
            dmatrix_simple_test = cbind(data.frame(constant,high,edu,high_edu,female,high_female,conf,high_conf,prior,prior_squ,high_prior),data.frame(noise_matrix))
            dmatrix_complex_test = cbind(data.frame(constant,high,edu,high_edu,female,high_female,conf,high_conf,edu_conf,high_edu_conf,prior,prior_squ,high_prior,high_priorsqu,prior_priorsqu,high_prior_priorsqu),data.frame(noise_matrix))
            dmatrix_true_test = cbind(data.frame(constant,high,edu,high_edu,female,high_female,conf,high_conf,edu_conf,high_edu_conf,prior,prior_squ,high_prior,high_priorsqu,prior_priorsqu,high_prior_priorsqu),data.frame(noise_matrix))
            
            if (high_val == 0) {
              
              Y0_simple = as.matrix(dmatrix_simple_test) %*% beta_simple
              Y0_complex = as.matrix(dmatrix_complex_test) %*% beta_complex
              Y0_true = as.matrix(dmatrix_true_test) %*% beta_true
              
            }
            else {
              
              Y1_simple = as.matrix(dmatrix_simple_test) %*% beta_simple
              Y1_complex = as.matrix(dmatrix_complex_test) %*% beta_complex
              Y1_true = as.matrix(dmatrix_true_test) %*% beta_true
              
            }
          }
        }
      }
      results_simple_Y0 = cbind(results_simple_Y0, Y0_simple)
      results_simple_Y1 = cbind(results_simple_Y1, Y1_simple)
      
      results_complex_Y0 = cbind(results_complex_Y0, Y0_complex)
      results_complex_Y1 = cbind(results_complex_Y1, Y1_complex)
      
      results_true_Y0 = cbind(results_true_Y0, Y0_true)
      results_true_Y1 = cbind(results_true_Y1, Y1_true)
      
      counter = counter + 2
    }
    
    ### Calculate CATE MSE
    
    tau_hat_simple = results_simple_Y1 - results_simple_Y0
    tau_hat_complex = results_complex_Y1 - results_complex_Y0
    tau_true = results_true_Y1 - results_true_Y0
    
    MSE_CATE_simple = colMeans((tau_hat_simple - tau_true)^2)
    MSE_CATE_complex = colMeans((tau_hat_complex - tau_true)^2)
      
    EMSE_CATE_simple = mean(MSE_CATE_simple,na.rm=T)
    EMSE_CATE_complex = mean(MSE_CATE_complex,na.rm=T)
      
    MSE_CATE_simple_sd = sd(MSE_CATE_simple,na.rm=T)
    MSE_CATE_complex_sd = sd(MSE_CATE_complex,na.rm=T)
      
    ### Calculate ATE MSE 
    
    ATE_simple = colMeans(results_simple_Y1) - colMeans(results_simple_Y0)
    ATE_complex = colMeans(results_complex_Y1) - colMeans(results_complex_Y0)
    ATE_true = colMeans(results_true_Y1) - colMeans(results_true_Y0)
    
    MSE_ATE_simple = (ATE_simple - ATE_true)^2
    MSE_ATE_complex = (ATE_complex - ATE_true)^2
     
    EMSE_ATE_simple = mean(MSE_ATE_simple,na.rm=TRUE)  
    EMSE_ATE_complex = mean(MSE_ATE_complex,na.rm=TRUE)
    
    MSE_ATE_simple_sd = sd(MSE_ATE_simple,na.rm=T)
    MSE_ATE_complex_sd = sd(MSE_ATE_complex,na.rm=T)
    
    col = rbind(col,EMSE_CATE_simple,MSE_CATE_simple_sd,EMSE_CATE_complex,MSE_CATE_complex_sd,EMSE_ATE_simple,MSE_ATE_simple_sd,EMSE_ATE_complex,MSE_ATE_complex_sd)
    final_table = cbind(final_table,col)
    
  }
  
  row.names(final_table)[1] = "Noise"  
  return(final_table)
  
}




