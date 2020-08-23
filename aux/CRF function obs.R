### Start with causal forest 

crf_sim_obs <- function(n_obs=c(100), m_runs=10,dgp_complex=T) {
  
  n_sim = length(n_obs)
  
  ### Set up results data_frame
  
  final_table = c()
  
  ### Run simulation study for different n 
  
  for (k in seq(1,n_sim,1)) {
    
    n = n_obs[k]
    
    print(n)
    col = c()
    col = cbind(col,n)
    
    #beta_true = matrix(data=c(3.5,1.7,0.3,-0.9,-0.5, 0.8,-0.5,0.3,-0.5,0.8,0.15,0.3,-0.013,-0.000037,0.0001392), nrow=15)
    beta_true = matrix(data=c(3.5,1.7,0,-0.9,-0.5, 0.8,-0.5,0,-0.5,0.8,0.15,0.3,-0.013,-0.000037,0.0001392), nrow=15)
    
    set.seed(123)
    N = n 
    M = m_runs
    
    counter = 123
    
    
    for (i in seq(1:M)) {
      
      set.seed(counter)
      
      print(i)
      
      # Set up results dataframes
      
      if (i == 1) {
        
        results_tau_crf = data.frame(x=rep(NA,N))
        
        results_true_Y0 = data.frame(x=rep(NA,N))
        results_true_Y1 = data.frame(x=rep(NA,N))
        results_ATE_crf = data.frame(x=rep(NA,1))
        
      }
      
      # Draw datasample (test and training)
      
      for (j in seq(1:2)) {
        
        set.seed(counter+j)
        
        if (dgp_complex == TRUE) {
          
          constant = rep(1,N)
          epsilon = rnorm(N,sd=4.5)
          
          high = rbinom(N,1,0.5)
          prior = rgamma(n=N, scale=9.6812789, shape=0.8885924)
          prior_squ = prior*prior
          
          ##### Draw female, education and conf from a multivariate distribution ##### 
          
          ### First draw female and education (multivariate binomial distribution, correlation structure taken from data)
          
          commprob = matrix(c(sum(data$atleast_bachelor)/dim(data)[1],dim(data[data$female == 1 & data$atleast_bachelor == 1,])[1]/dim(data)[1]
                              ,dim(data[data$female == 1 & data$atleast_bachelor == 1,])[1]/dim(data)[1],sum(data$female)/dim(data)[1]),nrow=2,ncol=2)
          mprob = c(sum(data$female)/dim(data)[1],sum(data$atleast_bachelor)/dim(data)[1])
          
          sample_edu_fem = data.frame(rmvbin(n=N,margprob=as.matrix(mprob),commonprob=commprob))
          names(sample_edu_fem) = c("female","atleast_bachelor")
          
          ### Let conf depend on gender (correlation structure taken from data)
          
          cat_vector = as.numeric(names(table(data[data$female == 0,]$z_conf_prior)))
          
          prob_male = c(data.frame(table(data[data$female == 0,]$z_conf_prior)/length(data[data$female == 0,]$z_conf_prior))$Freq)
          prob_female = c(data.frame(table(data[data$female == 1,]$z_conf_prior)/length(data[data$female == 1,]$z_conf_prior))$Freq)
          
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
          edu = sample_edu_fem_conf$atleast_bachelor
          conf = sample_edu_fem_conf$conf
          
          #stockret = 3.5*constant + 1.7*high + 0.3*edu - 0.9*high*edu - 0.5*female + 0.8*high*female - 0.5*conf + 0.3*high*conf -0.5*edu*conf+0.8*high*edu*conf+0.15*prior+0.3*high*prior-0.013*high*prior_squ-0.000037*prior*prior_squ+0.0001392*high*prior*prior_squ + epsilon  
          stockret = 3.5*constant + 1.7*high + 0*edu - 0.9*high*edu - 0.5*female + 0.8*high*female - 0.5*conf + 0*high*conf -0.5*edu*conf+0.8*high*edu*conf+0.15*prior+0.3*high*prior-0.013*high*prior_squ-0.000037*prior*prior_squ+0.0001392*high*prior*prior_squ + epsilon  
          
        }
        
        else {
          
          constant = rep(1,N)
          epsilon = rnorm(N,sd=4.5)
          high = rbinom(N,1,0.5)
          female = rbinom(N,1,0.79)
          edu = rbinom(N,1,0.47)
          prior = rgamma(n=N, scale=9.6812789, shape=0.8885924)
          prior_squ = prior*prior
          conf = rnorm(N,mean=mean(data$z_conf_prior),sd=1)
          
          #stockret = 3.5*constant + 1.7*high + 0.3*edu - 0.9*high*edu - 0.5*female + 0.8*high*female - 0.5*conf + 0.3*high*conf -0.5*edu*conf+0.8*high*edu*conf+0.15*prior+0.3*high*prior-0.013*high*prior_squ-0.000037*prior*prior_squ+0.0001392*high*prior*prior_squ + epsilon  
          stockret = 3.5*constant + 1.7*high + 0*edu - 0.9*high*edu - 0.5*female + 0.8*high*female - 0.5*conf + 0*high*conf -0.5*edu*conf+0.8*high*edu*conf+0.15*prior+0.3*high*prior-0.013*high*prior_squ-0.000037*prior*prior_squ+0.0001392*high*prior*prior_squ + epsilon  
          
        }
        
        if (j == 1) {
          
          # Set up training sample for crf
          
          covariates_train = data.frame(edu,female,conf,prior,prior_squ)
          high_train = high
          stockret_train = stockret
          
          # Fit CRF with training data
          
          cf_fitting <- causal_forest(
            
            X <- as.matrix(covariates_train),
            Y <- as.matrix(stockret_train),
            W <- as.matrix(high_train),
            
            
            seed <- 123,
            tune.parameters = "all",
            num.trees = 2000
            
          )
        }
        
        else {
          
          # Set up test sample for crf 
          
          covariates_test = data.frame(edu,female,conf,prior,prior_squ)
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
            
            dmatrix_true_test = data.frame(constant,high,edu,high_edu,female,high_female,conf,high_conf,edu_conf,high_edu_conf,prior,high_prior,high_priorsqu,prior_priorsqu,high_prior_priorsqu)
            
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
    
    ### Calculate average of MSE and ATE
    
    results_tau_true = results_true_Y1 - results_true_Y0
    
    MSE_crf = colMeans((results_tau_crf - results_tau_true)^2)
    MSE_crf_avg = mean(MSE_crf,na.rm = T)
    MSE_crf_sd = sd(MSE_crf,na.rm = T)
    
    ATE_crf_avg = rowMeans(results_ATE_crf,na.rm =TRUE)
    ATE_crf_sd = sd(results_ATE_crf, na.rm=TRUE)
    
    ATE_true = colMeans(results_true_Y1) - colMeans(results_true_Y0)
    ATE_true_avg = mean(ATE_true,na.rm=T)
    ATE_true_sd = sd(ATE_true, na.rm=T)
    
    col = rbind(col,MSE_crf_avg,MSE_crf_sd,ATE_crf_avg,ATE_crf_sd,ATE_true_avg,ATE_true_sd)
    
    final_table = cbind(final_table,col)
    
  }
  
  return(final_table)
  
}


### OLS model

ols_sim_obs <- function(n_obs=c(100), m_runs=10,dgp_complex=T) {
  
  n_sim = length(n_obs)
  
  ### Set up results data_frame
  
  final_table = c()
  
  ### Run simulation study for different n 
  
  for (k in seq(1,n_sim,1)) {
    
    n = as.integer(n_obs[k])
    
    print(n)
    col = c()
    col = cbind(col,n)
    
    #beta_true = matrix(data=c(3.5,1.7,0.3,-0.9,-0.5, 0.8,-0.5,0.3,-0.5,0.8,0.15,0.3,-0.013,-0.000037,0.0001392), nrow=15)
    beta_true = matrix(data=c(3.5,1.7,0,-0.9,-0.5, 0.8,-0.5,0,-0.5,0.8,0.15,0.3,-0.013,-0.000037,0.0001392), nrow=15)
    
    
    set.seed(123)
    N = n 
    M = m_runs
    
    counter = 123
    
    
    for (i in seq(1:M)) {
      
      set.seed(counter)
      
      print(i)
      
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
          prior = rgamma(n=N, scale=9.6812789, shape=0.8885924)
          prior_squ = prior*prior
          
          ##### Draw female, education and conf from a multivariate distribution ##### 
          
          ### First draw female and education (multivariate binomial distribution, correlation structure taken from data)
          
          commprob = matrix(c(sum(data$atleast_bachelor)/dim(data)[1],dim(data[data$female == 1 & data$atleast_bachelor == 1,])[1]/dim(data)[1]
                              ,dim(data[data$female == 1 & data$atleast_bachelor == 1,])[1]/dim(data)[1],sum(data$female)/dim(data)[1]),nrow=2,ncol=2)
          mprob = c(sum(data$female)/dim(data)[1],sum(data$atleast_bachelor)/dim(data)[1])
          
          sample_edu_fem = data.frame(rmvbin(n=N,margprob=as.matrix(mprob),commonprob=commprob))
          names(sample_edu_fem) = c("female","atleast_bachelor")
          
          ### Let conf depend on gender (correlation structure taken from data)
          
          cat_vector = as.numeric(names(table(data[data$female == 0,]$z_conf_prior)))
          
          prob_male = c(data.frame(table(data[data$female == 0,]$z_conf_prior)/length(data[data$female == 0,]$z_conf_prior))$Freq)
          prob_female = c(data.frame(table(data[data$female == 1,]$z_conf_prior)/length(data[data$female == 1,]$z_conf_prior))$Freq)
          
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
          edu = sample_edu_fem_conf$atleast_bachelor
          conf = sample_edu_fem_conf$conf
          
          #stockret = 3.5*constant + 1.7*high + 0.3*edu - 0.9*high*edu - 0.5*female + 0.8*high*female - 0.5*conf + 0.3*high*conf -0.5*edu*conf+0.8*high*edu*conf+0.15*prior+0.3*high*prior-0.013*high*prior_squ-0.000037*prior*prior_squ+0.0001392*high*prior*prior_squ + epsilon  
          stockret = 3.5*constant + 1.7*high + 0*edu - 0.9*high*edu - 0.5*female + 0.8*high*female - 0.5*conf + 0*high*conf -0.5*edu*conf+0.8*high*edu*conf+0.15*prior+0.3*high*prior-0.013*high*prior_squ-0.000037*prior*prior_squ+0.0001392*high*prior*prior_squ + epsilon  
          
        }
        
        else {
          
          constant = rep(1,N)
          epsilon = rnorm(N,sd=4.5)
          high = rbinom(N,1,0.5)
          female = rbinom(N,1,0.79)
          edu = rbinom(N,1,0.47)
          prior = rgamma(n=N, scale=9.6812789, shape=0.8885924)
          prior_squ = prior*prior
          conf = rnorm(N,mean=mean(data$z_conf_prior),sd=1)
          
          #stockret = 3.5*constant + 1.7*high + 0.3*edu - 0.9*high*edu - 0.5*female + 0.8*high*female - 0.5*conf + 0.3*high*conf -0.5*edu*conf+0.8*high*edu*conf+0.15*prior+0.3*high*prior-0.013*high*prior_squ-0.000037*prior*prior_squ+0.0001392*high*prior*prior_squ + epsilon  
          stockret = 3.5*constant + 1.7*high + 0*edu - 0.9*high*edu - 0.5*female + 0.8*high*female - 0.5*conf + 0*high*conf -0.5*edu*conf+0.8*high*edu*conf+0.15*prior+0.3*high*prior-0.013*high*prior_squ-0.000037*prior*prior_squ+0.0001392*high*prior*prior_squ + epsilon  
          
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
          
          dmatrix_simple_train = data.frame(constant,high,edu,high_edu,female,high_female,conf,high_conf,prior,prior_squ,high_prior)
          dmatrix_complex_train = data.frame(constant,high,edu,high_edu,female,high_female,conf,high_conf,edu_conf,high_edu_conf,prior,prior_squ,high_prior,high_priorsqu,prior_priorsqu,high_prior_priorsqu)
          stockret_train = stockret
          
          beta_simple = beta_hat(y=stockret_train,design_matrix = dmatrix_simple_train)
          beta_complex = beta_hat(y=stockret_train,design_matrix = dmatrix_complex_train)
          beta_true = matrix(data=c(3.5,1.7,0.3,-0.9,-0.5, 0.8,-0.5,0.3,-0.5,0.8,0.15,0.3,-0.013,-0.000037,0.0001392), nrow=15)
          
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
            
            dmatrix_simple_test = data.frame(constant,high,edu,high_edu,female,high_female,conf,high_conf,prior,prior_squ,high_prior)
            dmatrix_complex_test = data.frame(constant,high,edu,high_edu,female,high_female,conf,high_conf,edu_conf,high_edu_conf,prior,prior_squ,high_prior,high_priorsqu,prior_priorsqu,high_prior_priorsqu)
            dmatrix_true_test = data.frame(constant,high,edu,high_edu,female,high_female,conf,high_conf,edu_conf,high_edu_conf,prior,high_prior,high_priorsqu,prior_priorsqu,high_prior_priorsqu)
            
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
    
    ### Calculate average of MSE and ATE
    
    tau_hat_simple = results_simple_Y1 - results_simple_Y0
    tau_hat_complex = results_complex_Y1 - results_complex_Y0
    tau_true = results_true_Y1 - results_true_Y0
    
    MSE_simple = colMeans((tau_hat_simple - tau_true)^2)
    MSE_complex = colMeans((tau_hat_complex - tau_true)^2)
    
    ATE_simple = colMeans(results_simple_Y1) - colMeans(results_simple_Y0)
    ATE_complex = colMeans(results_complex_Y1) - colMeans(results_complex_Y0)
    ATE_true = colMeans(results_true_Y1) - colMeans(results_true_Y0)
    
    MSE_simple_avg = mean(MSE_simple,na.rm=T)
    MSE_complex_avg = mean(MSE_complex,na.rm=T)
    
    MSE_simple_sd = sd(MSE_simple,na.rm=T)
    MSE_complex_sd = sd(MSE_complex,na.rm=T)
    
    ATE_simple_avg = mean(ATE_simple,na.rm=T)
    ATE_complex_avg = mean(ATE_complex,na.rm=T)
    ATE_true_avg = mean(ATE_true,na.rm=T)
    
    ATE_simple_sd = sd(ATE_simple,na.rm=T)
    ATE_complex_sd = sd(ATE_complex,na.rm=T)
    ATE_true_sd = sd(ATE_true,na.rm=T)
    
    col = rbind(col,MSE_simple_avg,MSE_simple_sd,MSE_complex_avg,MSE_complex_sd,ATE_simple_avg,ATE_simple_sd,ATE_complex_avg,ATE_complex_sd,ATE_true_avg,ATE_true_sd )
    final_table = cbind(final_table,col)
    
  }
  
  return(final_table)
  
}
