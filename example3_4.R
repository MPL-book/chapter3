## A simulation example - generate data

# Write a function to simulate double truncated,
#partly interval censored survival function
trunc_sim_data <- function(nsample, beta, 
                           prob_left_0, prob_right_inf, 
                           right_cens_par1, right_cens_par2,
                           event_length){
  
  n = 0 #set current sample size to 0
  dat = NULL
  
  while(n < nsample){
    #generate covariates
    X <- cbind(rbinom(1, 1, 0.4), rnorm(1, 0, 1)) 
    
    #generate event time
    neg.log.u <- -log(runif(1))
    mu_term <- exp(X%*%beta)
    y_i <- as.numeric((neg.log.u/(mu_term))^(1/3))
    
    #generate left truncation time
    L_i <- rbinom(1, 1, 1-prob_left_0) * runif(1, 0.1, 0.5) 
    
    #generate right truncation time
    R_i <- (Inf^rbinom(1, 1, prob_right_inf)) + L_i + runif(1, 0.5, 3) 
    
    #generate event observation times
    obs_n_i <- rpois(1, 40) + 1
    int_obs_i <- runif(obs_n_i, 0.2, 0.4)
    t_obs_i <- cumsum(int_obs_i)
    
    #generate independent dropout time
    cens_i <- runif(1, right_cens_par1, right_cens_par2)
    
    #truncate observation times at dropout
    t_obs_i_trunc <- c(t_obs_i[which(t_obs_i < cens_i)], cens_i)
    TL <- TR <- censor_type <- NA
    
    if(y_i > max(t_obs_i_trunc)){
      #right censored
      TL <- max(t_obs_i_trunc)
      TR <- Inf
      censor_type <- 0
    }else if(y_i < min(t_obs_i_trunc)){
      #left censored
      TL <- 0
      TR <- min(t_obs_i_trunc)
      censor_type <- 2
    }else{
      tiL_temp <- t_obs_i_trunc[max(which(t_obs_i_trunc < y_i))]
      tiR_temp <- t_obs_i_trunc[min(which(t_obs_i_trunc > y_i))]
      
      if((tiR_temp - tiL_temp) < event_length){
        #event time (if interval is small enough)
        TL <- y_i
        TR <- y_i
        censor_type <- 1
      }else{
        #interval censored
        TL <- tiL_temp
        TR <- tiR_temp
        censor_type <- 3
      }
    }
    
    if(censor_type == 1 & L_i < TL & TR < R_i){
      #if event time is inside the truncation interval 
      #then individual is in the sample
      dat <- rbind(dat, c(y_i, L_i, R_i, X, TL, TR)) 
      n <- n+1 #increase current sample size by 1
    }else if(censor_type == 0 & L_i < TL & TR < R_i){
      #if right censoring time is inside the truncation interval 
      #then individual is in the sample
      dat <- rbind(dat, c(y_i, L_i, R_i, X, TL, TR)) 
      n <- n+1 #increase current sample size by 1
    }else if(censor_type == 2 & L_i < TL & TR < R_i){
      #if left censoring time is inside the truncation interval 
      #then individual is in the sample
      dat <- rbind(dat, c(y_i, L_i, R_i, X, TL, TR))
      n <- n+1 #increase current sample size by 1
    }else if(censor_type == 3 & L_i < TL & TR < R_i){
      #if censoring interval is inside the truncation interval 
      #then individual is in the sample
      dat <- rbind(dat, c(y_i, L_i, R_i, X, TL, TR))
      n <- n+1 #increase current sample size by 1
    }
  }
  
  trunc_data <- data.frame(dat)
  colnames(trunc_data) <- c("y_i", "trunc_time_L", "trunc_time_R", "x", "X2", "censor_time_L", "censor_time_R")
  
  #censoring type
  trunc_data$censor_type <- NA
  trunc_data$censor_type[which(is.infinite(trunc_data$censor_time_R))] <- 0
  trunc_data$censor_type[which(trunc_data$censor_time_L == 0)] <- 2
  trunc_data$censor_type[which(trunc_data$censor_time_L == trunc_data$censor_time_R)] <- 1
  trunc_data$censor_type[is.na(trunc_data$censor_type)] <- 3
  
  #truncation type
  trunc_data$trunc_type = NA
  trunc_data$trunc_type[which(trunc_data$trunc_time_L == 0 & trunc_data$trunc_time_R == Inf)] <- 1
  trunc_data$trunc_type[which(trunc_data$trunc_time_L != 0 & trunc_data$trunc_time_R == Inf)] <- 2
  trunc_data$trunc_type[which(trunc_data$trunc_time_L == 0 & trunc_data$trunc_time_R != Inf)] <- 3
  trunc_data$trunc_type[which(trunc_data$trunc_time_L != 0 & trunc_data$trunc_time_R != Inf)] <- 4
  
  return(trunc_data)
}

#generate double truncated data with no censoring
trunc_data_event100  <- trunc_sim_data(nsample = 200,
                                       beta = c(0.5, -0.5), 
                                       prob_left_0 = 0, prob_right_inf = 0, 
                                       right_cens_par1 = 10, right_cens_par2 = 100, event_length = 100)


#generate double truncated data with only interval censoring
trunc_data_interval  = trunc_sim_data(nsample = 200,
                                      beta = c(0.5, -0.5), 
                                      prob_left_0 = 0, prob_right_inf = 0, 
                                      right_cens_par1 = 10, right_cens_par2 = 100, event_length = 0.2)




