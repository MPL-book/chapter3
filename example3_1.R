## Simulate left truncated and right censored survival times

trunc_sim_data3.1 <- function(nsample, beta, 
                              prob_left_0, 
                              right_cens_par1, right_cens_par2){
  
  n <- 0 #set current sample size to 0
  dat <- NULL
  n_all <- 0 #set total number of simulated data sets to 0
  
  while(n < nsample){
    #generate covariates
    X <- cbind(rbinom(1, 1, 0.7), runif(1, 1, 2)) 
    n_all <- n_all + 1
    
    #generate event time
    neg.log.u <- -log(runif(1))
    mu_term <- exp(X%*%beta)
    y_i <- as.numeric((((4*neg.log.u/(2*mu_term)) + 1)^(1/2) - 1)^(1/2))
    
    #generate left truncation time
    L_i <- rbinom(1, 1, 1-prob_left_0) * runif(1, 1, 3)
    
    #generate right truncation time
    R_i <- Inf
    
    #generate independent right censoring time
    cens_i <- runif(1, right_cens_par1, right_cens_par2)
    
    TL <- TR <- censor_type <- NA
    
    if(y_i > cens_i){
      #right censored
      TL <- cens_i
      TR <- Inf
      censor_type <- 0
    }else{
      TL <- y_i
      TR <- y_i
      censor_type <- 1
    }
    
    if(censor_type == 1 & L_i < TL & TR <= R_i){
      #if event time is inside the truncation interval 
      #then individual is in the sample
      dat <- rbind(dat, c(y_i, cens_i, L_i, R_i, X, TL, TR)) 
      n <- n+1 #increase current sample size by 1
    }else if(censor_type == 0 & L_i < TL & TR <= R_i){
      #if right censoring time is inside the truncation interval 
      #then individual is in the sample
      dat <- rbind(dat, c(y_i, cens_i, L_i, R_i, X, TL, TR)) 
      n <- n+1 #increase current sample size by 1
    }
  }
  
  trunc_data <- data.frame(dat)
  colnames(trunc_data) <- c("y_i", "cens_i", "Ltrunc", "Rtrunc", "x1", "X2", "Lcen", "Rcen")
  
  #censoring type
  trunc_data$censor_type <- NA
  trunc_data$censor_type[which(is.infinite
                               (trunc_data$Rcen))] <- 0
  trunc_data$censor_type[which(trunc_data$Lcen == trunc_data$Rcen)] <- 1
  
  #truncation type
  trunc_data$trunc_type = NA
  trunc_data$trunc_type[which(trunc_data$Ltrunc == 0 & trunc_data$Rtrunc == Inf)] <- 1
  trunc_data$trunc_type[which(trunc_data$Ltrunc != 0 & trunc_data$Rtrunc == Inf)] <- 2
  
  out <- list(trunc_data = trunc_data, n_all = n_all)
  return(out)
}

trunc_data_event40 <- trunc_sim_data3.1(nsample = 200,
                                        beta = c(1, -0.5), prob_left_0 = 0, right_cens_par1 = 0.2,
                                        right_cens_par2 = 1.4)
data <- trunc_data_event40$trunc_data

format(head(data), digits = 4)

#check censoring level
table(data$censor_type)

#check truncation types (2 = left truncated)
table(data$trunc_type)

#check truncation proportion
(trunc_data_event40$n_all - 200)/trunc_data_event40$n_all


