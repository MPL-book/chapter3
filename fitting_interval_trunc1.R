
CoxMPL_Trunc = function(tL_cens, tR_cens, cov, tL_trunc, tR_trunc, max.iter = c(5, 1000), n.knots.h0 = 4 , step_size=c(1,1)){
  
  kappa = 1/0.6
  
  #censoring type indicators
  
  cens_type = rep(NA, length(tL_cens))
  cens_type[which(tL_cens == tR_cens)] = 1
  cens_type[which(tR_cens == Inf)] = 0
  cens_type[which(tL_cens == 0)] = 2
  cens_type[which(tL_cens != tR_cens)] = 3
  
  cens_type_ind_Right = sum(cens_type == 0)
  cens_type_ind_Right = as.numeric(cens_type_ind_Right > 0)
  
  #truncation type indicators
  trunc_type = rep(NA, length(tL_cens))
  trunc_type[which(tL_trunc == 0 & tR_trunc == Inf)] = 1
  trunc_type[which(tL_trunc > 0 & tR_trunc == Inf)] = 2
  trunc_type[which(tL_trunc == 0 & tR_trunc < Inf)] = 3
  trunc_type[which(tL_trunc > 0 & tR_trunc < Inf)] = 4
  
  #choose knots for baseline hazard function
  times = c(tL_cens, tR_cens[which(!is.infinite(tR_cens))], tL_trunc, tR_trunc[which(!is.infinite(tR_trunc))])
  
  cen_times = c(tL_cens, tR_cens)
  
  #note in position of knots, we are placing a knot at t=34 (the end of the left truncation window)
  int.knots = c(34, quantile(cen_times, seq(0.3, 0.9, length.out=n.knots.h0)))
  bound.knots = c(-1e-3, max(times) + 1e-3)
  
  psi_e = mSpline(tL_cens, knots = int.knots, Boundary.knots = bound.knots, degree = 3, intercept = FALSE)[which(cens_type == 1),]
  Psi_e = mSpline(tL_cens, knots = int.knots, Boundary.knots = bound.knots, degree = 3, intercept = FALSE, integral = TRUE)[which(cens_type == 1),]
  
  Psi_c_r = mSpline(tL_cens, knots = int.knots, Boundary.knots = bound.knots, degree = 3, intercept = FALSE, integral = TRUE)[which(cens_type == 0),]
  Psi_c_l = mSpline(tR_cens, knots = int.knots, Boundary.knots = bound.knots, degree = 3, intercept = FALSE, integral = TRUE)[which(cens_type == 2),]
  
  Psi_c_i1 = mSpline(tL_cens, knots = int.knots, Boundary.knots = bound.knots, degree = 3, intercept = FALSE, integral = TRUE)[which(cens_type == 3),]
  Psi_c_i2 = mSpline(tR_cens, knots = int.knots, Boundary.knots = bound.knots, degree = 3, intercept = FALSE, integral = TRUE)[which(cens_type == 3),]
  
  lambda = 1e-4
  
  #initial quantities
  n = length(tL_cens)
  p = ncol(cov)
  m = ncol(psi_e)
  df.theta = n
  
  R = theta_penalty_f(3, int.knots, bound.knots)
  
  beta = matrix(0, ncol = 1, nrow = p)
  theta = matrix(0.5, ncol = 1, nrow = m)
  
  ##get quantities for likelihood function
  
  #compute S(L) and S(R) for everybody
  
  Psi_trunc_L_left = mSpline(tL_trunc, knots = int.knots, Boundary.knots = bound.knots, degree = 3, intercept = FALSE, integral = TRUE)[which(trunc_type == 2),]
  Psi_trunc_R_left = matrix(1, ncol = ncol(Psi_trunc_L_left), nrow = nrow(Psi_trunc_L_left))
  x_trunc_left = cov[which(trunc_type == 2),]
  
  Psi_trunc_L_interval = mSpline(tL_trunc, knots = int.knots, Boundary.knots = bound.knots, degree = 3, intercept = FALSE, integral = TRUE)[which(trunc_type == 4),]
  Psi_trunc_R_interval = mSpline(tR_trunc, knots = int.knots, Boundary.knots = bound.knots, degree = 3, intercept = FALSE, integral = TRUE)[which(trunc_type == 4),]
  x_trunc_interval = cov[which(trunc_type == 4),]
  
  Psi_trunc_L_right = mSpline(tL_trunc, knots = int.knots, Boundary.knots = bound.knots, degree = 3, intercept = FALSE, integral = TRUE)[which(trunc_type == 3),]
  Psi_trunc_R_right = mSpline(tR_trunc, knots = int.knots, Boundary.knots = bound.knots, degree = 3, intercept = FALSE, integral = TRUE)[which(trunc_type == 3),]
  x_trunc_right = cov[which(trunc_type == 3),]
  
  
  for(iter in 1:max.iter[1]){
    for(it in 1:max.iter[2]){
      
      S_trunc_L_left = exp(-(Psi_trunc_L_left %*% theta) * exp(x_trunc_left %*% beta))
      S_trunc_R_left = exp(-(Psi_trunc_R_left %*% theta) * exp(x_trunc_left %*% beta))
      
      S_trunc_L_right = exp(-(Psi_trunc_L_right %*% theta) * exp(x_trunc_right %*% beta))
      S_trunc_R_right = exp(-(Psi_trunc_R_right %*% theta) * exp(x_trunc_right %*% beta))
      
      S_trunc_L_interval = exp(-(Psi_trunc_L_interval %*% theta) * exp(x_trunc_interval %*% beta))
      S_trunc_R_interval = exp(-(Psi_trunc_R_interval %*% theta) * exp(x_trunc_interval %*% beta))
      
      log_S_trunc_diff = sum(log(S_trunc_L_left + 1e-6)) + sum(log(S_trunc_L_right - S_trunc_R_right + 1e-5)) + 
        sum(log(S_trunc_L_interval - S_trunc_R_interval + 1e-6))
      
      #likelihood for event times
      
      x_e = cov[which(cens_type == 1),]
      h0t_e = psi_e %*% theta
      H0t_e = Psi_e %*% theta
      Ht_e = H0t_e * exp(x_e %*% beta)
      St_e = exp(-Ht_e)
      
      #likelihood for right censored times
      x_rightC = cov[which(cens_type == 0 ),]
      S_cens_right = exp(-( Psi_c_r %*% theta) * exp(x_rightC %*% beta))
      
      #likelihood for left censored times
      x_c_left = cov[which(cens_type == 2),]
      S_cens_left = exp(-( Psi_c_l %*% theta) * exp(x_c_left %*% beta))
      
      #likelihood for interval censored times
      x_c_interval = cov[which(cens_type == 3),]
      S_cens_i1 = exp(-( Psi_c_i1 %*% theta) * exp(x_c_interval %*% beta))
      S_cens_i2 = exp(-( Psi_c_i2 %*% theta) * exp(x_c_interval %*% beta))
      
      log_lik = sum(log(h0t_e) + x_e %*% beta - Ht_e) + sum(log(S_cens_right)) + 
        sum(log(1-S_cens_left)) + sum(log(S_cens_i1 - S_cens_i2 +1e-6)) - 
        log_S_trunc_diff - lambda*t(theta) %*% R %*% theta
      
      #### estimate beta
      #truncation hazard functions
      H_trunc_L_left = -log(S_trunc_L_left)
      H_trunc_R_left = -log(S_trunc_R_left)
      
      H_trunc_L_right = -log(S_trunc_L_right)
      H_trunc_R_right = -log(S_trunc_R_right + 1e-6)
      
      H_trunc_L_interval = -log(S_trunc_L_interval)
      H_trunc_R_interval = -log(S_trunc_R_interval)
      
      #censoring hazard functions
      H_cens_right = -log(S_cens_right)
      H_cens_left = -log(S_cens_left)
      H_cens_i1 = -log(S_cens_i1)
      H_cens_i2 = -log(S_cens_i2)
      
      beta_score = t(x_e) %*% (c(1 - Ht_e)) -
        t(x_rightC) %*% H_cens_right + 
        t(x_c_left) %*% c(H_cens_left * S_cens_left/(1-S_cens_left)) -
        t(x_c_interval) %*% c((S_cens_i1*H_cens_i1 - S_cens_i2*H_cens_i2)/(S_cens_i1-S_cens_i2)) +
        t(x_trunc_left) %*% c((S_trunc_L_left*H_trunc_L_left)/(S_trunc_L_left)) + 
        t(x_trunc_right) %*% c((S_trunc_L_right*H_trunc_L_right - S_trunc_R_right*H_trunc_R_right)/(S_trunc_L_right-S_trunc_R_right)) + 
        t(x_trunc_interval) %*% c((S_trunc_L_interval*H_trunc_L_interval - S_trunc_R_interval*H_trunc_R_interval)/(S_trunc_L_interval-S_trunc_R_interval))
      
      beta_hess = -t(x_e) %*% diag(c(Ht_e)) %*% x_e - 
        t(x_rightC) %*% diag(c(H_cens_right)) %*% x_rightC + 
        t(x_c_left) %*% diag(c((-  S_cens_left*H_cens_left^2 - (S_cens_left^2)*H_cens_left)/(1-S_cens_left)^2)) %*% x_c_left - 
        t(x_c_interval) %*% diag(c(S_cens_i1*H_cens_i1/(S_cens_i1-S_cens_i2))) %*% x_c_interval -
        t(x_c_interval) %*% diag(c(S_cens_i1*S_cens_i2*(H_cens_i1-H_cens_i1)^2/(S_cens_i1-S_cens_i2)^2)) %*% x_c_interval + 
        t(x_trunc_left) %*% diag(c(S_trunc_L_left*H_trunc_L_left/(S_trunc_L_left-0))) %*% x_trunc_left + 
        t(x_trunc_right) %*% diag(c(S_trunc_L_right*H_trunc_L_right/(S_trunc_L_right-S_trunc_R_right))) %*% x_trunc_right -
        t(x_trunc_right) %*% diag(c(S_trunc_R_right*H_trunc_R_right/(S_trunc_L_right-S_trunc_R_right))) %*% x_trunc_right +
        t(x_trunc_right) %*% diag(c(S_trunc_L_right*S_trunc_R_right*(H_trunc_L_right-H_trunc_R_right)^2/(S_trunc_L_right-S_trunc_R_right)^2)) %*% x_trunc_right + 
        t(x_trunc_interval) %*% diag(c(S_trunc_L_interval*H_trunc_L_interval/(S_trunc_L_interval-S_trunc_R_interval))) %*% x_trunc_interval -
        t(x_trunc_interval) %*% diag(c(S_trunc_R_interval*H_trunc_R_interval/(S_trunc_L_interval-S_trunc_R_interval))) %*% x_trunc_interval +
        t(x_trunc_interval) %*% diag(c(S_trunc_L_interval*S_trunc_R_interval*(H_trunc_L_interval-H_trunc_R_interval)^2/(S_trunc_L_interval-S_trunc_R_interval)^2)) %*% x_trunc_interval
        
      
      beta_hess_neg = - beta_hess
      
      beta_old = beta
      
      beta = beta_old + solve(beta_hess_neg) %*% (beta_score)
      
      
      ##re calculate log likelihood
      
      S_trunc_L_left = exp(-(Psi_trunc_L_left %*% theta) * exp(x_trunc_left %*% beta))
      S_trunc_R_left = exp(-(Psi_trunc_R_left %*% theta) * exp(x_trunc_left %*% beta))
      
      S_trunc_L_right = exp(-(Psi_trunc_L_right %*% theta) * exp(x_trunc_right %*% beta))
      S_trunc_R_right = exp(-(Psi_trunc_R_right %*% theta) * exp(x_trunc_right %*% beta))
      
      S_trunc_L_interval = exp(-(Psi_trunc_L_interval %*% theta) * exp(x_trunc_interval %*% beta))
      S_trunc_R_interval = exp(-(Psi_trunc_R_interval %*% theta) * exp(x_trunc_interval %*% beta))
      
      log_S_trunc_diff = sum(log(S_trunc_L_left + 1e-6)) + sum(log(S_trunc_L_right - S_trunc_R_right + 1e-5)) + 
        sum(log(S_trunc_L_interval - S_trunc_R_interval + 1e-6))
      
      #likelihood for event times
      
      x_e = cov[which(cens_type == 1),]
      h0t_e = psi_e %*% theta
      H0t_e = Psi_e %*% theta
      Ht_e = H0t_e * exp(x_e %*% beta)
      St_e = exp(-Ht_e)
      
      #likelihood for right censored times
      x_rightC = cov[which(cens_type == 0 ),]
      S_cens_right = exp(-( Psi_c_r %*% theta) * exp(x_rightC %*% beta))
      
      #likelihood for left censored times
      x_c_left = cov[which(cens_type == 2),]
      S_cens_left = exp(-( Psi_c_l %*% theta) * exp(x_c_left %*% beta))
      
      #likelihood for interval censored times
      x_c_interval = cov[which(cens_type == 3),]
      S_cens_i1 = exp(-( Psi_c_i1 %*% theta) * exp(x_c_interval %*% beta))
      S_cens_i2 = exp(-( Psi_c_i2 %*% theta) * exp(x_c_interval %*% beta))
      
      log_lik_old = log_lik
      log_lik = sum(log(h0t_e) + x_e %*% beta - Ht_e) + sum(log(S_cens_right)) + 
        sum(log(1-S_cens_left)) + sum(log(S_cens_i1 - S_cens_i2 + 1e-5)) - 
        log_S_trunc_diff - lambda*t(theta) %*% R %*% theta
      
      if((log_lik < log_lik_old) & (step_size[1] == 1)){
        
        i = 0
        omega1 = 1/kappa
        
        
        while(log_lik<log_lik_old){
          beta = beta_old + omega1 * solve(beta_hess_neg) %*% (beta_score)
          
          #update log-likelihood
          
          S_trunc_L_left = exp(-(Psi_trunc_L_left %*% theta) * exp(x_trunc_left %*% beta))
          S_trunc_R_left = exp(-(Psi_trunc_R_left %*% theta) * exp(x_trunc_left %*% beta))
          
          S_trunc_L_right = exp(-(Psi_trunc_L_right %*% theta) * exp(x_trunc_right %*% beta))
          S_trunc_R_right = exp(-(Psi_trunc_R_right %*% theta) * exp(x_trunc_right %*% beta))
          
          S_trunc_L_interval = exp(-(Psi_trunc_L_interval %*% theta) * exp(x_trunc_interval %*% beta))
          S_trunc_R_interval = exp(-(Psi_trunc_R_interval %*% theta) * exp(x_trunc_interval %*% beta))
          
          log_S_trunc_diff = sum(log(S_trunc_L_left + 1e-6)) + sum(log(S_trunc_L_right - S_trunc_R_right + 1e-5)) + 
            sum(log(S_trunc_L_interval - S_trunc_R_interval + 1e-6))
          
          #likelihood for event times
          
          x_e = cov[which(cens_type == 1),]
          h0t_e = psi_e %*% theta
          H0t_e = Psi_e %*% theta
          Ht_e = H0t_e * exp(x_e %*% beta)
          St_e = exp(-Ht_e)
          
          #likelihood for right censored times
          x_rightC = cov[which(cens_type == 0 ),]
          S_cens_right = exp(-( Psi_c_r %*% theta) * exp(x_rightC %*% beta))
          
          #likelihood for left censored times
          x_c_left = cov[which(cens_type == 2),]
          S_cens_left = exp(-( Psi_c_l %*% theta) * exp(x_c_left %*% beta))
          
          #likelihood for interval censored times
          x_c_interval = cov[which(cens_type == 3),]
          S_cens_i1 = exp(-( Psi_c_i1 %*% theta) * exp(x_c_interval %*% beta))
          S_cens_i2 = exp(-( Psi_c_i2 %*% theta) * exp(x_c_interval %*% beta))
          
          log_lik = sum(log(h0t_e) + x_e %*% beta - Ht_e) + sum(log(S_cens_right)) + 
            sum(log(1-S_cens_left)) + sum(log(S_cens_i1 - S_cens_i2 + 1e-5)) - 
            log_S_trunc_diff - lambda*t(theta) %*% R %*% theta
          
          #update value of omega1
          if(omega1>=1e-2){
            omega1 = omega1/kappa
          }else if(omega1<1e-2 & omega1>=1e-5){
            omega1 = omega1*(5e-2)
          }else if(omega1<1e-5){
            omega1 = omega1*(1e-5)
          }
          i = i+1
          if(i>1000){break}
          
          
        }
      }
      
      
      ##estimate theta
      S_trunc_L_left = exp(-(Psi_trunc_L_left %*% theta) * exp(x_trunc_left %*% beta))
      
      
      TwoLRtheta = as.numeric(2*lambda)*(R%*%theta)
      
      theta_score = t(1/h0t_e) %*% psi_e - t(exp(x_e %*% beta)) %*% Psi_e - 
        t(exp(x_rightC %*% beta)) %*% Psi_c_r +
        t(exp(x_c_left %*% beta)*S_cens_left/(1-S_cens_left)) %*% Psi_c_l - 
        t(exp(x_c_interval %*% beta)*S_cens_i1/(S_cens_i1-S_cens_i2)) %*% Psi_c_i1 + 
        t(exp(x_c_interval %*% beta)*S_cens_i2/(S_cens_i1-S_cens_i2)) %*% Psi_c_i2 +
        t(exp(x_trunc_left %*% beta)*S_trunc_L_left/(S_trunc_L_left-S_trunc_R_left)) %*% Psi_trunc_L_left -
        t(exp(x_trunc_left %*% beta)*S_trunc_R_left/(S_trunc_L_left-S_trunc_R_left)) %*% Psi_trunc_R_left + 
        t(exp(x_trunc_right %*% beta)*S_trunc_L_right/(S_trunc_L_right-S_trunc_R_right)) %*% Psi_trunc_L_right -
        t(exp(x_trunc_right %*% beta)*S_trunc_R_right/(S_trunc_L_right-S_trunc_R_right)) %*% Psi_trunc_R_right +
        t(exp(x_trunc_interval %*% beta)*S_trunc_L_interval/(S_trunc_L_interval-S_trunc_R_interval)) %*% Psi_trunc_L_interval -
        t(exp(x_trunc_interval %*% beta)*S_trunc_R_interval/(S_trunc_L_interval-S_trunc_R_interval)) %*% Psi_trunc_R_interval
        
      theta_score_neg = t(exp(x_e %*% beta)) %*% Psi_e +
        t(exp(x_rightC %*% beta)) %*% Psi_c_r +
        t(exp(x_c_interval %*% beta)*S_cens_i1/(S_cens_i1-S_cens_i2)) %*% Psi_c_i1 + 
        t(exp(x_trunc_left %*% beta)*S_trunc_R_left/(S_trunc_L_left-S_trunc_R_left)) %*% Psi_trunc_R_left + 
        t(exp(x_trunc_right %*% beta)*S_trunc_R_right/(S_trunc_L_right-S_trunc_R_right)) %*% Psi_trunc_R_right + 
        t(exp(x_trunc_interval %*% beta)*S_trunc_R_interval/(S_trunc_L_interval-S_trunc_R_interval)) %*% Psi_trunc_R_interval + 
        t(TwoLRtheta*(TwoLRtheta>0)) + 1e-4
      
      B_matrix = diag(c(as.numeric(theta) / theta_score_neg))
      theta_old = theta
      theta = theta_old + B_matrix %*% t(theta_score)
      
      ##re calculate log likelihood
      
      S_trunc_L_left = exp(-(Psi_trunc_L_left %*% theta) * exp(x_trunc_left %*% beta))
      S_trunc_R_left = exp(-(Psi_trunc_R_left %*% theta) * exp(x_trunc_left %*% beta))
      
      S_trunc_L_right = exp(-(Psi_trunc_L_right %*% theta) * exp(x_trunc_right %*% beta))
      S_trunc_R_right = exp(-(Psi_trunc_R_right %*% theta) * exp(x_trunc_right %*% beta))
      
      S_trunc_L_interval = exp(-(Psi_trunc_L_interval %*% theta) * exp(x_trunc_interval %*% beta))
      S_trunc_R_interval = exp(-(Psi_trunc_R_interval %*% theta) * exp(x_trunc_interval %*% beta))
      
      log_S_trunc_diff = sum(log(S_trunc_L_left + 1e-6)) + sum(log(S_trunc_L_right - S_trunc_R_right + 1e-5)) + 
        sum(log(S_trunc_L_interval - S_trunc_R_interval + 1e-6))
      
      #likelihood for event times
      
      x_e = cov[which(cens_type == 1),]
      h0t_e = psi_e %*% theta
      H0t_e = Psi_e %*% theta
      Ht_e = H0t_e * exp(x_e %*% beta)
      St_e = exp(-Ht_e)
      
      #likelihood for right censored times
      x_rightC = cov[which(cens_type == 0 ),]
      S_cens_right = exp(-( Psi_c_r %*% theta) * exp(x_rightC %*% beta))
      
      #likelihood for left censored times
      x_c_left = cov[which(cens_type == 2),]
      S_cens_left = exp(-( Psi_c_l %*% theta) * exp(x_c_left %*% beta))
      
      #likelihood for interval censored times
      x_c_interval = cov[which(cens_type == 3),]
      S_cens_i1 = exp(-( Psi_c_i1 %*% theta) * exp(x_c_interval %*% beta))
      S_cens_i2 = exp(-( Psi_c_i2 %*% theta) * exp(x_c_interval %*% beta))
      
      log_lik_old = log_lik
      log_lik = sum(log(h0t_e) + x_e %*% beta - Ht_e) + sum(log(S_cens_right)) + 
        sum(log(1-S_cens_left)) + sum(log(S_cens_i1 - S_cens_i2 + 1e-5)) - 
        log_S_trunc_diff - lambda*t(theta) %*% R %*% theta
      
      if((log_lik < log_lik_old) & step_size[2] == 1){
        
        i = 0
        omega1 = 1/kappa
        
        
        while(log_lik<log_lik_old){
          theta = theta_old + omega1 * B_matrix %*% t(theta_score)
          
          #update log-likelihood
          
          S_trunc_L_left = exp(-(Psi_trunc_L_left %*% theta) * exp(x_trunc_left %*% beta))
          S_trunc_R_left = exp(-(Psi_trunc_R_left %*% theta) * exp(x_trunc_left %*% beta))
          
          S_trunc_L_right = exp(-(Psi_trunc_L_right %*% theta) * exp(x_trunc_right %*% beta))
          S_trunc_R_right = exp(-(Psi_trunc_R_right %*% theta) * exp(x_trunc_right %*% beta))
          
          S_trunc_L_interval = exp(-(Psi_trunc_L_interval %*% theta) * exp(x_trunc_interval %*% beta))
          S_trunc_R_interval = exp(-(Psi_trunc_R_interval %*% theta) * exp(x_trunc_interval %*% beta))
          
          log_S_trunc_diff = sum(log(S_trunc_L_left + 1e-6)) + sum(log(S_trunc_L_right - S_trunc_R_right + 1e-5)) + 
            sum(log(S_trunc_L_interval - S_trunc_R_interval + 1e-6))
          
          #likelihood for event times
          
          x_e = cov[which(cens_type == 1),]
          h0t_e = psi_e %*% theta
          H0t_e = Psi_e %*% theta
          Ht_e = H0t_e * exp(x_e %*% beta)
          St_e = exp(-Ht_e)
          
          #likelihood for right censored times
          x_rightC = cov[which(cens_type == 0 ),]
          S_cens_right = exp(-( Psi_c_r %*% theta) * exp(x_rightC %*% beta))
          
          #likelihood for left censored times
          x_c_left = cov[which(cens_type == 2),]
          S_cens_left = exp(-( Psi_c_l %*% theta) * exp(x_c_left %*% beta))
          
          #likelihood for interval censored times
          x_c_interval = cov[which(cens_type == 3),]
          S_cens_i1 = exp(-( Psi_c_i1 %*% theta) * exp(x_c_interval %*% beta))
          S_cens_i2 = exp(-( Psi_c_i2 %*% theta) * exp(x_c_interval %*% beta))
          
          log_lik = sum(log(h0t_e) + x_e %*% beta - Ht_e) + sum(log(S_cens_right)) + 
            sum(log(1-S_cens_left)) + sum(log(S_cens_i1 - S_cens_i2 + 1e-5)) - 
            log_S_trunc_diff - lambda*t(theta) %*% R %*% theta
          
          #update value of omega1
          if(omega1>=1e-2){
            omega1 = omega1/kappa
          }else if(omega1<1e-2 & omega1>=1e-5){
            omega1 = omega1*(5e-2)
          }else if(omega1<1e-5){
            omega1 = omega1*(1e-5)
          }
          i = i+1
          if(i>1000){break}
          
          
        }
      }
      
      print(c(it, beta, theta, log_lik))
      
      if(all(abs(c(beta - beta_old, theta- theta_old)) < 1e-6)){
        break
        
        }
      
    }
    
    Hess = matrix(0, nrow = (p+m), ncol = (p+m))
    
    ##beta beta
    H_trunc_L_left = -log(S_trunc_L_left)
    H_trunc_R_left = -log(S_trunc_R_left)
    
    H_trunc_L_right = -log(S_trunc_L_right)
    H_trunc_R_right = -log(S_trunc_R_right + 1e-6)
    
    H_trunc_L_interval = -log(S_trunc_L_interval)
    H_trunc_R_interval = -log(S_trunc_R_interval)
    
    H_cens_right = -log(S_cens_right)
    H_cens_left = -log(S_cens_left)
    H_cens_i1 = -log(S_cens_i1)
    H_cens_i2 = -log(S_cens_i2)
    
    beta_hess = -t(x_e) %*% diag(c(Ht_e)) %*% x_e - 
      t(x_rightC) %*% diag(c(H_cens_right)) %*% x_rightC + 
      t(x_c_left) %*% diag(c((S_cens_left*H_cens_left -  S_cens_left*H_cens_left^2 - (S_cens_left^2)*H_cens_left)/(1-S_cens_left)^2)) %*% x_c_left - 
      t(x_c_interval) %*% diag(c(S_cens_i1*H_cens_i1/(S_cens_i1-S_cens_i2 + 1e-6))) %*% x_c_interval + 
      t(x_c_interval) %*% diag(c(S_cens_i2*H_cens_i2/(S_cens_i1-S_cens_i2 + 1e-6))) %*% x_c_interval -
      t(x_c_interval) %*% diag(c(S_cens_i1*S_cens_i2*(H_cens_i1-H_cens_i2)^2/(S_cens_i1-S_cens_i2 + 1e-6)^2)) %*% x_c_interval + 
      t(x_trunc_left) %*% diag(c(S_trunc_L_left*H_trunc_L_left/(S_trunc_L_left-0))) %*% x_trunc_left + 
      t(x_trunc_right) %*% diag(c(S_trunc_L_right*H_trunc_L_right/(S_trunc_L_right-S_trunc_R_right))) %*% x_trunc_right -
      t(x_trunc_right) %*% diag(c(S_trunc_R_right*H_trunc_R_right/(S_trunc_L_right-S_trunc_R_right))) %*% x_trunc_right +
      t(x_trunc_right) %*% diag(c(S_trunc_L_right*S_trunc_R_right*(H_trunc_L_right-H_trunc_R_right)^2/(S_trunc_L_right-S_trunc_R_right)^2)) %*% x_trunc_right + 
      t(x_trunc_interval) %*% diag(c(S_trunc_L_interval*H_trunc_L_interval/(S_trunc_L_interval-S_trunc_R_interval + 1e-5))) %*% x_trunc_interval -
      t(x_trunc_interval) %*% diag(c(S_trunc_R_interval*H_trunc_R_interval/(S_trunc_L_interval-S_trunc_R_interval + 1e-5))) %*% x_trunc_interval +
      t(x_trunc_interval) %*% diag(c(S_trunc_L_interval*S_trunc_R_interval*(H_trunc_L_interval-H_trunc_R_interval)^2/(S_trunc_L_interval-S_trunc_R_interval + 1e-5)^2)) %*% x_trunc_interval
    
    Hess[1:p, 1:p] = beta_hess
    
    ##beta theta
    beta_theta_hess = - t(x_e) %*% (c(exp(x_e %*% beta)) * Psi_e) -
      t(x_rightC) %*% (c(exp(x_rightC %*% beta)) * Psi_c_r) +
      t(x_c_left) %*% (c(exp(x_c_left %*% beta)*(S_cens_left - S_cens_left*H_cens_left - S_cens_left^2)/(1-S_cens_left)) * Psi_c_l) - 
      t(x_c_interval) %*% (c(exp(x_c_interval %*% beta)*(S_cens_i1)/(S_cens_i1-S_cens_i2 + 1e-4)) * Psi_c_i1) +
      t(x_c_interval) %*% (c(exp(x_c_interval %*% beta)*(S_cens_i2)/(S_cens_i1-S_cens_i2 + 1e-4)) * Psi_c_i2) -
      t(x_c_interval) %*% (c(exp(x_c_interval %*% beta)*(S_cens_i1*S_cens_i2*(H_cens_i1-H_cens_i2))/(S_cens_i1-S_cens_i2 + 1e-4)^2) * (Psi_c_i1-Psi_c_i2)) +
      t(x_trunc_left) %*% (c(exp(x_trunc_left %*% beta)*(S_trunc_L_left)/(S_trunc_L_left-S_trunc_R_left)) * Psi_trunc_L_left) + 
      t(x_trunc_right) %*% (c(exp(x_trunc_right %*% beta)*(S_trunc_L_right)/(S_trunc_L_right-S_trunc_R_right + 1e-4)) * Psi_trunc_L_right) -
      t(x_trunc_right) %*% (c(exp(x_trunc_right %*% beta)*(S_trunc_R_right)/(S_trunc_L_right-S_trunc_R_right + 1e-4)) * Psi_trunc_R_right) + 
      t(x_trunc_right) %*% (c(exp(x_trunc_right %*% beta)*(S_trunc_L_right*S_trunc_R_right*(H_trunc_L_right-H_trunc_R_right))/(S_trunc_L_right-S_trunc_R_right + 1e-4)^2) * (Psi_trunc_L_right-Psi_trunc_R_right)) + 
      t(x_trunc_interval) %*% (c(exp(x_trunc_interval %*% beta)*(S_trunc_L_interval)/(S_trunc_L_interval-S_trunc_R_interval)) * Psi_trunc_L_interval) -
      t(x_trunc_interval) %*% (c(exp(x_trunc_interval %*% beta)*(S_trunc_R_interval)/(S_trunc_L_interval-S_trunc_R_interval)) * Psi_trunc_R_interval) + 
      t(x_trunc_interval) %*% (c(exp(x_trunc_interval %*% beta)*(S_trunc_L_interval*S_trunc_R_interval*(H_trunc_L_interval-H_trunc_R_interval))/(S_trunc_L_interval-S_trunc_R_interval)^2) * (Psi_trunc_L_interval-Psi_trunc_R_interval))
    
    Hess[1:p, (p+1):(p+m)] = beta_theta_hess
    Hess[(p+1):(p+m), 1:p] = t(beta_theta_hess)
    
    theta_theta_hess = -t(psi_e) %*% diag(c(1/(h0t_e^2))) %*% psi_e -
      t(Psi_c_l) %*% diag(c(S_cens_left*exp(x_c_left %*% beta)^2/(1-S_cens_left)^2)) %*% Psi_c_l -
      t(Psi_c_i1 - Psi_c_i2) %*% diag(c(S_cens_i1*S_cens_i2*exp(x_c_interval %*% beta)^2/(S_cens_i1-S_cens_i2 + 1e-4)^2)) %*% (Psi_c_i1 - Psi_c_i2) +
      t(Psi_trunc_L_left - Psi_trunc_R_left) %*% diag(c(S_trunc_L_left*0*exp(x_trunc_left %*% beta)^2/(S_trunc_L_left-S_trunc_R_left + 1e-4)^2)) %*% (Psi_trunc_L_left - Psi_trunc_R_left) +
      t(Psi_trunc_L_right - Psi_trunc_R_right) %*% diag(c(S_trunc_L_right*S_trunc_R_right*exp(x_trunc_right %*% beta)^2/(S_trunc_L_right-S_trunc_R_right + 1e-4)^2)) %*% (Psi_trunc_L_right - Psi_trunc_R_right) +
      t(Psi_trunc_L_interval - Psi_trunc_R_interval) %*% diag(c(S_trunc_L_interval*S_trunc_R_interval*exp(x_trunc_interval %*% beta)^2/(S_trunc_L_interval-S_trunc_R_interval + 1e-4)^2)) %*% (Psi_trunc_L_interval - Psi_trunc_R_interval)
    
    Hess[(p+1):(p+m), (p+1):(p+m)] = theta_theta_hess
    
    Hess_pen = Hess
    Hess_pen[(p+1):(p+m), (p+1):(p+m)] = Hess_pen[(p+1):(p+m), (p+1):(p+m)] - lambda*R
    lambdaR_star = matrix(0, nrow = nrow(Hess), ncol = ncol(Hess))
    lambdaR_star[(p+1):(p+m), (p+1):(p+m)] = - lambda*R
    
    H_pen_inv = solve(-(Hess_pen))
    
    df.theta.old = df.theta
    df.theta = sum(diag( -H_pen_inv %*% -(lambdaR_star)))
    df.theta = m - df.theta
    theta.sigma2 = t(theta) %*% R %*% theta / df.theta
    lambda = c(1/(2*theta.sigma2))
    
    
    if(abs(df.theta.old - df.theta) < 1){
      break
    }
    
  }
  
  #Q
  Q_asy = matrix(0, nrow = n, ncol = (p + m))
  
  #beta
  #Q_asy[which(cens_type == 1), 1:p] = (x_e) * (c(1 - Ht_e)) + ((x_trunc_right) * c((S_trunc_L_right*H_trunc_L_right - S_trunc_R_right*H_trunc_R_right)/(S_trunc_L_right-S_trunc_R_right)))[cens_type == 1]
  #Q_asy[which(cens_type == 3), 1:p] = -(x_c_interval) * c((S_cens_i1*H_cens_i1 - S_cens_i2*H_cens_i2)/(S_cens_i1-S_cens_i2)) + ((x_trunc_right) * c((S_trunc_L_right*H_trunc_L_right - S_trunc_R_right*H_trunc_R_right)/(S_trunc_L_right-S_trunc_R_right)))[cens_type == 3]
  
  #Q_asy[which(cens_type == 1), (p+1):(p+m)] = c(1/h0t_e) * psi_e - c(exp(x_e %*% beta)) * Psi_e + 
  #  (c(exp(x_trunc_right %*% beta)*S_trunc_L_right/(S_trunc_L_right-S_trunc_R_right)) * Psi_trunc_L_right)[cens_type == 1,] -
  #  (c(exp(x_trunc_right %*% beta)*S_trunc_R_right/(S_trunc_L_right-S_trunc_R_right)) * Psi_trunc_R_right)[cens_type == 1,]
  #Q_asy[which(cens_type == 3), (p+1):(p+m)] = -c(exp(x_c_interval %*% beta)*S_cens_i1/(S_cens_i1-S_cens_i2)) * Psi_c_i1 + 
  #  c(exp(x_c_interval %*% beta)*S_cens_i2/(S_cens_i1-S_cens_i2)) * Psi_c_i2 +
  #  (c(exp(x_trunc_right %*% beta)*S_trunc_L_right/(S_trunc_L_right-S_trunc_R_right)) * Psi_trunc_L_right)[cens_type == 3,] -
  #  (c(exp(x_trunc_right %*% beta)*S_trunc_R_right/(S_trunc_L_right-S_trunc_R_right)) * Psi_trunc_R_right)[cens_type == 3,]
  
  Q_asy = t(Q_asy)%*%Q_asy
  
  #calculate asymptotic variance
  M2 = -(Hess + lambdaR_star)
  M2_q = (Q_asy + lambdaR_star)
  
  pos = rep(1, nrow(M2))
  for(u in 1:m){
    if(theta[u]< (0.5) & theta_score[u] < (-0.1) | (theta[u]< (1e-8))){
      pos[p+u] = 0
    }
  }
  
  corr = M2q_inv = M2_inv = matrix(0, nrow(M2), nrow(M2))
  diag(corr) = pos
  #corr[!pos,]=0
  
  M2_inv[(pos==1), (pos==1)] = solve(M2[(pos==1), (pos==1)])
  #M2q_inv[(pos==1), (pos==1)] = solve(M2_q[(pos==1), (pos==1)])
  
  A_omega = corr %*% M2_inv %*% t(corr)
  cov_H = A_omega %*% (-(Hess)) %*% A_omega
  
  #A_omega_q = corr %*% M2q_inv %*% t(corr)
  #cov_Q = A_omega_q %*% ((Q_asy)) %*% A_omega_q
  
  
  se_H = sqrt(diag(cov_H))
  se_M2 = sqrt(diag(M2_inv))
  #se_Q = sqrt(diag(cov_Q))
  
  pars = list(beta = beta, theta = theta)
  se = list(se_H = se_H, se_M2 = se_M2)
  cov = list(M2_inv = M2_inv, cov_H = cov_H)
  smoothing = list(lambda = lambda, df.theta = df.theta)
  scores = list(beta_score = beta_score, theta_score = theta_score)
  knots = list(int.knots = int.knots, bound.knots = bound.knots)
  out = list(pars = pars, it = it, se = se, cov = cov, smoothing = smoothing, scores = scores, knots = knots)
 
  return(out)
  
}









theta_penalty_f = function(order, IntKnt, bryKnt, norm = 2){
  if(norm == 2){
    
    ordSp = order
    dgrSp = ordSp - 1
    numIntKnt = length(IntKnt)
    numSp = numIntKnt+ordSp
    
    minTime = min(bryKnt)
    maxTime = max(bryKnt)
    
    R=matrix(0, nrow=numSp, ncol=numSp) 
    xknots = c(rep(minTime, ordSp), IntKnt, rep(maxTime, ordSp))
    for (ii in 1:numSp){
      for (jj in ii:numSp){
        if (jj - ii<ordSp){
          kntset = xknots[xknots>=xknots[jj] & xknots<=xknots[ii+ordSp]]
          kntsum = 0
          for (kk in 1:(length(kntset)-1)){
            kntsum = kntsum + mSpline(kntset[kk], knots=IntKnt, degree=dgrSp, intercept=T, Boundary.knots=bryKnt, 
                                      derivs=dgrSp)[ii]*mSpline(kntset[kk], knots=IntKnt, degree=dgrSp, intercept=T, 
                                                                Boundary.knots=bryKnt,derivs=dgrSp)[jj]*(kntset[kk+1]-kntset[kk])
          }
          R[ii, jj] = kntsum
        }
      }
    }
    R[lower.tri(R, diag = FALSE)] = t(R)[lower.tri(R, diag = FALSE)]
    
  }
  
  
  return(R)
  
  
}






