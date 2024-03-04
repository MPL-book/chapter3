trunc_pointwise_CI_St = function(mod, v, x){
  
  Psi.fit = mSpline(v, degree = 3, knots = mod$knots$int.knots, 
                    Boundary.knots = mod$knots$bound.knots, 
                    integral = T, intercept = F)
  theta = mod$pars$theta
  est_Ht = (Psi.fit %*% theta)^exp(x[1]*c(mod$pars$beta)[1])
  est_St = exp(-est_Ht)
  
  lower = upper = NULL
  
  for(t in 1:length(v)){
    
    logit_St = log(est_St[t]/(1-est_St[t]))
    
    dlogitSt_dBeta = -x * est_St[t] * est_Ht[t] * (1/est_St[t] + 1/(1-est_St[t]))
    
    dlogitSt_dTheta = Psi.fit[t,] * c(-est_St[t] * exp(x%*%mod$pars$beta) * (1/est_St[t] + 1/(1-est_St[t])))
    
    dlogitSt_dEta = c(dlogitSt_dBeta, dlogitSt_dTheta)
    
    var_logit_St = t(dlogitSt_dEta) %*% mod$cov$cov_H %*% dlogitSt_dEta
    
    lower_logit = logit_St - 1.96*sqrt(var_logit_St)
    upper_logit = logit_St + 1.96*sqrt(var_logit_St)
    
    lower = c(lower, exp(lower_logit)/(1 + exp(lower_logit)))
    upper = c(upper, exp(upper_logit)/(1 + exp(upper_logit)))
    
  }

  out = data.frame(v, est_St, lower, upper)
  return(out)
  
}






