## An R simulation study

#create empty matrix to save results
save_lefttrunc_sim <- matrix(0, nrow = 300, ncol = 8)

for(s in 1:300){
  #simulate data
  trunc_data_event40  <- trunc_sim_data3.1(nsample = 1000,
                                           beta = c(1, -0.5), 
                                           prob_left_0 = 0, 
                                           right_cens_par1 = 0.55, right_cens_par2 = 1.5)
  
  #fit profile likelihood method
  trunc.y <- trunc_data_event40$trunc_data[,c(3,7,9)]
  design <- trunc_data_event40$trunc_data[,5:6]
  prof.trunc <- mpleltrc(trunc.y, design, 100)
  
  save_lefttrunc_sim[s,1] <- prof.trunc[[1]][1]
  save_lefttrunc_sim[s,2] <- prof.trunc[[1]][2]
  save_lefttrunc_sim[s,3] <- prof.trunc[[5]][1]
  save_lefttrunc_sim[s,4] <- prof.trunc[[5]][2]
  
  #fit partial likelihood method
  pl.trunc <- coxph(Surv(trunc_data_event40$trunc_data$Lcen, 
                         trunc_data_event40$trunc_data$censor_type) ~ 
                      trunc_data_event40$trunc_data$x + trunc_data_event40$trunc_data$X2)
  
  save_lefttrunc_sim[s,5] <- pl.trunc$coefficients[1]
  save_lefttrunc_sim[s,6] <- pl.trunc$coefficients[2]
  save_lefttrunc_sim[s,7] <- sqrt(pl.trunc$var[1,1])
  save_lefttrunc_sim[s,8] <- sqrt(pl.trunc$var[2,2])
  
}
