# Read in dataset
parkinsons <- read.csv("pdearly_8_12_09.csv")

# Remove rows with missing observations
parkinsons <- parkinsons[-which(is.na(parkinsons$Sampling.Age)),]

# Load package/code for fitting MPL double truncation model
source("interval_trunc_code/fitting_interval_trunc1.R")
source("interval_trunc_code/pointwiseCIs.R")

#### Assume exact event time
parkinsons$t_L <- parkinsons$Onset.Age
parkinsons$t_R <- parkinsons$Onset.Age

# Refine truncation times
parkinsons$trunc_time_L <- parkinsons$Sampling.Age - 8
parkinsons$trunc_time_R <- parkinsons$Sampling.Age

#Univariate analysis for X10398
design <- cbind(as.numeric(parkinsons$X10398 == "G"))
design <- as.matrix(design)

parkinsons.mod1 <- CoxMPL_Trunc(parkinsons$t_L, parkinsons$t_R, 
             design, parkinsons$trunc_time_L, parkinsons$trunc_time_R, 
             max.iter = c(10, 2000), n.knots.h0 = 4, step_size=c(1,1))
parkinsons.mod1$pars$beta
parkinsons.mod1$se$se_H[1]
parkinsons.mod1$pars$beta - 1.96*parkinsons.mod1$se$se_H[1]
parkinsons.mod1$pars$beta + 1.96*parkinsons.mod1$se$se_H[1]


#Univariate analysis for PGC1A
design <- cbind(as.numeric(parkinsons$PGC1A_GLY482SER == "AG"),
               as.numeric(parkinsons$PGC1A_GLY482SER == "G"))
design <- as.matrix(design)

parkinsons.mod2 <- CoxMPL_Trunc(parkinsons$t_L, parkinsons$t_R, 
                               design, parkinsons$trunc_time_L, parkinsons$trunc_time_R, 
                               max.iter = c(10, 2000), n.knots.h0 = 4, step_size=c(1,1))

parkinsons.mod2$pars$beta
parkinsons.mod2$se$se_H[1:2]
parkinsons.mod2$pars$beta - 1.96*parkinsons.mod2$se$se_H[1:2]
parkinsons.mod2$pars$beta + 1.96*parkinsons.mod2$se$se_H[1:2]


#Multivariate analysis
design <- cbind(as.numeric(parkinsons$X10398 == "G"),
               as.numeric(parkinsons$PGC1A_GLY482SER == "AG"),
               as.numeric(parkinsons$PGC1A_GLY482SER == "G"))
design <- as.matrix(design)

parkinsons.mod3 <- CoxMPL_Trunc(parkinsons$t_L, parkinsons$t_R, 
                               design, parkinsons$trunc_time_L, parkinsons$trunc_time_R, 
                               max.iter = c(10, 2000), n.knots.h0 = 4, step_size=c(1,1))
parkinsons.mod3$pars$beta
parkinsons.mod3$se$se_H[1:3]
parkinsons.mod3$pars$beta - 1.96*parkinsons.mod3$se$se_H[1:3]
parkinsons.mod3$pars$beta + 1.96*parkinsons.mod3$se$se_H[1:3]



###Refit assuming interval censoring
parkinsons$t_L <- parkinsons$Onset.Age
parkinsons$t_R <- parkinsons$Onset.Age + 0.99
parkinsons$trunc_time_R <- parkinsons$Sampling.Age + 0.99
parkinsons$trunc_time_L <- parkinsons$Sampling.Age - 8


#Univariate analysis for X10398
design <- cbind(as.numeric(parkinsons$X10398 == "G"))
design <- as.matrix(design)

parkinsons.mod1 <- CoxMPL_Trunc(parkinsons$t_L, parkinsons$t_R, 
                               design, parkinsons$trunc_time_L, parkinsons$trunc_time_R, 
                               max.iter = c(10, 2000), n.knots.h0 = 5, step_size=c(1,1))
parkinsons.mod1$pars$beta
parkinsons.mod1$se$se_H[1]
parkinsons.mod1$pars$beta - 1.96*parkinsons.mod1$se$se_H[1]
parkinsons.mod1$pars$beta + 1.96*parkinsons.mod1$se$se_H[1]

#Univariate analysis for PGC1A
design <- cbind(as.numeric(parkinsons$PGC1A_GLY482SER == "AG"),
               as.numeric(parkinsons$PGC1A_GLY482SER == "G"))
design <- as.matrix(design)

parkinsons.mod2 <- CoxMPL_Trunc(parkinsons$t_L, parkinsons$t_R, 
                               design, parkinsons$trunc_time_L, parkinsons$trunc_time_R, 
                               max.iter = c(10, 2000), n.knots.h0 = 5, step_size=c(1,1))

parkinsons.mod2$pars$beta
parkinsons.mod2$se$se_H[1:2]
parkinsons.mod2$pars$beta - 1.96*parkinsons.mod2$se$se_H[1:2]
parkinsons.mod2$pars$beta + 1.96*parkinsons.mod2$se$se_H[1:2]

# Multivariate analysis
design <- cbind(as.numeric(parkinsons$X10398 == "G"),
               as.numeric(parkinsons$PGC1A_GLY482SER == "AG"),
               as.numeric(parkinsons$PGC1A_GLY482SER == "G"))
design <- as.matrix(design)

parkinsons.mod4 <- CoxMPL_Trunc(parkinsons$t_L, parkinsons$t_R, 
                               design, parkinsons$trunc_time_L, parkinsons$trunc_time_R, 
                               max.iter = c(10, 2000), n.knots.h0 = 5, step_size=c(1,1))
parkinsons.mod4$pars$beta
parkinsons.mod4$se$se_H[1:3]
parkinsons.mod4$pars$beta - 1.96*parkinsons.mod4$se$se_H[1:3]
parkinsons.mod4$pars$beta + 1.96*parkinsons.mod4$se$se_H[1:3]


# Create predicted survival plots

CI.AA <- trunc_pointwise_CI_St(parkinsons.mod4, seq(0,56,0.1), c(0,0,0))

plot(CI.AA$est_St ~ CI.AA$v, type = "l", ylim = c(0,1),
     xlab = "Age (years)", ylab = "Predicted disease-free survival",
     main = "SNP10398 = A; PGC-1a = A")
polygon(x = c(c(0,CI.AA$v), rev(c(0.01,CI.AA$v))), 
        y = c(c(1,CI.AA$lower), rev(c(1.01,CI.AA$upper))), 
        border = NA, col = adjustcolor("grey", alpha.f=0.5) )
lines(CI.AA$est_St ~ CI.AA$v)


CI.GA <- trunc_pointwise_CI_St(parkinsons.mod4, seq(0,56,0.1), c(1,0,0))

plot(CI.GA$est_St ~ CI.GA$v, type = "l", ylim = c(0,1),
     xlab = "Age (years)", ylab = "Predicted disease-free survival",
     main = "SNP10398 = G; PGC-1a = A")
polygon(x = c(c(0,CI.GA$v), rev(c(0.01,CI.GA$v))), 
        y = c(c(1,CI.GA$lower), rev(c(1.01,CI.GA$upper))), 
        border = NA, col = adjustcolor("grey", alpha.f=0.5) )
lines(CI.GA$est_St ~ CI.GA$v)



CI.AAG <- trunc_pointwise_CI_St(parkinsons.mod4, seq(0,56,0.1), c(0,1,0))

plot(CI.AAG$est_St ~ CI.AAG$v, type = "l", ylim = c(0,1),
     xlab = "Age (years)", ylab = "Predicted disease-free survival",
     main = "SNP10398 = A; PGC-1a = AG")
polygon(x = c(c(0,CI.AAG$v), rev(c(0.01,CI.AAG$v))), 
        y = c(c(1,CI.AAG$lower), rev(c(1.01,CI.AAG$upper))), 
        border = NA, col = adjustcolor("grey", alpha.f=0.5) )
lines(CI.AAG$est_St ~ CI.AAG$v)


CI.GAG <- trunc_pointwise_CI_St(parkinsons.mod4, seq(0,56,0.1), c(1,1,0))

plot(CI.GAG$est_St ~ CI.GAG$v, type = "l", ylim = c(0,1),
     xlab = "Age (years)", ylab = "Predicted disease-free survival",
     main = "SNP10398 = G; PGC-1a = AG")
polygon(x = c(c(0,CI.GAG$v), rev(c(0.01,CI.GAG$v))), 
        y = c(c(1,CI.GAG$lower), rev(c(1.01,CI.GAG$upper))), 
        border = NA, col = adjustcolor("grey", alpha.f=0.5) )
lines(CI.GAG$est_St ~ CI.GAG$v)



CI.AG <- trunc_pointwise_CI_St(parkinsons.mod4, seq(0,56,0.1), c(0,0,1))

plot(CI.AG$est_St ~ CI.AG$v, type = "l", ylim = c(0,1),
     xlab = "Age (years)", ylab = "Predicted disease-free survival",
     main = "SNP10398 = A; PGC-1a = G")
polygon(x = c(c(0,CI.AG$v), rev(c(0.01,CI.AG$v))), 
        y = c(c(1,CI.AG$lower), rev(c(1.01,CI.AG$upper))), 
        border = NA, col = adjustcolor("grey", alpha.f=0.5) )
lines(CI.AG$est_St ~ CI.AG$v)


CI.GG <- trunc_pointwise_CI_St(parkinsons.mod4, seq(0,56,0.1), c(1,0,1))

plot(CI.GG$est_St ~ CI.GG$v, type = "l", ylim = c(0,1),
     xlab = "Age (years)", ylab = "Predicted disease-free survival",
     main = "SNP10398 = G; PGC-1a = G")
polygon(x = c(c(0,CI.GG$v), rev(c(0.01,CI.GG$v))), 
        y = c(c(1,CI.GG$lower), rev(c(1.01,CI.GG$upper))), 
        border = NA, col = adjustcolor("grey", alpha.f=0.5) )
lines(CI.GG$est_St ~ CI.GG$v)




