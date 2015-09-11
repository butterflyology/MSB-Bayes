# Bayesian estimation of population size for Neonympha mitchellii mitchellii

library("R2jags")
library("coda")
library("mcmcplots")

set.seed(862476501)
setwd("~/Desktop/Projects/MSB-Bayes")

# save(list = ls(), file = "MSB-Bayes.RData")
# load("MSB-Bayes.RData")
(sessID <- sessionInfo())



MSB <- read.csv("Neonympha-Bayes.csv", header = TRUE)
str(MSB)
summary(MSB)


MSB.array <- array(NA, dim = c(27, 1, 5))
for(k in 1:5){
	select <- MSB$day == k
	MSB.array[, , k] <- as.matrix(MSB)[select, 3]
}
MSB.array
str(MSB.array)

day.max <- apply(MSB.array, c(1, 3), max, na.rm = TRUE)
day.max
site.max <- apply(day.max, 1, max, na.rm = TRUE)
table(site.max)
plot(table(site.max))
table(site.max > 0)

max1 <- apply(MSB.array, c(1, 3), max)
obs.max.sum <- apply(max1, 2, sum, na.rm = TRUE)
obs.max.sum


#####
##### Poisson model
#####

Poisson.model <- function(){
# Priors
for (k in 1:5){
   alpha.lam[k] ~ dnorm(0, 0.01)
   p[k] ~ dunif(0, 1)
   }

# Likelihood
# Ecological model for true abundance
for (k in 1:5){                          # Loop over days (5)
   lambda[k] <- exp(alpha.lam[k])
   for (i in 1:R){                       # Loop over R sites (23)
      N[i,k] ~ dpois(lambda[k])          # Abundance

      # Observation model for replicated counts
      for (j in 1:T){                    # Loop over temporal reps (1)
         y[i,j,k] ~ dbin(p[k], N[i,k])   # Detection

         # Assess model fit using Chi-squared discrepancy
         # Compute fit statistic E for observed data
         eval[i,j,k] <- p[k] * N[i,k]   	# Expected values
         E[i,j,k] <- pow((y[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k] + 0.5)
         # Generate replicate data and compute fit stats for them
         y.new[i,j,k] ~ dbin(p[k], N[i,k])
         E.new[i,j,k] <- pow((y.new[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k] + 0.5)

         } #j 
      } #i 
   } #k 

# Derived and other quantities
for (k in 1:5){
   totalN[k] <- sum(N[,k])	# Total pop. size across all sites
   mean.abundance[k] <- exp(alpha.lam[k])
   }
fit <- sum(E[,,])
fit.new <- sum(E.new[,,])
}

MSB.dat <- list(y = MSB.array, R = nrow(MSB.array), T = ncol(MSB.array))

# Initial values
Nst.MP <- apply(MSB.array, c(1, 3), max) + 1
Nst.MP[is.na(Nst.MP)] <- 1
inits.MP <- function(){list(N = Nst.MP, alpha.lam = runif(5, -1, 1))}
params.MP <- c("totalN", "mean.abundance", "alpha.lam", "p", "fit", "fit.new")


MSB.Poisson <- jags(MSB.dat, inits = inits.MP, parameters = params.MP, n.thin = 10, n.iter = 1e5, n.burnin = 1e4, n.chains = 4, model.file = Poisson.model)
MSB.Poisson
str(MSB.Poisson)


plot(MSB.Poisson$BUGSoutput$sims.list$fit, MSB.Poisson$BUGSoutput$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE, col = rgb(0, 0, 0, 0.1), pch = 19)
abline(0, 1, lwd = 3, lty = 2, col = "red")
mean(MSB.Poisson$BUGSoutput$sims.list$fit.new > MSB.Poisson$BUGSoutput$sims.list$fit)
mean(MSB.Poisson$BUGSoutput$sims.list$fit) / mean(MSB.Poisson$BUGSoutput$sims.list$fit.new) # bad fit


#####
##### Zero-inflated Poisson model
#####

ZiP.model <- function(){
# Priors
omega ~ dunif(0, 1)
for (k in 1:5){
   alpha.lam[k] ~ dnorm(0, 0.01)
   p[k] ~ dunif(0, 1)
   }

# Likelihood
# Ecological model for true abundance
for (i in 1:R){                         # Loop over R sites (95)
   z[i] ~ dbern(omega)                  # Latent suitability state
   for (k in 1:5){                      # Loop over survey periods (seasons)
      N[i,k] ~ dpois(lam.eff[i,k])      # Latent abundance state
      lam.eff[i,k] <- z[i] * lambda[i,k]
      log(lambda[i,k]) <- alpha.lam[k]
      # Observation model for replicated counts
      for (j in 1:T){                    # Loop over temporal reps (1)
         y[i,j,k] ~ dbin(p[k], N[i,k])   # Detection
         # Assess model fit using Chi-squared discrepancy
         # Compute fit statistic for observed data
         eval[i,j,k] <- p[k] * N[i,k]
         E[i,j,k] <- pow((y[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k] + 0.5)
         # Generate replicate data and compute fit stats for them
         y.new[i,j,k] ~ dbin(p[k], N[i,k])
         E.new[i,j,k] <- pow((y.new[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k]+0.5)
         } #j
      } #k
   } #i

# Derived and other quantities
for (k in 1:5){
   totalN[k] <- sum(N[,k])	# Estimate total pop. size across all sites
   mean.abundance[k] <- exp(alpha.lam[k])
   }
fit <- sum(E[,,])
fit.new <- sum(E.new[,,])
}

# Initial values
Nst.ZiP <- apply(MSB.array, c(1, 3), max) + 1
Nst.ZiP[is.na(Nst.ZiP)] <- 1
inits.ZiP <- function(){list(N = Nst.ZiP, alpha.lam = runif(5, -1, 1))}


# Parameters
params.ZiP <- c("omega", "totalN", "alpha.lam", "p", "mean.abundance", "fit", "fit.new")

MSB.Poisson <- jags(MSB.dat, inits = inits.ZiP, parameters = params.ZiP, n.thin = 10, n.iter = 1e5, n.burnin = 1e4, n.chains = 4, model.file = ZiP.model) # needs tweaking to run


#####
##### Binomial mixture model with overdispersion
#####

BMM.model <- function(){
# Priors
for (k in 1:5){
    alpha.lam[k] ~ dnorm(0, 0.1)
    beta[k] ~ dnorm(0, 0.1)
   }

# Abundance site and detection site-by-day random effects
for (i in 1:R){
   eps[i] ~ dnorm(0, tau.lam)                    # Abundance noise
   }
tau.lam <- 1 / (sd.lam * sd.lam)
sd.lam ~ dunif(0, 3)
tau.p <- 1 / (sd.p * sd.p)
sd.p ~ dunif(0, 3)

# Likelihood
# Ecological model for true abundance
for (i in 1:R){                                 # Loop over R sites (27)
   for (k in 1:5){                              # Loop over days (7)
      N[i,k] ~ dpois(lambda[i,k])               # Abundance
      log(lambda[i,k]) <- alpha.lam[k] + eps[i]

      # Observation model for replicated counts
      for (j in 1:T){                           # Loop over temporal reps (1)
         y[i,j,k] ~ dbin(p[i,j,k], N[i,k])      # Detection
         p[i,j,k] <- 1 / (1 + exp(-lp[i,j,k])) 
         lp[i,j,k] ~ dnorm(beta[k], tau.p) # random delta defined implicitly

         # Assess model fit using Chi-squared discrepancy
         # Compute fit statistic for observed data
         eval[i,j,k] <- p[i,j,k] * N[i,k]
         E[i,j,k] <- pow((y[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k]+0.5)
         # Generate replicate data and compute fit stats for them
         y.new[i,j,k] ~ dbin(p[i,j,k], N[i,k])
         E.new[i,j,k] <- pow((y.new[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k]+0.5)
         } #j
         ik.p[i,k] <- mean(p[i,,k])
      } #k
   } #i

# Derived and other quantities
for (k in 1:5){
   totalN[k] <- sum(N[,k])   # Estimate total pop. size across all sites
   mean.abundance[k] <- mean(lambda[,k])
   mean.N[k] <- mean(N[,k])
   mean.detection[k] <- mean(ik.p[,k])
   }
fit <- sum(E[,,])
fit.new <- sum(E.new[,,])
}

inits.BMM <- function(){list(N = Nst.MP, alpha.lam = runif(5, -3, 3), beta = runif(5, -3, 3), sd.lam = runif(1, 0, 1), sd.p = runif(1, 0, 1))}

# Parameters monitored
params.BMM <- c("totalN", "alpha.lam", "beta", "sd.lam", "sd.p", "mean.abundance", "mean.N", "mean.detection", "fit", "fit.new")

MSB.BMM <- jags(MSB.dat, inits = inits.BMM, parameters = params.BMM, n.thin = 3e2, n.iter = 3.5e5, n.burnin = 5e4, n.chains = 3, model.file = BMM.model)
MSB.BMM 
str(MSB.BMM)

plot(MSB.BMM)
BMM.mcmc <- as.mcmc(MSB.BMM)
summary(BMM.mcmc)
# autocorr.plot(BMM.mcmc)
denplot(BMM.mcmc)
denplot(BMM.mcmc, parms = c("totalN"))



plot(MSB.BMM$BUGSoutput$sims.list$fit, MSB.BMM$BUGSoutput$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE, xlim = c(0, 150), ylim = c(0, 150), las = 1, pch = 19, col = rgb(0, 0, 0, 0.1))
abline(0, 1, lwd = 3, col = "black")
mean(MSB.BMM$BUGSoutput$sims.list$fit.new > MSB.BMM$BUGSoutput$sims.list$fit)
mean(MSB.BMM$BUGSoutput$mean$fit) / mean(MSB.BMM$BUGSoutput$mean$fit.new) # pretty goodfit

hist(MSB.BMM$BUGSoutput$sims.list$mean.abundance[, 4], breaks = 40, main = "Posterior abundance", xlab = "Mean abundance") # Note the skew

median(MSB.BMM$BUGSoutput$sims.list$mean.abundance[, 4])

max.day.count <- apply(MSB.array, c(1, 3), max, na.rm = TRUE)
max.day.count[max.day.count == "-Inf"] <- NA
mean.max.count <- apply(max.day.count, c(2), mean, na.rm = TRUE)
mean.max.count

par(mfrow = c(2, 1))
plot(1:5, mean.max.count, xlab = "Day", ylab = "Mean abundance by section", las = 1, ylim = c(0, 100), type = "b", main = "", frame.plot = FALSE, pch = 19, lwd = 2)
lines(1:5, MSB.BMM$BUGSoutput$summary[14:18, 5], type = "b", pch = 19, col = "blue", lwd = 2)
segments(1:5, MSB.BMM$BUGSoutput$summary[14:18, 3], 1:5, MSB.BMM$BUGSoutput$summary[14:18, 7], col = "blue")



dim(MSB.BMM$BUGSoutput$summary)
MSB.BMM$BUGSoutput$summary[14:18, 5]

