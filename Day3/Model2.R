# Bordeaux workshop May 2022
# Model 2: Joint longitudinal - Survival

library(INLA)
set.seed(1)
# data generation - one longitudinal marker
nsujet=500 # number of individuals

# Y1 (continuous)
b_0=0.2 # intercept
b_1=-0.1 # slope
b_e=0.1 # residual error

# S1 (survival)
s_1=0.2 # ctsX
s_2=-0.2 # binX
phi_b0=1 # random intercept association with survival
phi_b1=1 # random slope association with survival

gap=1 # gap between measurements
followup=5 # follow-up time
mestime=seq(0,followup,gap) # measurement times
time=rep(mestime, nsujet) # time column
nmesindiv=followup/gap+1 # number of individual measurements
nmesy= nmesindiv*nsujet # total number of measurements
id<-rep(1:nsujet, each=nmesindiv) # individual id

# random effects variance and covariance matrix
b0 <- 0.2
b1 <- 0.3
corb0b1 <- 0.5
covb0b1 <- b0*b1*corb0b1
Sigma <- matrix(c(b0^2, covb0b1, covb0b1, b1^2),ncol=2,nrow=2)
MVnorm <- mvtnorm::rmvnorm(nsujet, rep(0, 2), Sigma)
b_i0 <- rep(MVnorm[,1], each=nmesindiv) # random intercept Y1
b_i1 <- rep(MVnorm[,2], each=nmesindiv) # random slope Y1

# linear predictor
linPredY1 <- (b_i0+b_0) + (b_i1+b_1)*time
# continuous outcome Y1
Y1 <- rnorm(nmesy, linPredY1, b_e)
lon <- data.frame(Y1, id, time)
summary(lon) # dataset

ctsX=rnorm(nsujet,1, 0.5) # continuous covariate
binX=rbinom(nsujet,1, 0.5) # binary covariate
## generation of exponential death times
u <- runif(nsujet) # uniform distribution for survival times generation
baseScale=0.2
deathTimes <- -(log(u) / (baseScale * exp(ctsX*s_1 + binX*s_2 + MVnorm[,1]*phi_b0 + MVnorm[,2]*phi_b1)))
d <- as.numeric(deathTimes<followup) # deathtimes indicator
## censoring individuals at end of follow-up (not at random)
deathTimes[deathTimes>=followup]=followup
surv <- data.frame(id=1:nsujet,deathTimes, d, ctsX, binX) # survival times dataset

## removing longi measurements after death
ind <- rep(NA, nsujet*length(mestime))
for (i in 1:nsujet){
  for(j in 1:length(mestime)){
    if(lon[(i-1)*length(mestime)+j, "time"]<=surv[i,"deathTimes"]) ind[(i-1)*length(mestime)+j]=1
  }
}
lon <- lon[!is.na(ind),]
summary(surv)
summary(lon)









