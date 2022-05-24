# Bordeaux workshop May 2022
# Model 1: Multivariate longitudinal

library(INLA)
set.seed(1)
# data generation - two longitudinal markers
nsujet=500 # number of individuals

# Y1 (continuous)
b1_0=0.2 # intercept
b1_1=-0.1 # slope
b1_e=0.1 # residual error

# Y2 (counts)
b2_0=3 # intercept
b2_1=-0.1 # slope

gap=1 # gap between measurements
followup=5 # follow-up time
mestime=seq(0,followup,gap) # measurement times
time=rep(mestime, nsujet) # time column
nmesindiv=followup/gap+1 # number of individual measurements
nmesy= nmesindiv*nsujet # total number of measurements
id<-rep(1:nsujet, each=nmesindiv) # individual id

# random effects variance-covariance matrix
Sigma <- matrix(c(0.16, 0.03, 0.02, 0.04,
                  0.03, 0.09, 0.03, 0.00,
                  0.02, 0.03, 0.25, 0.08,
                  0.04, 0.00, 0.08, 0.16),ncol=4,nrow=4)

MVnorm <- mvtnorm::rmvnorm(nsujet, rep(0, 4), Sigma)
b1_i0 <- rep(MVnorm[,1], each=nmesindiv) # random intercept Y1
b1_i1 <- rep(MVnorm[,2], each=nmesindiv) # random slope Y1
b2_i0 <- rep(MVnorm[,3], each=nmesindiv) # random intercept Y2
b2_i1 <- rep(MVnorm[,4], each=nmesindiv) # random slope Y2

# linear predictor
linPredY1 <- (b1_i0+b1_0) + (b1_i1+b1_1)*time
linPredY2 <- (b2_i0+b2_0) + (b2_i1+b2_1)*time

# observed outcomes:
# continuous outcome Y1
Y1 <- rnorm(nmesy, linPredY1, b1_e)
# count outcome Y2
Y2 <- rpois(nmesy, exp(linPredY2))
lon <- data.frame(Y1, Y2, id, time)
summary(lon) # dataset







