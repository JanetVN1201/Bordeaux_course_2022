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


NL <- dim(lon)[1]
NS <- dim(surv)[1]

data <- list(
  IntL1 = c(rep(1, NL), rep(NA, NS)),
  TimeL1 = c(lon$time, rep(NA, NS)),
  IntS1 = c(rep(NA, NL), rep(1, NS)),
  ctsXS1 = c(rep(NA, NL), surv$ctsX),
  binXS1 = c(rep(NA, NL), surv$binX),
  b_i0L1 = c(lon$id, rep(NA, NS)),
  b_i1L1 = c(NS+lon$id, rep(NA, NS)),
  b_i0S1 = c(rep(NA, NL), surv$id),
  b_i1S1 = c(rep(NA, NL), NS+surv$id),
  Yjoint = list(Y1 = c(lon$Y1, rep(NA, NS)),
                 S1 = inla.surv(c(rep(NA, NL), surv$deathTimes), c(rep(NA, NL), surv$d)))
)

formula <- Yjoint ~ -1 + IntL1 + TimeL1 + IntS1 + ctsXS1 + binXS1 +
  f(b_i0L1, model="iidkd", order=2, n=2*NS, hyper=list(theta1=list(param=c(10,1,1,0)))) +
  f(b_i1L1, TimeL1, copy="b_i0L1") +
  f(b_i0S1, copy="b_i0L1", fixed=FALSE) +
  f(b_i1S1, copy="b_i0L1", fixed=FALSE)

M3 <- inla(formula, data=data, family=c("gaussian", "exponentialsurv"), control.inla=list(int.strategy="eb"))
summary(M3)

lambda <- inla.tmarginal(function(x) exp(x), marginal = M3$marginals.fixed$IntS1)
inla.zmarginal(lambda)

resErr <- inla.tmarginal(function(x) sqrt(1/x), marginal = M3$marginals.hyperpar$`Precision for the Gaussian observations`)
inla.zmarginal(resErr)


MC_samples <- inla.iidkd.sample(10^4, M3, "b_i0L1", return.cov = F)
VarCov <- matrix(unlist(MC_samples), nrow=2*2)
VarCovMeans <- matrix(rowMeans(VarCov), 2, 2); round(VarCovMeans, 3)
VarCovSD <- matrix(apply(VarCov, 1, sd), 2, 2);round(VarCovSD, 3)
VarCov025 <- matrix(apply(VarCov, 1, function(x) quantile(x, 0.025)), 2, 2) ; round(VarCov025, 3)
VarCov05 <- matrix(apply(VarCov, 1, function(x) quantile(x, 0.5)), 2, 2) ; round(VarCov05, 3)
VarCov975 <- matrix(apply(VarCov, 1, function(x) quantile(x, 0.975)), 2, 2) ; round(VarCov975, 3)



set.seed(1)
M3w <- inla(formula, data=data, family=c("gaussian", "weibullsurv"), control.inla=list(int.strategy="eb"),
            control.family=list(list(), list(variant=1)), num.threads=1)
summary(M3w)







