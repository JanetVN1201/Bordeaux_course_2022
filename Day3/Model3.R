# Bordeaux workshop May 2022
# Model 3: Joint longitudinal - Recurrent - Survival

library(INLA)
library(frailtypack)
set.seed(1)


data(colorectal)
data(colorectalLongi)
summary(colorectalLongi)
colorectalLongi$treatment <- as.integer(colorectalLongi$treatment)-1
colorectal$treatment <- as.integer(colorectal$treatment)-1
colorectalSurv <- subset(colorectal, new.lesions == 0)


# Weibull baseline hazard function
# Random effects as the link function, Gap timescale
model.weib.RE.gap <-trivPenal(Surv(gap.time, new.lesions) ~ cluster(id) + treatment + terminal(state),
                              formula.terminalEvent =~ treatment,
                              tumor.size ~ year * treatment, data = colorectal,
                              data.Longi = colorectalLongi, random = c("1", "year"), id = "id",
                              link = "Random-effects", left.censoring = -3.33, 
                              hazard = "Weibull", method.GH="Pseudo-adaptive",n.nodes=20)

model.weib.RE.gap





NL <- dim(colorectalLongi)[1]
NR <- dim(colorectal)[1]
NS <- dim(colorectalSurv)[1]
data <- list(
  IntL = c(rep(1, NL), rep(NA, NR), rep(NA, NS)),
  TimeL = c(colorectalLongi$year, rep(NA, NR), rep(NA, NS)),
  TrtL = c(colorectalLongi$treatment, rep(NA, NR), rep(NA, NS)),
  TimeTrtL = c(colorectalLongi$year*colorectalLongi$treatment, rep(NA, NR), rep(NA, NS)),
  IntR = c(rep(NA, NL), rep(1, NR), rep(NA, NS)),
  TrtR = c(rep(NA, NL), colorectal$treatment, rep(NA, NS)),
  IntS = c(rep(NA, NL), rep(NA, NR), rep(1, NS)),
  TrtS = c(rep(NA, NL), rep(NA, NR), colorectalSurv$treatment),
  b_i0L = c(colorectalLongi$id, rep(NA, NR), rep(NA, NS)),
  b_i1L = c(NS+colorectalLongi$id, rep(NA, NR), rep(NA, NS)),
  b_i0R = c(rep(NA, NL), colorectal$id, rep(NA, NS)),
  b_i1R = c(rep(NA, NL), NS+colorectal$id, rep(NA, NS)),
  b_i0S = c(rep(NA, NL), rep(NA, NR), colorectalSurv$id),
  b_i1S = c(rep(NA, NL), rep(NA, NR), NS+colorectalSurv$id),
  f_i0R = c(rep(NA, NL), colorectal$id, rep(NA, NS)),
  f_i0S = c(rep(NA, NL), rep(NA, NR), colorectalSurv$id),
  Yjoint = list(Y = inla.surv(time=c(exp(colorectalLongi$tumor.size), rep(NA, NR), rep(NA, NS)),
                              event=c(ifelse(colorectalLongi$tumor.size < -3.32, 2, 1))),
                R = inla.surv(time=c(rep(NA, NL), colorectal$gap.time, rep(NA, NS)),
                              event=c(rep(NA, NL), colorectal$new.lesions, rep(NA, NS))),
                S = inla.surv(time=c(rep(NA, NL), rep(NA, NR), colorectalSurv$time1),
                              event=c(rep(NA, NL), rep(NA, NR), colorectalSurv$state)))
)

# INLA fit
formula <- Yjoint ~ -1 + IntL + TimeL + TrtL + TimeTrtL + IntR + TrtR + IntS + TrtS +
  f(b_i0L, model="iidkd", order=2, n=NS*2, hyper = list(theta1 = list(param = c(10,1,1,0)))) +
  #f(b_i0L, model="iid", hyper=list(prec=list(prior="loggamma", param=c(0.01,0.01)))) +
  f(b_i1L, TimeL, copy="b_i0L") +
  f(b_i0R, copy="b_i0L", fixed=FALSE) +
  f(b_i1R, copy="b_i0L", fixed=FALSE) +
  f(b_i0S, copy="b_i0L", fixed=FALSE) +
  f(b_i1S, copy="b_i0L", fixed=FALSE) +
  f(f_i0R, model="iid")+#, hyper=list(prec=list(prior="loggamma", param=c(0.01,0.01)))) +
  f(f_i0S, copy="f_i0R", fixed=FALSE)

M3 <- inla(formula, data=data, family=c("lognormalsurv", "weibullsurv", "weibullsurv"), 
           control.family=list(list(),list(variant=1),list(variant=1)),
           control.inla = list(int.strategy="eb"))
summary(M3)

#control.compute=list(config=T) # to draw sample from the model


MC_samples <- inla.iidkd.sample(10^4, M3, "b_i0L", return.cov=T)
VarCov <- matrix(unlist(MC_samples), nrow = 2^2)
VarCovMeans = matrix(rowMeans(VarCov),2,2);round(VarCovMeans,2)
VarCovSD = matrix(apply(VarCov, 1, sd),2,2);round(VarCovSD,2)
VarCov025 = matrix(apply(VarCov, 1, function(x) quantile(x, 0.025)),2,2);round(VarCov025,2)
VarCov05 = matrix(apply(VarCov, 1, function(x) quantile(x, 0.5)),2,2);round(VarCov05,2)
VarCov975 = matrix(apply(VarCov, 1, function(x) quantile(x, 0.975)),2,2);round(VarCov975,2)


MAR <- inla.tmarginal(function(x) 1/exp(x), marginal = M3$internal.marginals.hyperpar$`Log precision for the lognormalsurv observations`)
inla.zmarginal(MAR)

MAR <- inla.tmarginal(function(x) 1/exp(x), marginal = M3$internal.marginals.hyperpar$`Log precision for f_i0R`)
inla.zmarginal(MAR)

MAR <- inla.tmarginal(function(x) exp(x), marginal = M3$marginals.fixed$IntS)
inla.zmarginal(MAR)
MAR <- inla.tmarginal(function(x) exp(x), marginal = M3$marginals.fixed$IntR)
inla.zmarginal(MAR)








