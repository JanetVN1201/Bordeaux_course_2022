
devtools::install_github('DenisRustand/INLAjoint')



library(INLAjoint)
library(JM) # This package contains the dataset

data(pbc2) # dataset
pbc2 <- pbc2[which(pbc2$id %in% c(1:20)),]
# extract some variable of interest without missing values
Longi <- na.omit(pbc2[, c("id", "years", "status","drug","age", 
                          "sex","year","serBilir","SGOT", "albumin", "edema",
                          "platelets", "alkaline","spiders", "ascites")])
Surv <- Longi[c(which(diff(as.numeric(Longi[,which(colnames(Longi)=="id")]))==1),
                length(Longi[,which(colnames(Longi)=="id")])),-c(7:10, 12:16)]
Surv$death <- ifelse(Surv$status=="dead",1,0) # competing event 1
Surv$trans <- ifelse(Surv$status=="transplanted",1,0) # competing event 2

summary(Longi)
summary(Surv)


M1 <- joint(formLong = serBilir ~ year + drug + (1|id),
      dataLong = Longi, id="id", timeVar = "year", family="lognormal")
summary(M1)


M2 <- joint(formLong = list(serBilir ~ year * drug + sex + (1+year|id),
                            platelets ~ year * sex + drug + (1+year|id),
                            spiders ~ year * drug + sex + (1+year|id)),
            dataLong = Longi, id="id", timeVar="year", corLong = T,
            family=c("lognormal", "poisson", "binomial"),
            control=list(int.strategy="eb"))
summary(M2)
summary(M2, sdcor=T)


M2 <- joint(formLong = list(serBilir ~ year * drug + sex + (1+year|id),
                            platelets ~ year * sex + drug + (1+year|id),
                            spiders ~ year * drug + sex + (1+year|id)),
            dataLong = Longi, id="id", timeVar="year", corLong = T,
            family=c("lognormal", "poisson", "binomial"),
            link=c("default", "default", "probit"),
            control=list(int.strategy="eb"))
summary(M2)


DTH <- inla.surv(time = Surv$years, event = Surv$death)
f1 <- function(x) x^2
f2 <- function(x) x^3


M3 <- joint(formSurv = DTH ~ drug,
            formLong = serBilir ~ (1 + year + f1(year) + f2(year))*drug +
              (1 + year + f1(year) + f2(year) | id), family="lognormal",
            dataLong = Longi, id="id", timeVar="year", basRisk = "rw1", NbasRisk = 25,
            control=list(int.strategy="eb"), assoc="CV_CS")
summary(M3)

# "" (empty string) = no association
# SRE = individual deviation at time t defined by random effects
# SRE_ind = each random effect associated to a parameter in survival
# CV = current value of the linear predictor (error free)
# CS = current slope
# CV_CS = current value and current slope


plotM3 <- plot(M3, sdcor=T)
plotM3$Outcomes$L1
plotM3$Outcomes$S1
plotM3$Covariances
plotM3$Random + scale_y_log10() + 
  geom_abline(slope=0, intercept=0, linetype="dashed")


M4 <- joint(formSurv = DTH ~ sex + drug,
            formLong = albumin ~ year * drug + (1 + year|id),
            dataLong = Longi, id="id", timeVar="year", assoc="CV",
            control=list(priorFixed=list(mean.intercept=0, prec.intercept=0.16,
                                         mean=0, prec=0.16),
                         priorAssoc=list(mean=0, prec=0.16)))
summary(M4)


# JMBayes
library(JMbayes)
M4JMB_lme <- lme(albumin ~ (1 + year)*drug,
                 random = ~ 1 + year |id, data = Longi)
M4JMB_cox <- coxph(Surv(Surv$years, Surv$death) ~ sex + drug,
                   data = Surv, x = TRUE)
JMpr = list(priorMean.alphas=0, priorTau.alphas = matrix(0.16))
M4JMB <- jointModelBayes(M4JMB_lme, M4JMB_cox, timeVar = "year", priors=JMpr)
# Computation time in the table includes LME + Cox + JM
Summary(M4JMB)

# rstanarm
library(rstanarm)
library(survival)
options(mc.cores = parallel::detectCores())
M4rstanarm <- stan_jm(
  formulaLong = list(albumin ~ (1 + year)*drug + (1 + year |id)),
  formulaEvent = Surv(years, death) ~ sex + drug,
  dataLong = Longi, dataEvent = Surv,
  time_var = "year",
  priorLong_intercept = normal(0, 2.5, autoscale=TRUE),
  priorLong = normal(0, 2.5),
  priorEvent_assoc = normal(0, 2.5),
  seed = 12345)
summary(M4rstanarm)


TSP <- inla.surv(time=Surv$years, event=Surv$trans)

M5 <- joint(formLong = serBilir ~ year * (drug + sex) + (1 + year|id),
            formSurv = list(DTH ~ sex + drug,
                            TSP ~ edema * sex), dataLong = Longi,
            id="id", timeVar="year", family="lognormal", basRisk=c("rw1", "rw1"),
            assoc=c("SRE", "SRE_ind"), control=list(int.strategy="eb"))

summary(M5)


M6 <- joint(formLong = list(serBilir ~ year * drug  + sex + (1 + year|id),
                            platelets ~ year * sex + drug + (1 + year|id),
                            spiders ~ year * drug + (1 + year|id)),
            formSurv = list(DTH ~ drug,
                            TSP ~ drug), dataLong = Longi, id="id", timeVar="year",
            assoc = list(c("CV", "CV"), c("SRE", ""), c("CV_CS", "CS")),
            family=c("lognormal", "poisson", "binomial"), basRisk=c("rw2", "rw2"),
            control=list(int.strategy="eb"))
summary(M6)


data(SurvMS)
H12 <- inla.surv(time=SurvMS[[1]]$Tstop, event=SurvMS[[1]]$status)
H13 <- inla.surv(time=SurvMS[[2]]$Tstop, event=SurvMS[[2]]$status)
H23 <- inla.surv(time=SurvMS[[3]]$Tstop, truncation = SurvMS[[3]]$Tstart, event=SurvMS[[3]]$status)

M7 <- joint(formSurv = list(H12 ~ X, H13 ~ X, H23 ~ X),
            basRisk = c("rw2", "rw2", "rw2"), dataSurv = SurvMS)
summary(M7)


data(LongMS)

M8 <- joint(formSurv = list(H12 ~ X, H13 ~ X, H23 ~ X),
            formLong = y ~ time + X + (1+time|id), timeVar = "time", id="id",
            assoc = c("CV", "CV", "CV"), dataLong = LongMS,
            basRisk = c("rw2", "rw2", "rw2"), dataSurv = SurvMS)
summary(M8)



Nsplines <- ns(Longi$year, knots = c(1, 4))
f1 <- function(x) predict(Nsplines, x)[,1]
f2 <- function(x) predict(Nsplines, x)[,2]
f3 <- function(x) predict(Nsplines, x)[,3]

M9 <- joint(formLong = list(serBilir ~ (1 + f1(year) + f2(year) + f3(year))*drug + (1 + f1(year) + f2(year) + f3(year)|id),
                            SGOT ~ (1 + f1(year) + f2(year) + f3(year))*drug + (1 + f1(year) + f2(year) + f3(year)|id),
                            albumin ~ (1 + year) * drug + (1+year|id),
                            platelets ~ (1 + f1(year) + f2(year) + f3(year))*drug + (1 + f1(year) + f2(year) + f3(year)|id),
                            spiders ~ (1 + year) * drug + (1+year|id)),
            formSurv = list(DTH ~ drug, TSP ~ drug),
            dataLong = Longi, id="id", timeVar="year", basRisk = c("rw2", "rw1"),
            family=c("lognormal", "lognormal", "gaussian", "poisson", "binomial"),
            assoc=list(c("CV_CS", "CV"), c("CV", ""), c("CV", "CV"), c("CV", "CV"), c("CV", "")),
            control=list(priorFixed=list(mean=0, prec=1, mean.intercept=0, prec.intercept=1),
                         priorAssoc = list(mean=0, prec=1),
                         priorRandom = list(r=10, R=1), int.strategy="eb"))
summary(M9)



names(inla.models()$link)
inla.models()$likelihood

inla.doc("loglog")






