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










