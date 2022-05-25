### example from
### https://becarioprecario.bitbucket.io/spde-gitbook/ch-risk.html

library(INLA)

## load the data 
data(Leuk)

### Survival time as year
Leuk$time <- Leuk$time / 365
round(sapply(Leuk[, c(1, 2, 5:8)], summary), 2)

### KM
library(survival)
km <- survfit(Surv(time, cens) ~ sex, Leuk) 

par(mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.5, 0.5, 0), las = 1)
plot(km, conf.int = TRUE, col = 2:1) 
legend('topright', c('female', 'male'), lty = 1, col = 2:1,
  bty = "n") 

### mesh 
loc <- cbind(Leuk$xcoord, Leuk$ycoord)
nwseg <- inla.sp2segment(nwEngland)

### non-convex hull boundaries around the polygon
bnd1 <- inla.nonconvex.hull(nwseg$loc, 0.03, 0.1, resol = 50)
bnd2 <- inla.nonconvex.hull(nwseg$loc, 0.25)

### make the mesh
mesh <- inla.mesh.2d(
    loc, boundary = list(bnd1, bnd2),
    max.edge = c(0.05, 0.2), cutoff = 0.02)

### visualize it 
plot(mesh, asp=1)
lines(nwseg$loc, col=2)
points(Leuk$x, Leuk$y, pch=8)

#### ANOTHER MESH
### intitial mesh inside the first boundary segment
mesh0 <- inla.mesh.2d(
    boundary=bnd1, 
    max.edge=0.03)
mesh0$n

plot(nwEngland)
plot(mesh0, add=TRUE)

### now, consider the locations of this inital ponts 
mesh <- inla.mesh.2d(
    loc=mesh0$loc[,1:2],
    max.edge=0.15, offset=0.25, cutoff=0.02)
mesh$n

### visualize it 
plot(mesh, asp=1)
plot(nwEngland, add=TRUE)
points(Leuk$x, Leuk$y, pch=8)


## projector matrix 
A <- inla.spde.make.A(mesh, loc)

## the spde model object (with priors) 
spde <- inla.spde2.pcmatern(mesh = mesh,
  prior.range = c(0.05, 0.01), # P(range < 0.05) = 0.01
  prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01


## model formulae (without spatial)
form0 <- inla.surv(time, cens) ~
    0 + a0 + sex + age + wbc + tpi 


### model formulae with the spatial effect 
form <- update(
    form0, . ~ . +
               f(spatial, model = spde))


### stack the data 
stk <- inla.stack(
  data = list(time = Leuk$time, cens = Leuk$cens), 
  A = list(A, 1), 
  effect = list(
    list(spatial = 1:spde$n.spde), 
    data.frame(a0 = 1, Leuk[, -c(1:4)]))) 


### inla result 
r <- inla(
    form, family = "weibullsurv",
    data = inla.stack.data(stk), 
    control.predictor = list(
        A = inla.stack.A(stk),
        compute = TRUE)) 


### fixed effects result 
round(r$summary.fixed, 4)

### hyperparameters summary result 
round(r$summary.hyperpar, 4)

### projector to a grid 
bbnw <- bbox(nwEngland)
bbnw
r0 <- diff(range(bbnw[1, ])) / diff(range(bbnw[2, ]))
r0
prj <- inla.mesh.projector(mesh, xlim = bbnw[1, ], 
  ylim = bbnw[2, ], dims = round(c(r0, 1)*300))

## NA's were not to plot 
spat.m <- inla.mesh.project(prj, r$summary.random$spatial$mean)
spat.sd <- inla.mesh.project(prj, r$summary.random$spatial$sd)
ov <- over(SpatialPoints(prj$lattice$loc), nwEngland)
spat.sd[is.na(ov)] <- NA
spat.m[is.na(ov)] <- NA

### plot the spatial risk
library(fields)
par(mfrow = c(1, 2), mar = c(0, 0, 0, 5))
plot(nwEngland)
image.plot(x = prj$x, y = prj$y, z = spat.m, add=TRUE)
plot(nwEngland)
image.plot(x = prj$x, y = prj$y, z = spat.sd, add = TRUE)

##################################################################################
## CPH 
m0 <- coxph(Surv(time, cens) ~ sex + age + wbc + tpi, Leuk)

## expand the data to fit the CPH 
form0
cph.leuk <- inla.coxph(form0,
  data = data.frame(a0 = 1, Leuk[, 1:8]),
  control.hazard = list(n.intervals = 25))


## fit
cph.res0 <- inla(form0, family = 'coxph', 
  data = data.frame(a0 = 1, Leuk[, c(1,2, 5:8)])) 


## add the spatial effect into the model formulae 
cph.formula <- update(cph.leuk$formula, 
  '. ~ . + f(spatial, model = spde)')


## The projector for the expanded CPH dataset 
cph.A <- inla.spde.make.A(mesh,
  loc = cbind(cph.leuk$data$xcoord, cph.leuk$data$ycoord))


## stack the CPH expanded dataset 
cph.stk <- inla.stack(
  data = c(list(E = cph.leuk$E), cph.leuk$data[c('y..coxph')]),
  A = list(cph.A, 1),
  effects = list(
    list(spatial = 1:spde$n.spde), 
      cph.leuk$data[c('baseline.hazard', 'a0', 
        'age', 'sex', 'wbc', 'tpi')]))

cph.data <- c(inla.stack.data(cph.stk), cph.leuk$data.list)


## CPH fit
cph.res <- inla(cph.formula, family = 'Poisson', 
  data = cph.data, E = cph.data$E, 
  control.predictor = list(A = inla.stack.A(cph.stk)))


## compare fitted fixed effects 
round(data.frame(surv = coef(summary(m0))[, c(1,3)], 
  r0 = cph.res0$summary.fixed[-1, 1:2], 
  r1 = cph.res$summary.fixed[-1, 1:2]), 4) 


## similarity of the fitted random effect 
s.m <- inla.mesh.project(prj, cph.res$summary.random$spatial$mean)
cor(as.vector(spat.m),  as.vector(s.m), use = 'p')

### as well the error 
s.sd <- inla.mesh.project(prj, cph.res$summary.random$spatial$sd)
cor(log(as.vector(spat.sd)), log(as.vector(s.sd)), use = 'p')

