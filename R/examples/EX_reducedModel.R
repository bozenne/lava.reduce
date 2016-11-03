#### regression ####

## simulation
m <- lvm()
m <- regression(m,y='y1',x='x'%++%1:2)
m <- regression(m,y='y1',x='z'%++%1:5)

set.seed(10)
d <- sim(m,150)

## reduced model 1
mR1 <- lvm.reduced()
mR1 <- regression(mR1,y='y1',x='x'%++%1:2)
mR1 <- regression(mR1,y='y1',x='z'%++%1:5, reduce = TRUE)

## reduced model 2
mR2 <- reduce(m)

## check estimation
emGS <- estimate(m, d)
em1 <- estimate(m, d, estimator = "gaussian1")
emR1 <- estimate(mR1, d)
coef(emGS) - coef(em1)
coef(em1) - coef(emR1)[names(coef(em1))]

emR2 <- estimate(mR1, d)

#### latent variable model ####
m <- lvm()
m <- regression(m,y=c('y1','y2','y3','y4'),x='eta')
m <- regression(m,y=c('y2','y3'),x='x'%++%1:5)
latent(m) <- ~eta
m <- regression(m,y=c('y1','y2'),x='z'%++%1:2)
covariance(m) <- y2~y1

# simul
set.seed(10)
d <- sim(m,100)

# reduced model
mR1 <- reduce(m, endo = c("y2"))
mR2 <- reduce(m)

## estimation
em <- estimate(m,d)
emR <- estimate(mR1, data = d)
coef(em)
coef(emR)[names(coef(em))]
