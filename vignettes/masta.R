## ---- echo=TRUE---------------------------------------------------------------
library(MASTA)

## -----------------------------------------------------------------------------
head(longitudinal)
table(longitudinal$code)

head(follow_up_time) ; 
nrow(follow_up_time) ;

head(survival)
nrow(survival)

## -----------------------------------------------------------------------------
system.time(Z <- fpca.combine(longitudinal, follow_up_time, K.select = "PropVar"))

## -----------------------------------------------------------------------------
system.time(b <- masta.fit(Z, survival, follow_up_time, Tend=1, cov_group = NULL, thresh = 0.7, PCAthresh = 0.9, seed = 100))
names(b)

