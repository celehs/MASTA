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

