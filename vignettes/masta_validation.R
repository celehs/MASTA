## -----------------------------------------------------------------------------
library(MASTA)
fpca <- fpca.combine(longitudinal, follow_up_time, K.select = "PropVar")
ft <- masta.fit(fpca, survival, follow_up_time, Tend=1, cov_group = NULL, thresh = 0.7, PCAthresh = 0.9, seed = 100)

## -----------------------------------------------------------------------------
head(new_longitudinal)

head(new_follow_up_time) ; 

head(new_survival)

## -----------------------------------------------------------------------------
val=masta_validation(ft, new_longitudinal, new_follow_up_time, new_survival)
names(val)

