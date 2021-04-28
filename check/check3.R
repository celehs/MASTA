library(MASTA)
#head(longitudinal)
#head(follow_up_time)
#head(survival)

#---------------------------------------------------
# Error 3:
#---------------------------------------------------
# Error in xj[i] : invalid subscript type 'list'
#In addition: Warning message:
#  In ft.e[pos, 2] * tmp$densities[locm, 1] :
#  
#  Error in xj[i] : invalid subscript type 'list' 
#---------------------------------------------------

#---- create data sample for test ---
nn = 200 
ncode = 4
fu_month = 40
fu_rate = 0.1

set.seed(123)
#-- follow-up-data
id = 1:nn
train_valid = rep(2, nn) ; train_valid[1:(nn/2)] = 1
fu_time = rexp(nn, rate=fu_rate)
fu_time = pmax(fu_time, fu_month)  ##pmax not pmin!
follow_up_time_data = data.frame(id = id, fu_time = fu_time,  train_valid = train_valid)

#-- longitudinal data ---
out=c()
junk = 1:fu_month
for (j in 1:ncode){
  for (i in 1:nn){
  rsize = rpois(1, lambda=20)
  tmp1=sort(sample(junk, size = rsize, replace=TRUE))
  tmp2=cbind(rep(j,rsize),rep(i,rsize),tmp1)
  out=rbind(out, tmp2)
  }
}
out = data.frame(out)
colnames(out)=c("code","id","time")
longitudinal_data = out
head(longitudinal_data)

system.time(Z <- fpca.combine(longitudinal_data, follow_up_time_data, K.select = "PropVar"))
  
