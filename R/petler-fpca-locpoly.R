### local linear regression
den_locpoly<-function(partition,t,N){
  n           = length(N)
  l           = length(partition)
  partition_x = {partition[-1]+partition[-l]}/2
  cumsumN     = c(0,cumsum(N))
  f_locpoly   = lapply(seq(1,n),function(i){
    a         = {cumsumN[i]+1}:cumsumN[i+1]
    ### the histogram density estimate for the ith subject
    f_H       = hist(t[a],breaks=partition,plot=FALSE)$counts/N[i]/diff(partition)
    ### leave one out cross validation, local linear regression
    GCV_loclinear = function(hh){
      D   = sapply(partition_x,function(x) partition_x-x)
      W   = dnorm(D/hh)
      S0  = apply(W,1,sum)
      S1  = apply(W*D,1,sum)
      S2  = apply(W*{D^2},1,sum)
      L   = W*{t(VTM(S2,l-1))-D*t(VTM(S1,l-1))}/t(VTM(S0*S2-S1^2,l-1))
      f_oneout = L%*%f_H
      return(sum({{f_H-f_oneout}/{1-mean(diag(L))}}^2))
    }
    ### find bandwidth that minimize generalized cross validation
    # here set lower limit to be the min(diff(partition_x)) to avoid error
    hh   = optimize(GCV_loclinear,lower=min(diff(partition_x)),upper=IQR(t)*2)$minimum
    ### use the obtained GCV bandwidth to get local linear regression estimates
    D   = t(sapply(t[a],function(x) partition_x-x))
    W   = dnorm(D/hh)
    S0  = apply(W,1,sum)
    S1  = apply(W*D,1,sum)
    S2  = apply(W*{D^2},1,sum)
    L   = W*{t(VTM(S2,l-1))-D*t(VTM(S1,l-1))}/t(VTM(S0*S2-S1^2,l-1))
    return(L%*%f_H)
  })
}

### local linear regression
den_locpoly2<-function(t,N){
  n           = length(N)
  nbreak      = numeric(n)
  cumsumN     = c(0,cumsum(N))
  f_locpoly   = lapply(seq(1,n),function(i){
    a         = {cumsumN[i]+1}:cumsumN[i+1]
    ### the histogram density estimate for the ith subject
    f_H       = hist(t[a],plot=FALSE)
    partition = f_H$breaks
    nbreak[i] = length(partition)
    f_H       = f_H$counts/N[i]/diff(partition)
    partition_x = {partition[-1]+partition[-length(partition)]}/2
    ### leave one out cross validation, local linear regression
    GCV_loclinear = function(hh){
      D   = sapply(partition_x,function(x) partition_x-x)
      W   = dnorm(D/hh)
      S0  = apply(W,1,sum)
      S1  = apply(W*D,1,sum)
      S2  = apply(W*{D^2},1,sum)
      L   = W*{t(VTM(S2,nbreak[i]-1))-D*t(VTM(S1,nbreak[i]-1))}/t(VTM(S0*S2-S1^2,nbreak[i]-1))
      f_oneout = L%*%f_H
      return(sum({{f_H-f_oneout}/{1-mean(diag(L))}}^2))
    }
    ### find bandwidth that minimize generalized cross validation
    # here set lower limit to be the min(diff(partition_x)) to avoid error
    hh   = optimize(GCV_loclinear,lower=min(diff(partition_x)),upper=IQR(t[a])*2)$minimum
    ### use the obtained GCV bandwidth to get local linear regression estimates
    D   = t(sapply(t[a],function(x) partition_x-x))
    W   = dnorm(D/hh)
    S0  = apply(W,1,sum)
    S1  = apply(W*D,1,sum)
    S2  = apply(W*{D^2},1,sum)
    L   = W*{t(VTM(S2,nbreak[i]-1))-D*t(VTM(S1,nbreak[i]-1))}/t(VTM(S0*S2-S1^2,nbreak[i]-1))
    return(L%*%f_H)
  })
}
