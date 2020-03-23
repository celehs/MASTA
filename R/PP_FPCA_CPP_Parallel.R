FPC_Kern_S <- function(x, t, N, h1, h2) {
  grp <- rep(1:length(N), N)
  M <- outer(t, x, "-")
  D1 <- dnorm(M, 0, h1)
  D2 <- dnorm(M, 0, h2)   
  S2 <- rowsum(D2, grp)
  list(f_mu = colSums(D1), 
       Cov_G = crossprod(S2) - crossprod(D2))
}

######################################################################
### Second version of the functional principal component analysis 
### on rare events (Wu et al.,	2013). 

########################################################################
VTM<-function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}

########################################################################
## FPCA approach by Wu et al (2013)
## n = num of patients
## t: observed event times of all the individuals, can have duplicates
## h: bandwidth
## N: vector, contains num of observed event of each patient
## index: index of the patient; default is NULL, i.e., patient are labelled from 1 to n.
## ngrid: number of grid points in estimating covariance function g
## subdivisions: number of subdivisions used in GCV
## propvar: proportion of variation used to select number of FPCs
## PPIC: if TRUE, the propvar will be ignored and PPIC will be used to select number of FPCS
## polybinx: if use the same partition (x) for the polynomial regression
########################################################################
############                    Parallel                    ############     
########################################################################
PP_FPCA_CPP_Parallel<-function(t, h1 = NULL, h2 = NULL, N, bw = "ucv", Tend = 1, # assume it is transformed to [0,1]
                      ngrid = 101, K.select = c('PropVar','PPIC'), Kmax = 10,
                      propvar = 0.9, ## For K.select == PropVar
                      density.method = c("kernel","local linear"), ## PPIC
                      polybinx = FALSE, derivatives = FALSE,
                      nsubs = 10, subsize = NULL,PPIC.sub = TRUE){
  ## eleminate patients with 0 observed event
  if(sum(N==0)>0){
    NN      = N
    N.index = N!=0
    N       = N[N.index]
    # cat("Note: patients with zero observed event have been eliminated from analysis!","\n")
  }
  n      = length(N) # number of patients with at least one observed event
  ## if h is null then set it to be the optimal bandwidth under GCV
  if(is.null(h1) & is.null(h2)){
    if(bw == "nrd0"){
      h1 = h2 = bw.nrd0(t)
    } else if(bw == "nrd"){
      h1 = h2 = bw.nrd(t)
    } else if(bw == "ucv"){ # leave one out cross validation
      h1 = h2 = bw.ucv(t)
    } else if(bw == "bcv"){ # biased cross validation
      h1 = h2 = bw.bcv(t)
    } else if(bw == "SJ-ste"){
      h1 = h2 = bw.SJ(t,method=c("ste"))
    } else if(bw == "SJ-dpi"){
      h1 = h2 = bw.SJ(t,method=c("dpi"))
    } else {
      h1 = h2 = bw.ucv(t)
    }
  }
  ## get a fine grid (x,y)
  # tmp  = range(t)
  x    = y = seq(0,Tend,length.out = ngrid)
  ## estimate the mean density f_mu and the intermediate value g that is used to 
  if(is.null(subsize)){
    subsize = floor(n/nsubs)
    subsize = c(rep(subsize,nsubs-1),n-subsize*(nsubs-1))
  }
  cumsumsub = cumsum(c(0,subsize))
  cumsumN   = c(0,cumsum(sapply(1:nsubs,function(i) sum(N[{cumsumsub[i]+1}:cumsumsub[i+1]]))))
  tmplist  = foreach(i=1:nsubs) %dopar% {
    tmp    = FPC_Kern_S(x,t[{cumsumN[i]+1}:cumsumN[i+1]],N[{cumsumsub[i]+1}:cumsumsub[i+1]],h1,h2)
    list(f_mu = as.numeric(tmp$f_mu)/sum(N), g=tmp$Cov_G/sum(N*(N-1)))
  }
  f_mu     = apply(simplify2array(lapply(tmplist,`[[`,1)),1,sum)
  G        = apply(simplify2array(lapply(tmplist,`[[`,2)),c(1,2),sum)
  G        = G-outer(f_mu,f_mu)
  
  G.eigen  = svd(G)
  delta    = sqrt(x[2]-x[1])
  baseline = as.numeric(t(G.eigen$u[,1:Kmax])%*%f_mu*delta) # second term in xi
  cumsumN2 = c(0,cumsum(N))
  
  # interpolate the eigenvectors to get eigenfunctions and then get the xi's
  xi       = foreach(i=1:nsubs,.combine = cbind) %dopar% {
    tmp2   = apply(G.eigen$u[,1:Kmax],2,function(s) approx(x = x, y = s, xout = t[{cumsumN[i]+1}:cumsumN[i+1]])$y)/delta 
    indexi = rep({1+cumsumsub[i]}:cumsumsub[i+1],N[{cumsumsub[i]+1}:cumsumsub[i+1]])
    -baseline+t(apply(tmp2,2,FUN=function(x) tapply(x,indexi,mean)))
  }
  attr(xi,"dimnames") = NULL
  
  ## select number of FPCs
  K.select = match.arg(K.select)
  if(K.select=="PropVar"){
    # method 1: proportion of variation >= 90%
    K       = cumsum(G.eigen$d)/sum(G.eigen$d)
    K       = min(which(K>=propvar))
    if(K>Kmax) K = Kmax
  } else if(K.select=="PPIC"){
    K       = Kmax
    # K         = sapply(seq(1,ngrid),function(i) sum(G.eigen$d[1:i]))/sum(G.eigen$d)
    # # cat("Prop of variation=",K,"\n")
    # K         = min(which(K>=0.99))
    if (missing(density.method)) density.method = "kernel" else density.method = match.arg(density.method)
    
    if(density.method=="local linear"){
      if(polybinx){
        f_locpoly = foreach(i=1:nsubs,.combine = c)%dopar%{
          den_locpoly(partition=x,t=t[{cumsumN[i]+1}:cumsumN[i+1]],N=N[{cumsumsub[i]+1}:cumsumsub[i+1]])
        }
      } else{
        f_locpoly = foreach(i=1:nsubs,.combine = c)%dopar%{
          den_locpoly2(t=t[{cumsumN[i]+1}:cumsumN[i+1]],N=N[{cumsumsub[i]+1}:cumsumsub[i+1]])
        }
      } 
    } else{
      f_locpoly = foreach(i=1:nsubs,.combine = c)%dopar%{
        sapply({cumsumsub[i]+1}:cumsumsub[i+1],function(j){
          tmp = density(t[{cumsumN2[j]+1}:cumsumN2[j+1]],bw="nrd")
          tmp$y[sapply(t[{cumsumN2[j]+1}:cumsumN2[j+1]],function(s) which.min(abs(s-tmp$x)))]
        })
      }
    }
    K         = foreach(i=1:nsubs,.combine=rbind)%dopar%{
      sapply(c(1:K),function(k) PPIC(K = k,f_locpoly = f_locpoly[{cumsumsub[i]+1}:cumsumsub[i+1]],
                                     t=t[{cumsumN[i]+1}:cumsumN[i+1]],N=N[{cumsumsub[i]+1}:cumsumsub[i+1]],
                                     f_mu=f_mu,G.eigen_v = G.eigen$u,xi = xi[,{cumsumsub[i]+1}:cumsumsub[i+1]],
                                     xgrid = x,delta = delta))
    } ## parallel different from non-parallel results
    K         = apply(K,2,sum)
    K         = which.min(K)
  }
  
  ## get density functions
  if(K==1){
    tmp     = foreach(i=1:nsubs,.combine=cbind)%dopar%{
      rawden = f_mu + outer(G.eigen$u[,1],xi[1,{cumsumsub[i]+1}:cumsumsub[i+1]])/delta # density functions
      ## make adjustments to get valid density functions (nonnegative+integrate to 1)
      apply(rawden,2,function(x){
        x[x<0] = 0     # non-negative
        x      = x/sum(x) # integrate to delta^2
        return(x)
      })/delta^2
    }
    scores  = data.frame(xi[1:max(K,Kmax),])
  } else{
    tmp     = foreach(i=1:nsubs,.combine=cbind)%dopar%{
      rawden = f_mu + {G.eigen$u[,1:K]/delta}%*%xi[1:K,{cumsumsub[i]+1}:cumsumsub[i+1]] # density functions
      # make adjustments to get valid density functions (nonnegative+integrate to 1)
      apply(rawden,2,function(x){
        x[x<0] = 0     # non-negative
        x      = x/sum(x) # integrate to delta^2
        return(x)
      })/delta^2
    }
    scores  = data.frame(t(xi[1:max(K,Kmax),]))
  }
  
  names(scores) = paste0("score_",seq(1:max(K,Kmax)))  # FPC scores
  basis   = data.frame(cbind(x,G.eigen$u[,1:max(K,Kmax)]/delta))
  names(basis) = c("x",paste0("basis_",seq(1:max(K,Kmax)))) # FPC eigenfunctions
  
  ## name the density functions 
  if(exists("N.index")){
    colnames(tmp)  = paste0("f",which(N.index))
    scores         = data.frame(id = as.character(which(N.index)), scores, 
                                stringsAsFactors =FALSE) # add id to scores
  } else{
    colnames(tmp)  = paste0("f",seq(1,n))
    scores         = data.frame(id = as.character(seq(1,n)), scores, 
                                stringsAsFactors = FALSE)
  }
  tmp     = as.data.frame(tmp)
  tmp     = data.frame(x=x,tmp)
  
  # ## return K and prop of var by now
  # tmp     = data.frame(tmp,info=c(K,sum(G.eigen$d[1:K])/sum(G.eigen$d),rep(0,ngrid-2)))
  ## get derivatives if derivatives = TRUE
  if(derivatives){
    
    tmp2  = foreach(i=1:nsubs,.combine = cbind)%dopar%{
      {tmp[-c(1:2),{cumsumsub[i]+1}:cumsumsub[i+1]+1]-tmp[-c(1:2+ngrid-2),{cumsumsub[i]+1}:cumsumsub[i+1]+1]}/diff(x,lag=2) 
    }
    
    tmp2  = as.data.frame(tmp2)
    tmp2  = data.frame(x=x[-c(1,ngrid)],tmp2)
    colnames(tmp2) = paste0("d",colnames(tmp))
    return(list(scores = scores, densities = tmp, derivatives = tmp2, 
                mean = data.frame(x = x, f_mu = f_mu),cov = G,
                basis = basis, baseline = baseline,K = K))
  }
  else{
    return(list(scores = scores, densities = tmp, 
                mean = data.frame(x = x, f_mu = f_mu), 
                cov  = G,
                basis = basis, baseline = baseline,K = K))
  }
}


########################################################################
############                  Non-Parallel                  ############     
########################################################################
PP_FPCA_CPP<-function(t, h1 = NULL, h2 = NULL, N, bw = "ucv", Tend = 1, # assume it is transformed to [0,1]
                      ngrid = 101, K.select = c('PropVar','PPIC'), Kmax = 10,
                      propvar = 0.9, ## For K.select == PropVar
                      density.method = c("kernel","local linear"), ## PPIC
                      polybinx = FALSE, derivatives = FALSE){
  ## eleminate patients with 0 observed event
  if(sum(N==0)>0){
    NN      = N
    N.index = N!=0
    N       = N[N.index]
    # cat("Note: patients with zero observed event have been eliminated from analysis!","\n")
  }
  n      = length(N) # number of patients with at least one observed event
  ## if h is null then set it to be the optimal bandwidth under GCV
  if(is.null(h1) & is.null(h2)){
    if(bw == "nrd0"){
      h1 = h2 = bw.nrd0(t)
    } else if(bw == "nrd"){
      h1 = h2 = bw.nrd(t)
    } else if(bw == "ucv"){ # leave one out cross validation
      h1 = h2 = bw.ucv(t)
    } else if(bw == "bcv"){ # biased cross validation
      h1 = h2 = bw.bcv(t)
    } else if(bw == "SJ-ste"){
      h1 = h2 = bw.SJ(t,method=c("ste"))
    } else if(bw == "SJ-dpi"){
      h1 = h2 = bw.SJ(t,method=c("dpi"))
    } else {
      h1 = h2 = bw.ucv(t)
    }
  }
  ## get a fine grid (x,y)
  # tmp  = range(t)
  x    = y = seq(0,Tend,length.out = ngrid)
  ## estimate the mean density f_mu and the intermediate value g that is used to 
  
  tmp      = FPC_Kern_S(x,t,N,h1,h2)
  f_mu     = as.numeric(tmp$f_mu)/sum(N)
  G        = tmp$Cov_G/sum(N*(N-1))-outer(f_mu,f_mu)
  
  G.eigen  = svd(G)
  delta    = sqrt(x[2]-x[1])
  baseline = as.numeric(t(G.eigen$u[,1:Kmax])%*%f_mu*delta) # second term in xi
  cumsumN2 = c(0,cumsum(N))
  
  # interpolate the eigenvectors to get eigenfunctions
  tmp2    = apply(G.eigen$u[,1:Kmax],2,function(s) approx(x = x, y = s, xout = t)$y)/delta ### check whether we need parallel this
  indx    = rep(1:n,N)
  xi      = -baseline + t(apply(tmp2,2,FUN=function(x) tapply(x,indx,mean)))# FPC scores, ith column for ith patient
  rm(indx)
  
  ## select number of FPCs
  K.select = match.arg(K.select)
  if(K.select=="PropVar"){
    # method 1: proportion of variation >= 90%
    K       = cumsum(G.eigen$d)/sum(G.eigen$d)
    cat("Propvar=",K,"\n")
    K       = min(which(K>=propvar))
    if(K>Kmax) K = Kmax
  } else if(K.select=="PPIC"){
    K       = Kmax
    if (missing(density.method)) density.method = "kernel" else density.method = match.arg(density.method)
    
    if(density.method=="local linear"){
      if(polybinx){
        f_locpoly = den_locpoly(partition=x,t=t,N=N)
      } else{
        f_locpoly = den_locpoly2(t=t,N=N)
      }
    } else{
      f_locpoly = sapply(1:n,function(j){
        tmp = density(t[{cumsumN2[j]+1}:cumsumN2[j+1]],bw="nrd")
        tmp$y[sapply(t[{cumsumN2[j]+1}:cumsumN2[j+1]],function(s) which.min(abs(s-tmp$x)))]
      })
    }
    K         = sapply(1:K,function(k) PPIC(K = k,f_locpoly = f_locpoly,t=t,N=N,f_mu=f_mu,G.eigen_v = G.eigen$u,xi = xi,xgrid = x,delta = delta))
    K         = which.min(K)
  }
  
  ## get density functions 
  if(K==1){
    tmp     = f_mu + outer(G.eigen$u[,1],xi[1,])/delta # density functions
    scores  = data.frame(xi[1:max(K,Kmax),])
  } else{
    tmp     = f_mu + {G.eigen$u[,1:K]/delta}%*%xi[1:K,] # density functions
    scores  = data.frame(t(xi[1:max(K,Kmax),]))
  }
  ## make adjustments to get valid density functions (nonnegative+integrate to 1)
  tmp     = apply(tmp,2,function(x){
    x[x<0] = 0     # non-negative
    x      = x/sum(x) # integrate to delta^2
    return(x)
  })
  tmp     = tmp/delta^2
  
  
  names(scores) = paste0("score_",seq(1:max(K,Kmax)))  # FPC scores
  basis   = data.frame(cbind(x,G.eigen$u[,1:max(K,Kmax)]/delta))
  names(basis) = c("x",paste0("basis_",seq(1:max(K,Kmax)))) # FPC eigenfunctions
  
  ## name the density functions 
  if(exists("N.index")){
    colnames(tmp)  = paste0("f",which(N.index))
    scores         = data.frame(id = as.character(which(N.index)), scores, 
                                stringsAsFactors =FALSE) # add id to scores
  } else{
    colnames(tmp)  = paste0("f",seq(1,n))
    scores         = data.frame(id = as.character(seq(1,n)), scores, 
                                stringsAsFactors = FALSE)
  }
  tmp     = as.data.frame(tmp)
  tmp     = data.frame(x=x,tmp)
  
  # ## return K and prop of var by now
  # tmp     = data.frame(tmp,info=c(K,sum(G.eigen$d[1:K])/sum(G.eigen$d),rep(0,ngrid-2)))
  ## get derivatives if derivatives = TRUE
  if(derivatives){
    tmp2  = {tmp[-c(1:2),-1]-tmp[-c(1:2+ngrid-2),-1]}/diff(x,lag=2) 
    tmp2  = as.data.frame(tmp2)
    tmp2  = data.frame(x=x[-c(1,ngrid)],tmp2)
    colnames(tmp2) = paste0("d",colnames(tmp))
    return(list(scores = scores, densities = tmp, derivatives = tmp2, 
                mean = data.frame(x = x, f_mu = f_mu),cov = G,
                basis = basis, baseline = baseline,K = K))
  }
  else{
    return(list(scores = scores, densities = tmp, 
                mean = data.frame(x = x, f_mu = f_mu), 
                cov = G,
                basis = basis, baseline = baseline,K = K))
  }
}

### Pseudo-poisson informtion criterion (PPIC) for selecting the
### number of FPCs
PPIC<-function(K,f_locpoly,t,N,f_mu,G.eigen_v,xi,xgrid,delta){
  if(K==1){
    tmp     = f_mu + outer(G.eigen_v[,1],xi[1,])/delta # density functions
    # tmp2    = f_mu_d + outer(tmp2[,1],xi[1,])/delta # derivatives of density functions
  } else{
    tmp     = f_mu + {G.eigen_v[,1:K]/delta}%*%xi[1:K,] # density functions
    # tmp2    = f_mu_d + {tmp2[,1:K]/delta}%*%xi[1:K,] # derivatives of density functions
  }
  ## make adjustments to get valid density functions (nonnegative+integrate to 1)
  tmp     = apply(tmp,2,function(x){
    x[x<0] = 0     # non-negative
    x      = x/sum(x) # integrate to delta^2
    return(x)
  })
  tmp     = tmp/delta^2
  n       = length(N)
  tmp     = lapply(seq(1,n),function(i) approx(xgrid,tmp[,i],t[(sum(N[1:i])-N[i]+1):sum(N[1:i])])$y)
  f_locpoly = unlist(f_locpoly)
  if(sum(f_locpoly<=0)>0) f_locpoly[f_locpoly<1e-8] = 1e-8
  tmp       = unlist(tmp)
  tmp[tmp<1e-8] = 1e-8
  D_K       = 2*sum(f_locpoly*log(f_locpoly/tmp)-(f_locpoly-tmp))+2*K
  return(D_K)
}

### local constant regression (NW) to get density estimates that are used in PPIC
den_locNW<-function(partition,t,N){
  n           = length(N)
  l           = length(partition)
  partition_x = {partition[-1]+partition[-l]}/2
  f_locNW     = lapply(seq(1,n),function(i){
    a         = (sum(N[1:i])-N[i]+1):sum(N[1:i])
    ### the histogram density estimate for the ith subject
    f_H       = hist(t[a],breaks=partition,plot=FALSE)$counts/N[i]/diff(partition)
    ### leave one out cross validation, local linear regression
    GCV_locNW = function(hh){
      S   = sapply(partition_x,function(x) x-partition_x)
      S   = dnorm(S/hh)
      S   = apply(S,1,function(x) x/sum(x))
      f_oneout = S%*%f_H
      return(sum({{f_H-f_oneout}/{1-mean(diag(S))}}^2))
    }
    ### find bandwidth that minimize generalized cross validation
    # here set lower limit to be the min(diff(partition_x)) to avoid error
    hh   = optimize(GCV_locNW,lower=min(diff(partition_x)),upper=IQR(t))$minimum
    ### use the obtained GCV bandwidth to get local linear regression estimates
    S   = sapply(t[a],function(x) x-partition_x)
    S   = dnorm(S/hh)
    S   = t(apply(S,2,function(x) x/sum(x)))
    return(S%*%f_H)
  })
}

### local constant
### allows for different bins for different patients
den_locNW2<-function(nbreak=NULL,t,N){
  n           = length(N)
  if(is.null(nbreak)) nbreak = pmax(floor(N/2),50)
  f_locNW   = lapply(seq(1,n),function(i){
    a         = (sum(N[1:i])-N[i]+1):sum(N[1:i])
    ### the histogram density estimate for the ith subject
    f_H       = hist(t[a],breaks=nbreak[i],plot=FALSE)
    partition = f_H$breaks
    f_H       = f_H$counts/N[i]/diff(partition)
    partition_x = {partition[-1]+partition[-length(partition)]}/2
    ### leave one out cross validation, local constant regression (NW)
    GCV_locNW = function(hh){
      S   = sapply(partition_x,function(x) x-partition_x)
      S   = dnorm(S/hh)
      S   = apply(S,1,function(x) x/sum(x))
      f_oneout = S%*%f_H
      return(sum({{f_H-f_oneout}/{1-mean(diag(S))}}^2))
    }
    ### find bandwidth that minimize generalized cross validation
    # here set lower limit to be the min(diff(partition_x)) to avoid error
    hh   = optimize(GCV_NW,lower=min(diff(partition_x)),upper=IQR(t[a]))$minimum
    ### use the obtained GCV bandwidth to get local linear regression estimates
    S   = sapply(t[a],function(x) x-partition_x)
    S   = dnorm(S/hh)
    S   = t(apply(S,2,function(x) x/sum(x)))
    return(S%*%f_H)
  })
}

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



########################################################################
############        Predict Scores for New Subjects         ############     
########################################################################
PP_FPCA_Pred<-function(t,N,mean.fun,eigen.fun,K,parallel=FALSE){
  delta   = mean.fun[2,1]-mean.fun[1,1]
  tmp     = as.numeric(t(eigen.fun[,-1])%*%mean.fun[,-1]*delta) # second term in xi
  
  tmp2    = apply(eigen.fun[,-1],2,function(s) approx(x = mean.fun$x, y = s, xout = t)$y) ### check whether we need parallel this
  indx    = rep(1:length(N),N)
  xi      = -tmp+t(apply(tmp2,2,FUN=function(x) tapply(x,indx,mean))) # FPC scores, ith column for ith patient
  
  if(K==1){
    tmp   = mean.fun[,-1] + outer(eigen.fun[,2],xi[1,]) # density functions
  } else{
    tmp   = as.numeric(mean.fun[,-1]) + as.matrix(eigen.fun[,1:K+1])%*%xi[1:K,] # density functions
  }
  
  tmp     = apply(tmp,2,function(x){
    x[x<0]= 0     # non-negative
    x     = x/sum(x) # integrate to delta^2
    return(x)
  })
  tmp     = tmp/delta
  
  tmp2    = {tmp[-c(1:2),]-tmp[-c(1:2+length(mean.fun[,1])-2),]}/diff(mean.fun[,1],lag=2) 
  
  return(list(scores = t(xi), densities = cbind(mean.fun[,1],tmp), 
              derivatives = cbind({mean.fun[-c(1:2),1]+
                  mean.fun[-c(1:2+length(mean.fun[,1])-2),1]}/2,tmp2)))
}


### make prediction for patients with or without codes
PP_FPCA_Pred2<-function(t,N,mean.fun,eigen.fun,K,parallel=FALSE){
  NZ      = N==0
  
  delta   = mean.fun[2,1]-mean.fun[1,1]
  baseline= as.numeric(t(eigen.fun[,-1])%*%mean.fun[,-1]*delta) # second term in xi
  
  tmp2    = apply(eigen.fun[,-1],2,function(s) approx(x = mean.fun$x, y = s, xout = t)$y) ### check whether we need parallel this
  indx    = rep(1:length(N[!NZ]),N[!NZ])
  xi      = -baseline+t(apply(tmp2,2,FUN=function(x) tapply(x,indx,mean))) # FPC scores, ith column for ith patient
  
  if(K==1){
    tmp   = mean.fun[,-1] + outer(eigen.fun[,2],xi[1,]) # density functions
  } else{
    tmp   = as.numeric(mean.fun[,-1]) + as.matrix(eigen.fun[,1:K+1])%*%xi[1:K,] # density functions
  }
  
  tmp     = apply(tmp,2,function(x){
    x[x<0]= 0     # non-negative
    x     = x/sum(x) # integrate to delta^2
    return(x)
  })
  tmp     = tmp/delta
  
  tmp2    = {tmp[-c(1:2),]-tmp[-c(1:2+length(mean.fun[,1])-2),]}/diff(mean.fun[,1],lag=2) 
  
  return(list(scores = t(xi), densities = cbind(mean.fun[,1],tmp), 
              derivatives = cbind({mean.fun[-c(1:2),1]+
                  mean.fun[-c(1:2+length(mean.fun[,1])-2),1]}/2,tmp2),
              baseline = baseline))
}
