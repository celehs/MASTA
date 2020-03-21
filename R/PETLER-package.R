#' @name PETLER-package
#' @aliases PETLER
#' @docType package
#' @title PETLER Package
#' @description
#' Implements an algorithm to identify the occurrence of event outcome from trajectories of several predictors.
#' @references
#' Wu, S., Müller, H., & Zhang, Z. (2013). FUNCTIONAL DATA ANALYSIS FOR POINT PROCESSES WITH RARE EVENTS. Statistica Sinica, 23(1), 1-23.
#' @useDynLib PETLER
#' @import survival data.table doParallel foreach parallel survC1 rootSolve splines glmnet gglasso rpart rpart.utils  
#' @importFrom Rcpp sourceCpp
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics hist lines plot
#' @importFrom utils read.csv read.table write.table
#' @importFrom stats IQR aggregate approx as.formula binomial bw.SJ bw.bcv bw.nrd bw.nrd0 bw.ucv coef cor density dnorm glm knots median optim optimize prcomp predict quantile runif sd stepfun uniroot var
NULL

#' @name petler.base
#' @title ...
#' @description ...
#' @param data ...
#' @param PPIC_K ...
#' @param n.grid ...
#' @param propvar ...
#' @param n_core ...
#' @export
petler.base <- function(data, PPIC_K = FALSE, n.grid = 401, propvar = 0.85, n_core = 4) {
  TrainSurv <- data.frame(data$TrainSurv)
  ValidSurv <- data.frame(data$ValidSurv)
  TrainCode <- data.frame(data$TrainCode)
  ValidCode <- data.frame(data$ValidCode)
  TrainN <- data.frame(data$TrainN)
  ValidN <- data.frame(data$ValidN)
  TrainSurv_pred_org <- TrainSurv[-(1:4)]
  ValidSurv_pred_org <- TrainSurv[-(1:4)]
  codes <- names(TrainCode)[grep("pred", names(TrainCode))]
  nn <- nrow(TrainSurv) # labeled
  NN <- length(unique(TrainCode$case)) - nn # unlabeled
  nnv <- nrow(ValidSurv) # validation
  Tend <- 1
  TrainPatNum <- unique(TrainCode$case)
  ValidPatNum <- unique(ValidCode$case)
  TrainFU <- unique.data.frame(TrainCode[, 1:2])$analysisfu
  ValidFU <- unique.data.frame(ValidCode[, 1:2])$analysisfu
  TrainCode$monthstd <- TrainCode$month / TrainCode$analysisfu # standardize follow up time
  ValidCode$monthstd <- ValidCode$month / ValidCode$analysisfu # standardize follow up time 
  #--TRAINING---
  K      = NULL
  ft.e   = ft.e.S = PKTS = NULL
  for(i in seq_along(codes)) {
    tmp2    = TrainN[,i+1]>0
    TrainNP = TrainN[tmp2,i+1]
    txt     = paste0("t=with(TrainCode,unlist(sapply(seq_along(month),function(j) rep(month[j],",codes[i],"[j]))))")
    eval(parse(text = txt))
    id      = (1:{nn+NN})[tmp2]
    id      = rep(id,TrainNP)
    PKTS    = cbind(PKTS,GetPK(id = id,t = t,tseq = sort(unique(TrainCode$month)),nn=nrow(TrainN)))
    txt  = paste0("t=with(TrainCode,unlist(sapply(seq_along(monthstd),function(j) rep(monthstd[j],",codes[i],"[j]))))")
    eval(parse(text = txt))
    tmp  = TrainCode[,c("case","month",codes[i])]
    tmp  = tmp[tmp[,3]>0,-3]
    t1   = tapply(tmp[,2],factor(tmp[,1],levels = TrainPatNum),FUN=min)
    t1   = t1[!is.na(t1)]
    registerDoParallel(cores = n_core)
    if(PPIC_K){
      tmp = PP_FPCA_CPP_Parallel(t, h1 = bw.nrd(t), h2 = bw.nrd(t)^{5/6}, TrainNP,
                                 bw = "nrd", ngrid = n.grid, Tend = Tend,
                                 K.select = "PPIC", derivatives = TRUE, nsubs=4)
    } else{
      tmp = PP_FPCA_CPP_Parallel(t, h1 = bw.nrd(t), h2 = bw.nrd(t)^{5/6},
                                 TrainNP, bw = "nrd", ngrid = n.grid,
                                 Tend = Tend, K.select = "PropVar",
                                 propvar = propvar, derivatives = TRUE, nsubs = 4)
    }
    ft.e.tmp = cbind(matrix(TrainFU,nrow=length(TrainFU),ncol=3),-tmp$baseline[1],log(1+TrainN[,i+1]))
    nns  = sum(tmp2)
    ft.e.tmp[tmp2,1]      = t1
    locm = apply(tmp$densities[,1:nns+1],2,which.max)
    ft.e.tmp[tmp2,2] = ft.e.tmp[tmp2,2]*tmp$densities[locm,1]
    ft.e.tmp[tmp2,3] = ft.e.tmp[tmp2,3]*tmp$derivatives[sapply(1:nns,function(i){
      which.max(tmp$derivatives[1:locm[i],i+1])
    }),1]
    ft.e.tmp[tmp2,4] = tmp$scores[,2]
    ft.e = cbind(ft.e,ft.e.tmp)
    ft.e.S.tmp = cbind(ft.e.tmp[,1],VTM(-tmp$baseline[1:4],length(TrainFU)),log(1+TrainN[,i+1]))
    ft.e.S.tmp[tmp2,2:5] = as.matrix(tmp$scores[,2:5])
    ft.e.S = cbind(ft.e.S,ft.e.S.tmp)
    K    = c(K,tmp$K)
    txt = paste0(codes[i],"_FPCA=list(mean=tmp$mean,basis=tmp$basis);rm(tmp)")
    eval(parse(text = txt))
  }
  names(K) = codes
  colnames(ft.e) = paste0(c("1stCode","Pk","ChP","1stScore","logN"), rep(c(seq_along(codes)), each=5))
  colnames(ft.e.S) = paste0(c("1stCode","1stScore","2ndScore","3rdScore","4thScore","logN"), rep(c(seq_along(codes)), each=6))
  colnames(PKTS) = codes
  row.names(ft.e) = row.names(ft.e.S) = row.names(PKTS) = TrainPatNum
  #--Validation---
  txt = paste0("mean.fun=list(",paste0(paste0(codes,"_FPCA$mean"),collapse = ","),")")
  eval(parse(text = txt))
  txt = paste0("basis.fun=list(",paste0(paste0(codes,"_FPCA$basis"),collapse = ","),")")
  eval(parse(text = txt))
  ft.e2 = ft.e.S2 = PKTS2 = NULL
  for(i in seq_along(codes)){
    tmp2    = ValidN[,i+1]>0
    ValidNP = ValidN[tmp2,i+1]
    ### PKs from Two-step procedure
    txt     = paste0("t=with(ValidCode,unlist(sapply(seq_along(month),function(j) rep(month[j],",codes[i],"[j]))))")
    eval(parse(text = txt))
    id      = (1:nnv)[tmp2]
    id      = rep(id,ValidNP)
    PKTS2   = cbind(PKTS2,GetPK(id = id,t = t,tseq = sort(unique(ValidCode$month)),nn=nrow(ValidN)))
    txt  = paste0("t=with(ValidCode,unlist(sapply(seq_along(monthstd),function(j) rep(monthstd[j],",codes[i],"[j]))))")
    eval(parse(text = txt))
    tmp  = ValidCode[,c("case","month",codes[i])]
    tmp  = tmp[tmp[,3]>0,-3]
    t1   = tapply(tmp[,2],factor(tmp[,1],levels = ValidPatNum),FUN=min)
    t1   = t1[!is.na(t1)]
    tmp   = PP_FPCA_Pred2(t,ValidNP,mean.fun[[i]],basis.fun[[i]],K[i])
    ft.e.tmp = cbind(matrix(ValidFU,nrow=length(ValidFU),ncol=3),-tmp$baseline[1],log(1+ValidN[,i+1]))
    nns  = sum(tmp2)
    ft.e.tmp[tmp2,1] = t1
    locm = apply(tmp$densities[,1:nns+1],2,which.max)
    ft.e.tmp[tmp2,2] = ft.e.tmp[tmp2,2]*tmp$densities[locm,1]
    ft.e.tmp[tmp2,3] = ft.e.tmp[tmp2,3]*tmp$derivatives[sapply(1:nns,function(i){
      which.max(tmp$derivatives[1:locm[i],i+1])
    }),1]
    ft.e.tmp[tmp2,4] = tmp$scores[,2]
    ft.e2 = cbind(ft.e2,ft.e.tmp)
    ft.e.S.tmp = cbind(ft.e.tmp[,1],VTM(-tmp$baseline[1:4],length(ValidFU)),log(1+ValidN[,i+1]))
    ft.e.S.tmp[tmp2,2:5] = as.matrix(tmp$scores[,2:5])
    ft.e.S2 = cbind(ft.e.S2,ft.e.S.tmp)
  }
  colnames(ft.e2)   = paste0(c("1stCode","Pk","ChP","1stScore","logN"), rep(c(seq_along(codes)), each=5))
  colnames(ft.e.S2) = paste0(c("1stCode","1stScore","2ndScore","3rdScore","4thScore","logN"), rep(c(seq_along(codes)), each=6))
  colnames(PKTS2)   = codes
  row.names(ft.e2)  = row.names(ft.e.S2) = row.names(PKTS2) = ValidPatNum
  TrainFt   = ft.e[1:nn, ]
  TrainSc   = ft.e.S[1:nn, ]
  ValidFt   = ft.e2
  ValidSc   = ft.e.S2
  TrainPK = PKTS[row.names(PKTS) %in% TrainSurv$case, ]
  ValidPK = PKTS2
  list(nn = nn, 
       codes = codes, 
       Tend = Tend,
       TrainSurv = TrainSurv,
       ValidSurv = ValidSurv,     
       TrainN = TrainN,
       ValidN = ValidN, 
       TrainFt = TrainFt, 
       ValidFt = ValidFt,
       TrainPK = TrainPK,
       ValidPK = ValidPK,
       TrainSurv_pred_org = TrainSurv_pred_org, 
       ValidSurv_pred_org = ValidSurv_pred_org)
}

#' @name petler
#' @title Main function implementing the PETLER algorithm
#' @description ...
#' @param object ...
#' @param cov_group ...
#' @param thresh ...
#' @param PCAthresh ...
#' @param seed ...
#' @param seed2 ...
#' @return A list with components:
#' @return \item{bgbbest_FromChengInit_BFGS}{Details of the fitted model}
#' @return \item{Cstat_BrierSc_ChengInit_BFGS}{Performance of the derived algorithm. C-statistics, etc.}
#' @return \item{group}{A vector of consecutive integers describing the grouping coefficients}
#' @return Several output files will be saved in the \code{outdir} directory.
#' @export
petler <- function(object, cov_group = NULL, thresh = 0.7, PCAthresh = 0.9, seed = 1234, seed2 = 100) {  
  TrainSurv <- object$TrainSurv
  ValidSurv <- object$ValidSurv 
  TrainN <- object$TrainN
  ValidN <- object$ValidN  
  TrainFt <- object$TrainFt
  ValidFt <- object$ValidFt
  TrainPK <- object$TrainPK
  ValidPK <- object$ValidPK
  TrainSurv_pred_org <- object$TrainSurv_pred_org
  ValidSurv_pred_org <- object$ValidSurv_pred_org
  nn <- object$nn 
  codes <- object$codes
  Tend <- obj$Tend  
  #--switch peak & change point if later one is bigger---
  for(i in seq_along(codes)){
    tmp = TrainFt[,5*(i-1)+2] < TrainFt[,5*(i-1)+3]
    if(sum(tmp)>0){
      aa  = TrainFt[tmp,5*(i-1)+2]
      TrainFt[tmp,5*(i-1)+2] = TrainFt[tmp,5*(i-1)+3]
      TrainFt[tmp,5*(i-1)+3] = aa
    }
    tmp = ValidFt[,5*(i-1)+2] < ValidFt[,5*(i-1)+3]
    if(sum(tmp)>0){
      aa  = ValidFt[tmp,5*(i-1)+2]
      ValidFt[tmp,5*(i-1)+2] = ValidFt[tmp,5*(i-1)+3]
      ValidFt[tmp,5*(i-1)+3] = aa
    }
  }
  # set pk & chp to first code time if 0
  for(i in seq_along(codes)){
    tmp = TrainFt[,5*(i-1)+2] == 0
    if(sum(tmp)>0){
      TrainFt[tmp,5*(i-1)+2] = TrainFt[tmp,5*(i-1)+1]
    }
    tmp = TrainFt[,5*(i-1)+3] == 0
    if(sum(tmp)>0){
      TrainFt[tmp,5*(i-1)+3] = TrainFt[tmp,5*(i-1)+1]
    }
    tmp = ValidFt[,5*(i-1)+2] == 0
    if(sum(tmp)>0){
      ValidFt[tmp,5*(i-1)+2] = ValidFt[tmp,5*(i-1)+1]
    }
    tmp = ValidFt[,5*(i-1)+3] == 0
    if(sum(tmp)>0){
      ValidFt[tmp,5*(i-1)+3] = ValidFt[tmp,5*(i-1)+1]
    }
  }
  aa     = c(1:{5*length(codes)})[c(1:{5*length(codes)})%%5%in%c(1,2,3)]
  TrainFt[,aa] = log(TrainFt[,aa])
  ValidFt[,aa] = log(ValidFt[,aa])
  FirstCode = exp(ValidFt[,seq(1,{5*length(codes)},5)])
  txt = paste0("lm(x~", paste(paste0("base_pred", 1:length(TrainSurv_pred_org)), collapse="+"),",data=TrainSurv)")
  func_1 = function(x){
    tmp = eval(parse(text = txt))
    list(resid=tmp$residuals,coef=tmp$coefficients)
  }
  TrainFt = data.frame(TrainFt)
  TrainFt = lapply(TrainFt,func_1)
  RegCoef = do.call(cbind,lapply(TrainFt,`[[`,2))
  TrainFt = do.call(cbind,lapply(TrainFt,`[[`,1))
  txt = paste0("ValidFt[,l]-as.matrix(cbind(1,ValidSurv[,c('", paste(paste0("base_pred",1:length(ValidSurv_pred_org)),collapse="','"), "')]))%*%RegCoef[,l]")
  ValidFt = sapply(1:ncol(ValidFt),function(l){
    eval(parse(text = txt))
  })
  colnames(ValidFt) = colnames(TrainFt)
  SX     = TrainSurv$sx
  SC     = TrainSurv$sc
  Delta  = TrainSurv$delta
  tseq   = sort(unique(SX[Delta==1]))
  G_fit  = survfit(Surv(SX,1-Delta)~1)
  G_tseq = sapply(tseq,function(s){
    tmp2 = G_fit$time<=s
    if(sum(tmp2)==0){
      1
    } else{
      G_fit$surv[max(which(tmp2))]
    }
  })
  G_SX  = summary(G_fit,times=SX,extend=TRUE)$surv
  G_SX[sort(SX,index.return=TRUE)$ix] = G_SX
  G_SX[G_SX==0] = min(G_SX[G_SX>0])
  G_tseq[G_tseq==0] = min(G_tseq[G_tseq>0])
  TrainZC = sapply(seq_along(codes),function(i) sum(TrainN[,i+1]==0))
  names(TrainZC) = codes
  TrainZC/nrow(TrainN)
  ValidZC     = sapply(seq_along(codes),function(i) sum(ValidN[,i+1]==0))
  names(ValidZC) = codes
  ValidZC/nrow(ValidN)
  codes0  = which({TrainZC+ValidZC}/{nrow(TrainN)+nrow(ValidN)} > thresh)
  codes0  = c(sapply(codes0, function(x) (x-1)*5+2:4))
  codesK  = lapply(1:length(codes),function(x){
    tmp = (x-1)*5+1:5
    tmp = tmp[!tmp%in%codes0]
  })
  Z = as.matrix(cbind(TrainSurv[,paste0("base_pred", 1:length(TrainSurv_pred_org))],TrainFt))
  meanZ    = apply(Z,2,mean)
  sdZ      = apply(Z,2,sd)
  for(m in 1:length(TrainSurv_pred_org)){
    aa = unique(TrainSurv[,paste0("base_pred",m)])
    if(length(aa)==2){
      #--categorical -> do not standerdize
      meanZ[paste0("base_pred",m)] = 0
      sdZ[paste0("base_pred",m)]   = 1
    }
  }
  Z = {Z-VTM(meanZ,nrow(Z))}/VTM(sdZ,nrow(Z))
  # PCA within each code to pre-select features
  TrainFt_PCA = lapply(1:length(codes),function(i){
    if(length(codesK[[i]])>0){
      tmp  = Z[,codesK[[i]]+length(TrainSurv_pred_org)] 
      tmp2 = prcomp(tmp)
      PCA_K= sum(cumsum(tmp2$sdev)/sum(tmp2$sdev)<PCAthresh)+1
      tmp3 = tmp2$x[,1:PCA_K]
      list(PCA_K=PCA_K,PCAs = tmp3,PCAobj = tmp2)
    }else{
      list(PCA_K=0,PCAs = NULL,PCAobj=NULL)
    }
  })
  PCA_K    = do.call(c,lapply(TrainFt_PCA,`[[`,1))
  PCAnames = paste0(rep(codes,PCA_K),unlist(sapply(PCA_K,function(x) 1:x)), "_",rep(seq_along(codes),PCA_K))
  PCAobjs  = lapply(TrainFt_PCA,`[[`,3)
  TrainFt_PCA = do.call(cbind,lapply(TrainFt_PCA,`[[`,2))
  colnames(TrainFt_PCA) = PCAnames
  Z        = cbind(Z[,(1:length(TrainSurv_pred_org))],TrainFt_PCA)
  #--group: a vector of consecutive integers describing the grouping of the coefficients
  group = do.call(c,sapply(1:length(PCA_K),function(i) rep(i,PCA_K[i])))
  #--add grouping for covariaates and update the group vevtor--
  if(is.null(cov_group)){
    cov_group = seq(1:length(TrainSurv_pred_org))
  }
  group = c(cov_group, group+cov_group[length(cov_group)])
  corZExtreme(Z,0.7)
  degree = 1
  knots  = quantile(SX[Delta==1],prob=seq(0.1,0.9,0.1))
  Boundary.knots = c(0,max(SX))
  allknots = c(Boundary.knots[1],knots,Boundary.knots[2])
  q      = length(knots)+degree+1
  set.seed(seed)
  func1 = function(bb, Z, Delta, G_SX, SX){estbb.data(bb=bb, Z=Z, Delta=Delta, G_SX=G_SX, SX=SX)$est}
  func2 = function(bb, Z, Delta, G_SX, SX){estbb.data(bb=bb, Z=Z, Delta=Delta, G_SX=G_SX, SX=SX)$Jacob}
  bb    = multiroot(f=func1, start=runif(ncol(Z),-0.1,0.1), jacfunc=func2, Z=Z, Delta=Delta, G_SX=G_SX, SX=SX)
  bb_Cheng = bb.init = bb$root
  ht    = sapply(seq_along(tseq),function(i){
    est = function(a) mean({SX>=tseq[i]}/G_tseq[i]+g_fun(a+Z%*%bb_Cheng))-1
    if(est(-1e10)>0) -1e10 else uniroot(est,lower=-1e10,upper=1e10)$root
  })
  ht    = sort(ht)
  bbt_Cheng = cbind(tseq,ht)
  ht[1] = -30
  weit  = (2/{tseq[-1]+tseq[-length(tseq)]})^{1/4}
  bgi   = function(bg){
    tmpp = h.fun(tseq[-1],knots,Boundary.knots,bg)
    tmp  = {ht[-1]-log(tmpp)}^2
    dtmp = -2*{ht[-1]-log(tmpp)}*dh.fun(tseq[-1],knots,Boundary.knots,bg)/tmpp
    fn   = sum({tmp[-1]+tmp[-length(tmp)]}/2*diff(tseq[-1])*weit[-1])
    gr   = apply({dtmp[-1,]+dtmp[-length(tmp),]}/2*diff(tseq[-1])*weit[-1],2,sum)
    return(list(fn=fn,gr=gr))
  }
  bg.init.Cheng = optim(par = c(-12,5,6,rep(7,7),-3.5)+runif(q,-0.5,0.5),
                        fn = function(bg) bgi(bg)$fn,
                        gr= function(bg) bgi(bg)$gr,method = "BFGS",
                        control = list(maxit = 30000))
  bg.init.Cheng = bg.init.Cheng$par
  bbt.init = log(diff(c(0,exp(ht))))
  tmp = NPMLE.est(bb.init,bbt.init,SX,Z,Delta) ## either use the bb.init from Chen1995 or an arbitrary initial
  bb_NPMLE  = tmp$bb
  bbt_NPMLE = tmp$bbt
  ht    = bbt_NPMLE[,2]
  weit  = (2/{tseq[-1]+tseq[-length(tseq)]})^{1/4}
  bgi   = function(bg){
    tmpp = h.fun(tseq[-1],knots,Boundary.knots,bg)
    tmp  = {ht[-1]-log(tmpp)}^2
    dtmp = -2*{ht[-1]-log(tmpp)}*dh.fun(tseq[-1],knots,Boundary.knots,bg)/tmpp
    fn   = sum({tmp[-1]+tmp[-length(tmp)]}/2*diff(tseq[-1])*weit[-1])
    gr   = apply({dtmp[-1,]+dtmp[-length(tmp),]}/2*diff(tseq[-1])*weit[-1],2,sum)
    return(list(fn=fn,gr=gr))
  }
  aa = bg.init.Cheng
  bg.init.NPMLE = optim(par = aa,fn = function(bg) bgi(bg)$fn,
                        gr= function(bg) bgi(bg)$gr,method = "BFGS",
                        control = list(maxit = 30000))
  bg.init.NPMLE = bg.init.NPMLE$par
  bgbm.init  = c(bg.init.Cheng,bb_Cheng)
  lam        = 0
  bgbm.optim = optim(par = bgbm.init,fn=function(x) SurvLoglik(bg=x[1:q],bb=x[-c(1:q)],knots,Boundary.knots,Delta,SX,Z)$likelihood,
                     gr=function(x) SurvLoglik(bg=x[1:q],bb=x[-c(1:q)],knots,Boundary.knots,Delta,SX,Z,likelihood = FALSE,gradient = TRUE)$gradient,
                     method = "BFGS",control = list(maxit=100000))
  lam.glasso = sort(seq(0,0.003,1e-6),decreasing = TRUE)
  bgbm.optim.glasso = gglasso.Approx.BSNP(bgbm.optim$par[1:q],bgbm.optim$par[-c(1:q)],knots,Boundary.knots,Delta,SX,Z,0,lam.glasso,group)
  mo21      = which.min(bgbm.optim.glasso$AIC.LSA) # AIC
  mo22      = which.min(bgbm.optim.glasso$BIC.LSA) # BIC
  mo23      = which.min(bgbm.optim.glasso$AIC.Orig) # AIC on original model
  mo24      = which.min(bgbm.optim.glasso$BIC.Orig) # BIC on original model
  # non-zero parameters
  nzpar     = cbind(bgbm.optim.glasso$beta[,c(mo21,mo22,mo23,mo24)])!=0
  nzpar     = rbind(matrix(TRUE,nrow=q+4,ncol=4),nzpar[-c(1:4),])
  bgbm.optim.glasso = sapply(1:4,function(i){
    tmp = rep(0,q+ncol(Z))
    tmp[nzpar[,i]] = optim(par = bgbm.optim$par[nzpar[,i]],fn=function(x) SurvLoglik.nz(bg=x[1:q],bb=x[-c(1:q)],nzpar[-c(1:q),i],knots,Boundary.knots,Delta,SX,Z)$likelihood,
                           gr=function(x) SurvLoglik.nz(bg=x[1:q],bb=x[-c(1:q)],nzpar[-c(1:q),i],knots,Boundary.knots,Delta,SX,Z,likelihood = FALSE,gradient = TRUE)$gradient,
                           method = "BFGS",control = list(maxit=10000))$par
    tmp
  })
  bgbm.optim.glasso = list(bgbb = bgbm.optim.glasso,
                           lam.glasso = lam.glasso[c(mo21,mo22,mo23,mo24)])
  colnames(bgbm.optim.glasso$bgbb) = names(bgbm.optim.glasso$lam.glasso) = c("AIC","BIC","AIC.Orig","BIC.Orig")
  bgbbest  = cbind(bgbm.init,bgbm.optim$par,bgbm.optim.glasso$bgbb)
  tree.fit = treefit(Delta,Z[,-c(1:length(TrainSurv_pred_org))]) #??????
  set.seed(seed2)
  train    = sample(1:nn,nn*0.75)
  logi.fit = logifit(Delta,Z,train,colnames(Z)[1:length(TrainSurv_pred_org)])
  vars     = logi.fit$vars
  vars     = sapply(vars,function(x) substr(x,1,nchar(x)-3),USE.NAMES = FALSE)
  vars     = unique(vars)
  wei      = GetWei(TrainPK,vars,Delta,SX)
  SX     = ValidSurv$sx
  SC     = ValidSurv$sc
  Delta  = ValidSurv$delta
  ## same PCA basis as training set
  Z      = as.matrix(cbind(ValidSurv[,paste0("base_pred", 1:length(ValidSurv_pred_org))],ValidFt))
  Z      = {Z-VTM(meanZ,nrow(Z))}/VTM(sdZ,nrow(Z))
  ValidFt_PCA = do.call(cbind,sapply(1:length(codes),function(i){
    if(length(codesK[[i]])>0){
      #tmp  = matrix(Z[,codesK[[i]]+4],ncol=length(codesK[[i]]))
      tmp  = matrix(Z[,codesK[[i]]+length(ValidSurv_pred_org)],ncol=length(codesK[[i]]))
      colnames(tmp) = colnames(Z)[codesK[[i]]+length(ValidSurv_pred_org)]
      tmp3 = predict(PCAobjs[[i]],tmp)[,1:PCA_K[i]]
      list(PCAs = tmp3)
    } else{
      list(PCAs = NULL)
    }
  }))
  colnames(ValidFt_PCA) = PCAnames
  Z        = cbind(Z[,1:length(ValidSurv_pred_org)],ValidFt_PCA)
  endF     = quantile(SX,0.9)
  t.grid   = tseq[seq(9,length(tseq),10)]
  G_fit  = survfit(Surv(SX,1-Delta)~1)
  G_tseq = sapply(tseq,function(s){
    tmp2 = G_fit$time<=s
    if(sum(tmp2)==0){
      1
    } else{
      G_fit$surv[max(which(tmp2))]
    }
  })
  G_SX  = summary(G_fit,times=SX,extend=TRUE)$surv
  G_SX[sort(SX,index.return=TRUE)$ix] = G_SX
  G_SX[G_SX==0] = min(G_SX[G_SX>0])
  G_tseq[G_tseq==0] = min(G_tseq[G_tseq>0])
  tmp   = GetPrecAll.Comb(bgbbest,bb_Cheng,bbt_Cheng,
                          bb_NPMLE,bbt_NPMLE[-1,],SX,SX,SC,Delta,Z,tseq,t.grid,
                          knots,Boundary.knots,G_SX,G_tseq,endF,Tend,
                          tree.fit,FirstCode,
                          logi.fit,vars,wei$wei,wei$adj,as.matrix(ValidPK))
  colnames(tmp) = c("Init","MLE","AIC","BIC",
                    "AIC.Orig","BIC.Orig","Cheng",
                    "NPMLE","TreeTN","TreeAll","logi")
  tmp2 = mean(BrierScore.KM2(tseq[tseq<=endF],SX,SX,SC,Delta,tseq,G_SX,G_tseq)[,2])
  tmp[2,] = 1-tmp[2,]/tmp2
  cstats = tmp
  output=list()
  output$bgbbest_FromChengInit_BFGS   = bgbbest
  output$Cstat_BrierSc_ChengInit_BFGS = cstats
  output$group = group
  output
}
