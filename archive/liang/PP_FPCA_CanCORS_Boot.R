### Baseline variables added

######################################################################
rm(list = ls())
####### source PP_FPCA function
####### Function for loading/installing necessary package
pkgTest <- function(x){
  if (!require(x,character.only = TRUE)){
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}
pkgTest('plyr')
pkgTest('Rcpp')
pkgTest('RcppArmadillo')
pkgTest('foreach')
pkgTest('survival')
pkgTest('rootSolve')
pkgTest('splines')
pkgTest('survC1')
pkgTest('gglasso')
pkgTest('glmnet')
pkgTest('parallel')
pkgTest('rpart')
pkgTest('rpart.utils')
pkgTest('pROC')
pkgTest('abind')

##### get beta using glm, glmnet, etc...
type = "lung"; LabelOnly = FALSE; 
StdFollowUp = TRUE; # ScaleBack = TRUE; be careful about scale back, the score would change
ReadInFt = TRUE; Scores = FALSE; 

############## run the above with std=FALSE, type="crc"/"lung", scores = FALSE/TRUE

### set directory
wkdir  = "~/CanCORS_CRN/"

source(paste0(wkdir,"Likelihood_BS_GrpLasso.R"))
source(paste0(wkdir,"NPMLE.R"))
source(paste0(wkdir,"Cheng1995.R"))

### read in data
ValidSurv = read.csv(paste0(wkdir,"CleanData_new_std/ValidSurv.csv"),
                     stringsAsFactors = FALSE) # validation; survival data
TrainSurv = read.csv(paste0(wkdir,"CleanData_new_std/TrainSurv.csv"),
                     stringsAsFactors = FALSE) # training (labeled); survival data

TrainCode = read.csv(paste0(wkdir,"CleanData_new_std/TrainCode.csv"),
                     stringsAsFactors = FALSE)
ValidCode = read.csv(paste0(wkdir,"CleanData_new_std/ValidCode.csv"),
                     stringsAsFactors = FALSE)

TrainN = read.csv(paste0(wkdir,"CleanData_new_std/TrainN.csv"),
                  stringsAsFactors = FALSE)
ValidN = read.csv(paste0(wkdir,"CleanData_new_std/ValidN.csv"),
                  stringsAsFactors = FALSE)

### Get codes and sample size
codes     = names(TrainCode)[4:14]
nn        = nrow(TrainSurv)  ## labeled
NN        = length(unique(TrainCode$case))-nn ## unlabeled
nnv       = nrow(ValidSurv)  ## validation

TrainPatNum = unique(TrainCode$case)
ValidPatNum = unique(ValidCode$case)

##### standardize follow up time or not
dirpath = paste0("All_",
                 ifelse(StdFollowUp,"StdFollowUp_","UnStdFollowUp_"),
                 ifelse(type=="lung","Lung/","CRC/"))
if(StdFollowUp){
  TrainCode$monthstd = TrainCode[,"month"]/TrainCode[,"analysisfu"] # standardize follow up time
  ValidCode$monthstd = ValidCode[,"month"]/ValidCode[,"analysisfu"] # standardize follow up time
  Tend    = 1
  TrainFU = aggregate(TrainCode$analysisfu,list(TrainCode$case),max)
  TrainFU = TrainFU[match(TrainPatNum,TrainFU[,1]),2]
  ValidFU = aggregate(ValidCode$analysisfu,list(ValidCode$case),max)
  ValidFU = ValidFU[match(ValidPatNum,ValidFU[,1]),2]
} else{
  Tend    = max(TrainCode$month)
}

####### number of patients with no codes in the labeled set (training)+SEER
# corrplot(cor(TrainN[,-1]),title="Total number of code times",
#          mar=c(0,0,1,0),addCoef.col = "black",diag=FALSE,type="upper")  # Mostly postively correlated, high correlation
TrainZC = sapply(seq_along(codes),function(i) sum(TrainN[,i+1]==0))
names(TrainZC) = codes
TrainZC/nrow(TrainN)

####### number of patients with no codes in the labeled set (validation)
# corrplot(cor(ValidN[,-1]),title="Total number of code times",
#          mar=c(0,0,1,0),addCoef.col = "black",diag=FALSE,type="upper")  # Mostly postively correlated, high correlation
ValidZC     = sapply(seq_along(codes),function(i) sum(ValidN[,i+1]==0))
names(ValidZC) = codes
ValidZC/nrow(ValidN)


# delete nsmnoln (high corr nsm); nobs, ~ 95% has no codes
thresh  = 0.9
codesrm = c(which((TrainZC/nrow(TrainN) > thresh | ValidZC/nrow(ValidN)>thresh)),10)
codesrm = codes[codesrm]
codes   = codes[!codes%in%codesrm]
TrainCode = TrainCode[,!colnames(TrainCode)%in%codesrm]
ValidCode = ValidCode[,!colnames(ValidCode)%in%codesrm]
TrainN    = TrainN[,!colnames(TrainN)%in%paste0(codesrm,"_total")]
ValidN    = ValidN[,!colnames(ValidN)%in%paste0(codesrm,"_total")]
TrainZC   = TrainZC[!names(TrainZC)%in%codesrm]
ValidZC   = ValidZC[!names(ValidZC)%in%codesrm]

TrainPK = read.table(paste0(wkdir,"All_StdFollowUp_Lung/StdPKTS_Train.dat"),row.names = 1,header=TRUE)
TrainPK = TrainPK[row.names(TrainPK)%in%TrainSurv$PatNum,]
ValidPK = read.table(paste0(wkdir,"All_StdFollowUp_Lung/StdPKTS_Valid.dat"),row.names = 1,header=TRUE)
ValidPK = ValidPK[row.names(ValidPK)%in%ValidSurv$PatNum,]
FirstCode = read.table(paste0(wkdir,"CleanData_new_std/FirstCode.dat"),
                       header = TRUE)
FirstCode = as.matrix(FirstCode)
# colnames(FirstCode) = codes

###################################################################################
## Read in features (after PCA)
TrainFt_PCA = read.table(paste0(wkdir,"CleanData_new_std/TrainFt_PCA_Orth_Base.dat"),
                         row.names = 1,header=TRUE)
TrainFt_PCA = as.matrix(TrainFt_PCA)
ValidFt_PCA = read.table(paste0(wkdir,"CleanData_new_std/ValidFt_PCA_Orth_Base.dat"),
                         row.names = 1,header=TRUE)
ValidFt_PCA = as.matrix(ValidFt_PCA)
group  = as.numeric(sapply(colnames(TrainFt_PCA)[-c(1:4)],function(x) substr(x,nchar(x),nchar(x))))
group  = c(1,2,3,3,group+3)

#### B-spline setting
degree = 1
knots  = quantile(TrainSurv$SX[TrainSurv$Delta==1],prob=seq(0.1,0.9,0.1))
Boundary.knots = c(0,max(TrainSurv$SX))
allknots = c(Boundary.knots[1],knots,Boundary.knots[2])
q      = length(knots)+degree+1

bb_Cheng  = unlist(read.table(paste0(wkdir,"ParEst_new_std_orth_base/bb_Cheng.dat"),header = TRUE))
bbt_Cheng = read.table(paste0(wkdir,"ParEst_new_std_orth_base/bbt_Cheng.dat"),header = TRUE)

bgbb.init.Cheng = as.matrix(read.table(paste0(wkdir,"ParEst_new_std_orth_base/bgbbest_FromChengInit_BFGS.dat"),
                                       header = TRUE),ncol=6)

lam.glasso = sort(seq(0,0.003,1e-6),decreasing = TRUE)  
endF       = quantile(ValidSurv$SX,0.9)
tseq       = with(TrainSurv,sort(unique(SX[Delta==1])))
t.grid     = 0
FPR.grid   = seq(0.05,0.95,0.025)
simpath    = "Boot_new_std_orth_base/"
###################################################################################
simfunc.boot<-function(X){
  set.seed(X)
  ## Resample the training
  Train.Boot = sample(1:nn,nn,replace = TRUE)
  TrainSurv.Boot = TrainSurv[Train.Boot,]
  TrainPK.Boot   = TrainPK[Train.Boot,]
  Z              = TrainFt_PCA[Train.Boot,]
  
  ################################ survival data
  SX     = TrainSurv.Boot$SX
  SC     = TrainSurv.Boot$SC
  Delta  = TrainSurv.Boot$Delta
  
  ## Cheng 1995
  Cheng.init = Cheng.est(Delta,SX,Z)
  
  ## modify ht a little bit
  yy  = Cheng.init[[2]][,2] # ht
  llx = sum(yy==min(yy))
  if(llx>1){
    for(i in 1:{llx-1}) yy[i] = yy[i]-llx+i
  }
  ## get bg init
  xx   = Cheng.init[[2]][,1]
  weit = (2/{xx[-1]+xx[-length(xx)]})^{1/4}
  bgi  = function(bg){
    tmpp = h.fun(xx[-1],knots,Boundary.knots,bg)
    tmp  = {yy[-1]-log(tmpp)}^2
    dtmp = -2*{yy[-1]-log(tmpp)}*dh.fun(xx[-1],knots,Boundary.knots,bg)/tmpp
    fn   = sum({tmp[-1]+tmp[-length(tmp)]}/2*diff(xx[-1])*weit[-1])
    gr   = apply({dtmp[-1,]+dtmp[-length(tmp),]}/2*diff(xx[-1])*weit[-1],2,sum)
    return(list(fn=fn,gr=gr))
  }
  bg.init = optim(par = c(-12,5,6,rep(7,7),-3.5)+runif(q,-0.5,0.5),
                  fn = function(bg) bgi(bg)$fn,
                  gr= function(bg) bgi(bg)$gr,method = "BFGS",
                  control = list(maxit = 30000))$par
  ## NPMLE
  HZ_NPMLE = log(diff(c(0,exp(yy))))
  bb_NPMLE = NPMLE.est(bb_Cheng,HZ_NPMLE,SX,Z,Delta)
  bbt_NPMLE = bb_NPMLE$bbt
  bb_NPMLE = bb_NPMLE$bb
  
  
  # ours
  # non-zero parameters
  nzpar     = cbind(bgbb.init.Cheng)!=0
  aa        = matrix(c(bg.init,Cheng.init[[1]]),nrow=50,ncol=6)
  bgbbest   = sapply(2:6,function(i){
    tmp = rep(0,q+ncol(Z))
    tmp[nzpar[,i]] = optim(par = aa[nzpar[,i],i],fn=function(x) SurvLoglik.nz(bg=x[1:q],bb=x[-c(1:q)],nzpar[-c(1:q),i],knots,Boundary.knots,Delta,SX,Z)$likelihood,
                           gr=function(x) SurvLoglik.nz(bg=x[1:q],bb=x[-c(1:q)],nzpar[-c(1:q),i],knots,Boundary.knots,Delta,SX,Z,likelihood = FALSE,gradient = TRUE)$gradient,
                           method = "BFGS",control = list(maxit=10000))$par
    tmp
  })
  colnames(bgbbest) = c("MLE","AIC","BIC","AIC.Orig","BIC.Orig")
  write.table(bgbbest,paste0(wkdir,simpath,"bgbbest_comb_",X,".dat"))
  
  # plot(Cheng.init$HZ[,1],Cheng.init$HZ[,2],type="l",ylim=c(-20,10))
  # lines(bbt_NPMLE[,1],bbt_NPMLE[,2],col="red")
  # lines(bbt_NPMLE[,1],log(h.fun(bbt_NPMLE[,1],knots,Boundary.knots,bgbm.init[1:q])),col="blue")
  # lines(bbt_NPMLE[,1],log(h.fun(bbt_NPMLE[,1],knots,Boundary.knots,bgbbest[1:q,6])),col="magenta")
  
  bb.comb  = cbind(bgbbest[-c(1:q),],bb_NPMLE)
  bbt.comb = cbind(bbt_NPMLE[,1],
                   apply(bgbbest[1:q,],2,function(x) log(h.fun(bbt_NPMLE[,1],knots,
                                                               Boundary.knots,x))),
                   bbt_NPMLE[,2])
  
  colnames(bb.comb)=colnames(bbt.comb)[-1]=
    c("MLE","AIC","BIC","AIC.Orig","BIC.Orig","NPMLE")
  
  write.table(bb.comb,paste0(wkdir,simpath,"bb_comb_",X,".dat"))
  write.table(bbt.comb,paste0(wkdir,simpath,"bbt_comb_",X,".dat"))
  
  ## two-step procedures
  row.names(Z) = NULL
  tree.fit = treefit(Delta,Z[,-c(1:4)])
  train    = sample(1:nn,nn*0.75)
  logi.fit = logifit(Delta,Z,train,colnames(Z)[1:4])
  vars     = logi.fit$vars
  vars     = sapply(vars,function(x) substr(x,1,nchar(x)-3),USE.NAMES = FALSE)
  vars     = unique(vars)
  wei      = GetWei(TrainPK.Boot,vars,Delta,SX)
  
  ### validation
  Z      = ValidFt_PCA
  SC     = ValidSurv$SC
  SX     = ValidSurv$SX
  Delta  = ValidSurv$Delta
  
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
  colnames(tmp) = c("MLE","AIC","BIC",
                    "AIC.Orig","BIC.Orig","Cheng",
                    "NPMLE","TreeTN","TreeAll","logi")
  
  tmp2 = mean(BrierScore.KM2(tseq[tseq<=endF],SX,SX,SC,Delta,tseq,G_SX,G_tseq)[,2])
  tmp[2,] = 1-tmp[2,]/tmp2
  
  write.table(tmp,paste0(wkdir,simpath,"CB_",X,".dat"))
  
  tmp2 = GetPrecAll.Extra.Comb(bgbbest,bb_NPMLE,bbt_NPMLE[-1,],bb_Cheng,bbt_Cheng,
                              SX,SC,Delta,Z,tseq,FPR.grid,
                              knots,Boundary.knots,Tend,
                              tree.fit,FirstCode,
                              logi.fit,vars,wei$wei,wei$adj,as.matrix(ValidPK))
  write.table(tmp2$Prec.CV.Extra,
              paste0(wkdir,simpath,"PrecAll_ChengInit_BFGS_",X,".dat"))
  write.table(tmp2$Prec.CV.Extra.TwoStep,
              paste0(wkdir,simpath,"PrecAll_TwoStep_ChengInit_BFGS_",X,".dat"))
  write.table(tmp2$AUC.TwoStep,
              paste0(wkdir,simpath,"PrecAll_AUC_ChengInit_BFGS_",X,".dat"))
  
  list(tmp,tmp2)
}


a       = proc.time()
sim_par = mclapply(X = 1:400, FUN = simfunc.boot, mc.cores = 10)
proc.time()-a

tmp       = lapply(sim_par,`[[`,1)
tmp       = simplify2array(tmp)
write.table(tmp,paste0(combpath,"CB_comb.dat"),row.names = FALSE)

tmp       = lapply(sim_par,`[[`,2)
tmp2      = lapply(tmp,`[[`,1)
tmp2      = simplify2array(tmp2)
write.table(tmp2,paste0(combpath,"PrecAll_comb.dat"),row.names = FALSE)
tmp2      = lapply(tmp,`[[`,2)
tmp2      = simplify2array(tmp2)
write.table(tmp2,paste0(combpath,"PrecAll_TwoStep_comb.dat"),row.names = FALSE)
tmp2      = lapply(tmp,`[[`,3)
tmp2      = simplify2array(tmp2)
write.table(tmp2,paste0(combpath,"PrecAll_AUC_comb.dat"),row.names = FALSE)

