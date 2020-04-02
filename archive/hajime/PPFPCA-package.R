#' @name PPFPCA-package
#' @aliases  PPFPCA-package
#' @docType  package
#' @title PPFPCA
#' @description
#' It builds an algorithm to identify the occurrence of event outcome from trajectories of several predictors.
#' @author Liang Liang, Ming Yang, Tianxi Cai, Hajime Uno
#' Maintainer: Hajime Uno <huno at jimmy.harvard.edu>
#' @references
#' Wu, S., Müller, H., & Zhang, Z. (2013). FUNCTIONAL DATA ANALYSIS FOR POINT PROCESSES WITH RARE EVENTS. Statistica Sinica, 23(1), 1-23.
#' @useDynLib PPFPCA
#' @import RcppArmadillo foreign plyr corrplot RColorBrewer rpart rpart.utils
#'         doParallel foreach survival rootSolve splines survC1 gglasso glmnet parallel
#' @importFrom Rcpp sourceCpp
NULL



#' @name ppfpca
#' @aliases ppfpca
#' @title TBD
#' @description
#' This function builds an algorithm to identify the occurrence of event outcome from trajectories of several predictors.
#' @usage  ppfpca(datadir_org=NULL, datadir_base_func=NULL, outdir=NULL, read_base_func=TRUE,
#'                n.grid=401, PPIC_K=FALSE, propvar=0.85, n_core=4, StdFollowUp=TRUE,
#'                thresh=0.7, PCAthresh=0.9, seed=1234, seed2=100)
#' @param datadir_org a path for the directory where the original data files are saved. 
#' If NULL is specified (default), a directory named "./data_org" will be automatically specified.\cr
#'                    Note that 6 original data must be saved as the following file names;\cr
#'                       1. TrainSurv.csv: baseline survival data for training (labeled); 1st colum: patient id, 2nd colum: event indicator (1=event, 0=censoring),
#'                                         3rd colum: event time, 4th colum: follow-up time, 5th colum--: covariates.\cr
#'                       2. ValidSurv.csv: baseline survival data for validation; baseline survival data for training (labeled); 1st colum: patient id, 
#'                                         2nd colum: event indicator (1=event, 0=censoring), 3rd colum: event time, 4th colum: follow-up time, 
#'                                         5th colum--: covariates.\cr
#'                       3. TrainCode.csv: 1st colum: patient id, 2nd colum: follow-up time, 3rd colum: time label (month), 4th colum-- : predictors.\cr
#'                       4. ValidCode.csv: 1st colum: patient id, 2nd colum: follow-up time, 3rd colum: time label (month), 4th colum-- : predictors.\cr
#'                       5. TrainN.csv:    1st colum: patient id, 2nd colum-- : total number of counts for each predictor.\cr
#'                       6. ValidN.csv:    1st colum: patient id, 2nd colum-- : total number of counts for each predictor.\cr
#'
#' @param datadir_base_func a path for the directory where the base function data will be saved. If NULL is specified (default), a directory named "./data_base_func" will be automatically created under the current working directly.
#' @param outdir a path for the directory where output files will be saved. If NULL is specified (default), a directory named "./outdir" will be automatically created under the current working directly.
#' @param read_base_func a logical indicating whether to create base function data \code{FALSE} or to read base function data files you already created \code{TRUE}. Default is TRUE.
#' @param n.grid an integer value for grid points used in estimating covariance function g. Default is \code{401}.
#' @param PPIC_K a logical indicating whether you want to use Pseudo-Poisson Information Criterion to choose the number of principal components K (K.select="PPIC") \code{TRUE} or another criterion to choose K (K.select="PropVar") \code{FALSE} in the PP_FPCA_CPP_Parallel function (hidden). Default is \code{FALSE}.
#' @param cov_group a vector of consecutive integers describing the grouping only for covariates. When \code{NULL} is specified (default), each covariate will be in different group.
#' @param propvar a proportion of variation used to select number of FPCs. Default is \code{0.85}.
#' @param n_core an integer to specify the number of core using for parallel computing. Default is \code{4}.
#' @param StdFollowUp a logical indicating whether to use standardize follow-up time or not. Standardize follow-up time will be calculated as TrainCode$month/TrainCode$analysisfu).
#' @param thresh a default is \code{0.7}, which means if there are codes with >70\% patients no codes, only use first code time.
#' @param PCAthresh a threshold value for PCA. Default is \code{0.9}.
#' @param seed random seed used for the sampling. Default is \code{1234}.
#' @param seed2 random seed used for the sampling. Default is \code{100}.
#' @return A list with components:
#' @return \item{bgbbest_FromChengInit_BFGS}{Details of the fitted model}
#' @return \item{Cstat_BrierSc_ChengInit_BFGS}{Performance of the derived algorithm. C-statistics, etc.}
#' @return \item{group}{A vector of consecutive integers describing the grouping coefficients}
#' @return Several output files will be saved in the \code{outdir} directory.
#' @details For more details, please contact to the package manager.
#' @references Wu, S., Müller, H., & Zhang, Z. (2013). FUNCTIONAL DATA ANALYSIS FOR POINT PROCESSES WITH RARE EVENTS. Statistica Sinica, 23(1), 1-23.
#' @examples
#' #--------------------------------------------------------------------------
#' # 1. reading base function files already exist (it will take a few min)
#' #--------------------------------------------------------------------------
#' aa=ppfpca(read_base_func=TRUE);
#' print(aa);
#'
#' #-------------------------------------------------------------
#' # 2. creating base function files (it will take more than 2h)
#' #-------------------------------------------------------------
#' bb=ppfpca(read_base_func=FALSE);
#' print(bb);
NULL

#' @export
ppfpca <- function(datadir_org=NULL, datadir_base_func=NULL,
                   outdir=NULL, read_base_func=TRUE, n.grid=401, PPIC_K=FALSE, cov_group=NULL, propvar=0.85, n_core=4,
                   StdFollowUp=TRUE, thresh=0.7, PCAthresh=0.9, seed=1234, seed2=100){
  ##################################
  # Settings
  ##################################
  #--directory for original data--
  if(is.null(datadir_org)){
    datadir_org="./data_org/"
    if(dir.exists(datadir_org)==FALSE){
      dir.create(path=datadir_org)
    }
    if(file.exists(paste0(datadir_org,"TrainSurv.csv"))==FALSE){
      stop("Please save original data files under the directory ./data_org/")
    }
  }

  #--directory for base functions--
  if(is.null(datadir_base_func)){
    datadir_base_func="./data_base_func/"
    if(dir.exists(datadir_base_func)==FALSE){
      dir.create(path=datadir_base_func)
    }
  }

  #--directory for output--
  if(is.null(outdir)){
    outdir="./out/"
    if(dir.exists(outdir)==FALSE){
      dir.create(path=outdir)
    }
  }

  ################################
  ################################
  # func01: Read original data
  ################################
  ################################

  registerDoParallel(cores=n_core)

  TrainSurv = read.csv(paste0(datadir_org,"TrainSurv.csv"), stringsAsFactors = FALSE)
  ValidSurv = read.csv(paste0(datadir_org,"ValidSurv.csv"), stringsAsFactors = FALSE)

  TrainCode = read.csv(paste0(datadir_org,"TrainCode.csv"), stringsAsFactors = FALSE)
  ValidCode = read.csv(paste0(datadir_org,"ValidCode.csv"), stringsAsFactors = FALSE)

  TrainN = read.csv(paste0(datadir_org,"TrainN.csv"), stringsAsFactors = FALSE)
  ValidN = read.csv(paste0(datadir_org,"ValidN.csv"), stringsAsFactors = FALSE)


  #plot(with(TrainSurv,survfit(Surv(SX,Delta)~1),type="kaplan-meier"))
  #lines(with(ValidSurv,survfit(Surv(SX,Delta)~1),type="kaplan-meier"),col="red")


  #===================================================
  # Edit variable names for TrainCode and ValidCode
  #===================================================
  data_chg = c(1, 2)

  for(j in 1:length(data_chg)){
    if(data_chg[j]==1){
      D = TrainCode
      TrainCode_var_org  = names(D)
      TrainCode_pred_org = names(D)[4:ncol(D)]
    }else{
      D = ValidCode
      ValidCode_var_org  = names(D)
      ValidCode_pred_org = names(D)[4:ncol(D)]
    }

    k = ncol(D)
    colnames(D)[1:3] = c("case", "analysisfu", "month")
    colnames(D)[4:k] = paste0("pred", 1:(k-3))

    names(D)

    if(data_chg[j]==1){
      TrainCode = D
    }else{
      ValidCode = D
    }
  }


  #===================================================
  # Edit variable names for TrainN and ValidN
  #===================================================
  data_chg = c(1, 2)

  for(j in 1:length(data_chg)){
    if(data_chg[j]==1){
      D = TrainN
      TrainN_var_org  = names(D)
      TrainN_pred_org = names(D)[2:ncol(D)]
    }else{
      D = ValidN
      ValidN_var_org  = names(D)
      ValidN_pred_org = names(D)[2:ncol(D)]
    }

    k = ncol(D)
    colnames(D)[1] = "case"
    colnames(D)[2:k] = paste0("pred", 1:(k-1), "_total")

    names(D)

    if(data_chg[j]==1){
      TrainN = D
    }else{
      ValidN = D
    }
  }


  #============================================================
  # Edit variable names for TrainSurv and ValidSurv (baseline)
  #============================================================
  data_chg = c(1, 2)

  for(j in 1:length(data_chg)){
    if(data_chg[j]==1){
      D = TrainSurv
      TrainSurv_var_org  = names(D)
      TrainSurv_pred_org = names(D)[5:ncol(D)]
    }else{
      D = ValidSurv
      ValidSurv_var_org  = names(D)
      ValidSurv_pred_org = names(D)[5:ncol(D)]
    }

    k = ncol(D)
    colnames(D)[1:4]  = c("case", "delta", "sx", "sc")
    colnames(D)[5:k]  = paste0("base_pred", 1:(k-4))

    names(D)

    if(data_chg[j]==1){
      TrainSurv = D
    }else{
      ValidSurv = D
    }
  }


  #==========================
  # Resplit the data
  #==========================
  ### Get codes and sample size
  codes     = names(TrainCode)[grep("pred", names(TrainCode))]
  nn        = nrow(TrainSurv)  ## labeled
  NN        = length(unique(TrainCode$case))-nn ## unlabeled
  nnv       = nrow(ValidSurv)  ## validation

  TrainPatNum = unique(TrainCode$case)
  ValidPatNum = unique(ValidCode$case)

  ##### standardize follow up time or not
  dirpath = paste0("All_",
                   ifelse(StdFollowUp,"StdFollowUp/","UnStdFollowUp/"))

  if(dir.exists(paste0(datadir_base_func,dirpath))==FALSE){
    dir.create(path=paste0(datadir_base_func,dirpath))
  }

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


  ##############################################
  ##############################################
  # func02: Create (or read) base functions
  ##############################################
  ##############################################
  if(read_base_func==FALSE){

  #--TRAINING---
  K      = NULL
  ft.e   = ft.e.S = PKTS = NULL
  for(i in seq_along(codes)){
    cat(i,"\n")

    tmp2    = TrainN[,i+1]>0
    TrainNP = TrainN[tmp2,i+1]

    ### PKs from Two-step procedure
    txt     = paste0("t=with(TrainCode,unlist(sapply(seq_along(month),function(j) rep(month[j],",codes[i],"[j]))))")
    eval(parse(text = txt))
    id      = (1:{nn+NN})[tmp2]
    id      = rep(id,TrainNP)
    PKTS    = cbind(PKTS,GetPK(id = id,t = t,tseq = sort(unique(TrainCode$month)),nn=nrow(TrainN)))

    ### (i)
    txt  = paste0("t=with(TrainCode,unlist(sapply(seq_along(monthstd),function(j) rep(monthstd[j],",codes[i],"[j]))))")
    eval(parse(text = txt))
    tmp  = TrainCode[,c("case","month",codes[i])]
    tmp  = tmp[tmp[,3]>0,-3]
    t1   = tapply(tmp[,2],factor(tmp[,1],levels = TrainPatNum),FUN=min)
    t1   = t1[!is.na(t1)]

    ### (ii)
    if(PPIC_K){
      tmp = PP_FPCA_CPP_Parallel(t, h1 = bw.nrd(t), h2 = bw.nrd(t)^{5/6}, TrainNP,
                                 bw = "nrd", ngrid = n.grid, Tend = Tend,
                                 K.select = "PPIC", derivatives = TRUE, nsubs=4)
    } else{
      tmp = PP_FPCA_CPP_Parallel(t, h1 = bw.nrd(t), h2 = bw.nrd(t)^{5/6},
                                 TrainNP, bw = "nrd", ngrid = n.grid,
                                 Tend = Tend, K.select = "PropVar",
                                 propvar = propvar, derivatives = TRUE, nsubs=4)
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

    write.table(tmp$scores,paste0(datadir_base_func, dirpath, codes[i], "_scores.dat"),row.names = FALSE)
    write.table(tmp$densities,paste0(datadir_base_func, dirpath, codes[i], "_dens.dat"),row.names = FALSE)
    write.table(tmp$derivatives,paste0(datadir_base_func, dirpath,  codes[i],"_deriv.dat"),row.names = FALSE)
    write.table(tmp$mean,paste0(datadir_base_func, dirpath, codes[i],"_mean.dat"),row.names = FALSE)
    write.table(tmp$basis,paste0(datadir_base_func, dirpath, codes[i], "_basis.dat"),row.names = FALSE)
    write.table(tmp$baseline,paste0(datadir_base_func, dirpath, codes[i], "_baseline.dat"),row.names = FALSE)

    txt = paste0(codes[i],"_FPCA=list(mean=tmp$mean,basis=tmp$basis);rm(tmp)")
    eval(parse(text = txt))
  }

  names(K) = codes
  write.table(K,paste0(datadir_base_func, dirpath, "FPCAnums.dat"))
  colnames(ft.e) = paste0(c("1stCode","Pk","ChP","1stScore","logN"), rep(c(seq_along(codes)), each=5))
  colnames(ft.e.S) = paste0(c("1stCode","1stScore","2ndScore","3rdScore","4thScore","logN"), rep(c(seq_along(codes)), each=6))
  colnames(PKTS) = codes
  row.names(ft.e) = row.names(ft.e.S) = row.names(PKTS) = TrainPatNum
  write.table(ft.e,paste0(datadir_base_func, dirpath, ifelse(StdFollowUp,"Std","Org"), "Ft_Train.dat"))
  write.table(ft.e.S,paste0(datadir_base_func, dirpath, ifelse(StdFollowUp,"Std","Org"), "Sc_Train.dat"))
  write.table(PKTS,paste0(datadir_base_func, dirpath, ifelse(StdFollowUp,"Std","Org"), "PKTS_Train.dat"))

  cat("done Train\n")


  #--Validation---
  ############################################
  # for(i in codes){
  #   txt = paste0(i,"_FPCA=list(mean=read.table(paste0(datadir,dirpath,'",i,"','_mean.dat'),header=TRUE),
  #                basis=read.table(paste0(datadir,dirpath,'",i,"','_basis.dat'),header=TRUE))")
  #   eval(parse(text = txt))
  # }
  # K   = unlist(read.table(paste0(datadir,dirpath,"FPCAnums.dat")))
  txt = paste0("mean.fun=list(",paste0(paste0(codes,"_FPCA$mean"),collapse = ","),")")
  eval(parse(text = txt))
  txt = paste0("basis.fun=list(",paste0(paste0(codes,"_FPCA$basis"),collapse = ","),")")
  eval(parse(text = txt))

  ft.e2 = ft.e.S2 = PKTS2 = NULL
  for(i in seq_along(codes)){
    cat(i,"\n")
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
  write.table(ft.e2, paste0(datadir_base_func, dirpath, ifelse(StdFollowUp,"Std","Org"), "Ft_Valid.dat"))
  write.table(ft.e.S2, paste0(datadir_base_func, dirpath, ifelse(StdFollowUp,"Std","Org"), "Sc_Valid.dat"))
  write.table(PKTS2, paste0(datadir_base_func, dirpath, ifelse(StdFollowUp,"Std","Org"), "PKTS_Valid.dat"))

  TrainFt   = ft.e[1:nn,]
  TrainSc   = ft.e.S[1:nn,]
  ValidFt   = ft.e2
  ValidSc   = ft.e.S2

  }else{

    #---read existing base function files---
    ft.e    = read.table(paste0(datadir_base_func,dirpath,ifelse(StdFollowUp,"Std","Org"),
                                "Ft_Train.dat"),row.names = 1,header=TRUE)
    ft.e.S  = read.table(paste0(datadir_base_func,dirpath,ifelse(StdFollowUp,"Std","Org"),
                                "Sc_Train.dat"),row.names = 1,header=TRUE)
    TrainFt = ft.e[row.names(ft.e)%in%TrainSurv$case,]
    TrainSc = ft.e.S[row.names(ft.e.S)%in%TrainSurv$case,]
    PKTS = read.table(paste0(datadir_base_func,dirpath,ifelse(StdFollowUp,"Std","Org"),
                                "PKTS_Train.dat"),row.names = 1,header=TRUE)
    ValidFt = read.table(paste0(datadir_base_func,dirpath,ifelse(StdFollowUp,"Std","Org"),
                                "Ft_Valid.dat"),row.names = 1,header=TRUE)
    ValidSc = read.table(paste0(datadir_base_func,dirpath,ifelse(StdFollowUp,"Std","Org"),
                                "Sc_Valid.dat"),row.names = 1,header=TRUE)
    PKTS2 = read.table(paste0(datadir_base_func,dirpath,ifelse(StdFollowUp,"Std","Org"),
                                "PKTS_Valid.dat"),row.names = 1,header=TRUE)

  }

 #---
 TrainPK = PKTS[row.names(PKTS)%in%TrainSurv$case,]
 ValidPK = PKTS2



 ################################
 ################################
 # func03: main part
 ################################
 ################################
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

 ###########################################################
 #####  log scale
 aa     = c(1:{5*length(codes)})[c(1:{5*length(codes)})%%5%in%c(1,2,3)]
 TrainFt[,aa] = log(TrainFt[,aa])
 ValidFt[,aa] = log(ValidFt[,aa])
 FirstCode = exp(ValidFt[,seq(1,{5*length(codes)},5)])
 # write.table(FirstCode,paste0(datadir,"CleanData_new_std/FirstCode.dat"),
 #             row.names = FALSE)

 ###############################################################
 ###   Regress out baseline
 ###############################################################
 txt = paste0("lm(x~", paste(paste0("base_pred", 1:length(TrainSurv_pred_org)), collapse="+"),",data=TrainSurv)")

 func1 = function(x){
   tmp = eval(parse(text = txt))
   list(resid=tmp$residuals,coef=tmp$coefficients)
 }

 TrainFt = data.frame(TrainFt)
 TrainFt = lapply(TrainFt,func1)
 RegCoef = do.call(cbind,lapply(TrainFt,`[[`,2))
 TrainFt = do.call(cbind,lapply(TrainFt,`[[`,1))

 #--
 txt = paste0("ValidFt[,l]-as.matrix(cbind(1,ValidSurv[,c('", paste(paste0("base_pred",1:length(ValidSurv_pred_org)),collapse="','"), "')]))%*%RegCoef[,l]")
 ValidFt = sapply(1:ncol(ValidFt),function(l){
   eval(parse(text = txt))
 })
 colnames(ValidFt) = colnames(TrainFt)


 ###################################################################################
 ### get estimating equation for beta!!!!!!!!!!!!!!
 ################################ survival data
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


 #### get initial values using Cheng et al (1995,1997)



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



 # for codes with > 70% patients no codes, only use first code time
 codes0  = which({TrainZC+ValidZC}/{nrow(TrainN)+nrow(ValidN)} > thresh)
 codes0  = c(sapply(codes0, function(x) (x-1)*5+2:4))
 codesK  = lapply(1:length(codes),function(x){
   tmp = (x-1)*5+1:5
   tmp = tmp[!tmp%in%codes0]
 })

 Z = as.matrix(cbind(TrainSurv[,paste0("base_pred", 1:length(TrainSurv_pred_org))],TrainFt))


 ### standardize Z
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
 # write.table(Z,paste0(datadir,"CleanData_new_std/TrainFt_PCA_Orth_Base.dat"))
 ### try delete ChP, feel like Pk is estimated better


 #### B-spline setting
 degree = 1
 knots  = quantile(SX[Delta==1],prob=seq(0.1,0.9,0.1))
 Boundary.knots = c(0,max(SX))
 allknots = c(Boundary.knots[1],knots,Boundary.knots[2])
 q      = length(knots)+degree+1


 ########################
 #### Initial values
 ########################
 ####--Option I: Cheng et al (1995,1997) Not as good as NPMLE initial
 set.seed(seed)

 func1 = function(bb, Z, Delta, G_SX, SX){estbb.data(bb=bb, Z=Z, Delta=Delta, G_SX=G_SX, SX=SX)$est}
 func2 = function(bb, Z, Delta, G_SX, SX){estbb.data(bb=bb, Z=Z, Delta=Delta, G_SX=G_SX, SX=SX)$Jacob}

 bb    = multiroot(f=func1, start=runif(ncol(Z),-0.1,0.1), jacfunc=func2, Z=Z, Delta=Delta, G_SX=G_SX, SX=SX)
 bb_Cheng = bb.init = bb$root
 # write.table(bb_Cheng,paste0(datadir,"ParEst_new_std_orth_base/bb_Cheng.dat"),row.names = FALSE)
 ht    = sapply(seq_along(tseq),function(i){
   est = function(a) mean({SX>=tseq[i]}/G_tseq[i]+g_fun(a+Z%*%bb_Cheng))-1
   if(est(-1e10)>0) -1e10 else uniroot(est,lower=-1e10,upper=1e10)$root
 })
 ht    = sort(ht)
 bbt_Cheng = cbind(tseq,ht)
 # write.table(bbt_Cheng,paste0(datadir,"ParEst_new_std_orth_base/bbt_Cheng.dat"),row.names = FALSE)


 ### add more weights to origin to get better estimation there
 ### on original scale
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

 pdf(file=paste0(outdir, "bginit_FromCheng-plot1.pdf"))
 plot(tseq[-1],exp(ht[-1]),type="l")
 lines(tseq,h.fun(tseq,knots,Boundary.knots,bg.init.Cheng$par),col="red")
 dev.off()

 pdf(file=paste0(outdir, "bginit_FromCheng-plot2.pdf"))
 plot(tseq[-1],ht[-1],type="l")
 lines(tseq,log(h.fun(tseq,knots,Boundary.knots,bg.init.Cheng$par)),col="red")
 # write.table(bg.init.Cheng$par,paste0(datadir,"ParEst_new_std_orth_base/bginit_FromCheng.dat"),row.names = FALSE)
 dev.off()

 bg.init.Cheng = bg.init.Cheng$par

 ####--Option II: NPMLE (use Cheng as initial)
 bbt.init = log(diff(c(0,exp(ht))))
 # bbt.init = log(diff(c(0,h.fun(tseq,knots,Boundary.knots,bg.init.Cheng))))
 tmp = NPMLE.est(bb.init,bbt.init,SX,Z,Delta) ## either use the bb.init from Chen1995 or an arbitrary initial
 bb_NPMLE  = tmp$bb
 bbt_NPMLE = tmp$bbt
 # write.table(bb_NPMLE,paste0(datadir,"ParEst_new_std_orth_base/bb_NPMLE.dat"),row.names = FALSE)
 # write.table(bbt_NPMLE,paste0(datadir,"ParEst_new_std_orth_base/bbt_NPMLE.dat"),row.names = FALSE)

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


 pdf(file=paste0(outdir, "bginit_FromNPMLE-plot1.pdf"))
 plot(tseq,exp(ht),type="l")
 lines(bbt_NPMLE[,1],exp(bbt_NPMLE[,2]),col="blue")
 lines(tseq,h.fun(tseq,knots,Boundary.knots,bg.init.NPMLE),col="red")
 dev.off()

 pdf(file=paste0(outdir, "bginit_FromNPMLE-plot2.pdf"))
 plot(tseq,ht,type="l")
 lines(bbt_NPMLE[,1],bbt_NPMLE[,2],col="blue")
 lines(tseq,log(h.fun(tseq,knots,Boundary.knots,bg.init.NPMLE)),col="red")
 # write.table(bg.init.NPMLE,paste0(datadir,"ParEst_new_std_orth_base/bginit_FromNPMLE.dat"),row.names = FALSE)
 dev.off()

 # ####--Option III: arbitary initial but with some ridge
 # bg.init = rep(0.1,q)
 # bb.init = rep(0,ncol(Z))
 # lam     = 0

 ## optim with initial bg_bm
 bgbm.init  = c(bg.init.Cheng,bb_Cheng)
 # bgbm.init  = c(bg.init.NPMLE,bb_NPMLE)
 
 lam        = 0
 bgbm.optim = optim(par = bgbm.init,fn=function(x) SurvLoglik(bg=x[1:q],bb=x[-c(1:q)],knots,Boundary.knots,Delta,SX,Z)$likelihood,
                    gr=function(x) SurvLoglik(bg=x[1:q],bb=x[-c(1:q)],knots,Boundary.knots,Delta,SX,Z,likelihood = FALSE,gradient = TRUE)$gradient,
                    method = "BFGS",control = list(maxit=100000))


 ## optim with initial bg_bm
 lam.glasso = sort(seq(0,0.003,1e-6),decreasing = TRUE)
 bgbm.optim.glasso = gglasso.Approx.BSNP(bgbm.optim$par[1:q],bgbm.optim$par[-c(1:q)],knots,Boundary.knots,Delta,SX,Z,0,lam.glasso,group)
 
 
 # AIC
 mo21      = which.min(bgbm.optim.glasso$AIC.LSA)
 # BIC
 mo22      = which.min(bgbm.optim.glasso$BIC.LSA)
 # AIC on original model
 mo23      = which.min(bgbm.optim.glasso$AIC.Orig)
 # BIC on original model
 mo24      = which.min(bgbm.optim.glasso$BIC.Orig)
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

 write.table(bgbbest,paste0(outdir,"bgbbest_FromChengInit_BFGS.dat"),row.names = FALSE)
 #write.table(bgbbest,paste0(datadir,"ParEst_new_std_orth_base/bgbbest_FromNPMLEInit_BFGS.dat"),row.names = FALSE)

 tree.fit = treefit(Delta,Z[,-c(1:length(TrainSurv_pred_org))]) #??????
 set.seed(seed2)
 train    = sample(1:nn,nn*0.75)
 logi.fit = logifit(Delta,Z,train,colnames(Z)[1:length(TrainSurv_pred_org)])
 vars     = logi.fit$vars
 vars     = sapply(vars,function(x) substr(x,1,nchar(x)-3),USE.NAMES = FALSE)
 vars     = unique(vars)

 wei      = GetWei(TrainPK,vars,Delta,SX)



 #####################
 ### validation
 #####################
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
 # write.table(Z,paste0(datadir,"CleanData_new_std/ValidFt_PCA_Orth_Base.dat"))

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
 write.table(cstats,paste0(outdir,"Cstat_BrierSc_ChengInit_BFGS.dat"))
 # write.table(tmp,paste0(outdir,"Cstat_BrierSc_NPMLEInit_BFGS.dat"))



 ################################
 ################################
 # Output
 ################################
 ################################
  output=list()
  output$bgbbest_FromChengInit_BFGS   = bgbbest
  output$Cstat_BrierSc_ChengInit_BFGS = cstats
  output$group = group

  output
}
NULL
