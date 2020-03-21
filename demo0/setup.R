library(foreign)
library(plyr)
library(corrplot)
library(RColorBrewer)
library(rpart)
library(rpart.utils)
library(doParallel)
library(foreach)
library(rootSolve)
library(splines)
library(survC1)
library(gglasso)
library(glmnet)
library(parallel)

source("../R/Cheng1995-edit.R")
source("../R/Likelihood_BS_GrpLasso_RealData.R")
source("../R/NPMLE_RealData.R")
source("../R/PP_FPCA_CPP_Parallel.R")
source("../R/PETLER-package.R")

FPC_Kern_S <- function(x, t, N, h1, h2) {
  grp <- rep(1:length(N), N)
  M <- outer(t, x, "-")
  D1 <- dnorm(M, 0, h1)
  D2 <- dnorm(M, 0, h2)   
  S2 <- rowsum(D2, grp)
  list(f_mu = colSums(D1), 
       Cov_G = crossprod(S2) - crossprod(D2))
}

datadir_org=NULL
datadir_base_func=NULL
outdir=NULL
read_base_func=TRUE
n.grid=401
PPIC_K=FALSE
cov_group=NULL
propvar=0.85
n_core=4
StdFollowUp=TRUE
thresh=0.7
PCAthresh=0.9
seed=1234
seed2=100

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

read_base_func=FALSE
