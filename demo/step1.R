################################
################################
# func01: Read original data
################################
################################

registerDoParallel(cores = n_core)

TrainSurv = read.csv(paste0(datadir_org, "TrainSurv.csv"), stringsAsFactors = FALSE)
ValidSurv = read.csv(paste0(datadir_org, "ValidSurv.csv"), stringsAsFactors = FALSE)
TrainCode = read.csv(paste0(datadir_org, "TrainCode.csv"), stringsAsFactors = FALSE)
ValidCode = read.csv(paste0(datadir_org, "ValidCode.csv"), stringsAsFactors = FALSE)
TrainN = read.csv(paste0(datadir_org, "TrainN.csv"), stringsAsFactors = FALSE)
ValidN = read.csv(paste0(datadir_org, "ValidN.csv"), stringsAsFactors = FALSE)

#plot(with(TrainSurv,survfit(Surv(SX,Delta)~1),type="kaplan-meier"))
#lines(with(ValidSurv,survfit(Surv(SX,Delta)~1),type="kaplan-meier"),col="red")

#===================================================
# Edit variable names for TrainCode and ValidCode
#===================================================
data_chg = c(1, 2)
for(j in 1:length(data_chg)) {
  if (data_chg[j] == 1) {
    D = TrainCode
    TrainCode_var_org  = names(D)
    TrainCode_pred_org = names(D)[4:ncol(D)]
  } else {
    D = ValidCode
    ValidCode_var_org  = names(D)
    ValidCode_pred_org = names(D)[4:ncol(D)]
  }
  k = ncol(D)
  colnames(D)[1:3] = c("case", "analysisfu", "month")
  colnames(D)[4:k] = paste0("pred", 1:(k-3))
  # names(D)
  if(data_chg[j] == 1) {
    TrainCode = D
  } else {
    ValidCode = D
  }
}

#===================================================
# Edit variable names for TrainN and ValidN
#===================================================
data_chg = c(1, 2)
for(j in 1:length(data_chg)) {
  if (data_chg[j] == 1) {
    D = TrainN
    TrainN_var_org  = names(D)
    TrainN_pred_org = names(D)[2:ncol(D)]
  } else {
    D = ValidN
    ValidN_var_org  = names(D)
    ValidN_pred_org = names(D)[2:ncol(D)]
  }
  k = ncol(D)
  colnames(D)[1] = "case"
  colnames(D)[2:k] = paste0("pred", 1:(k - 1), "_total")
  # names(D)
  if (data_chg[j] == 1) {
    TrainN = D
  } else {
    ValidN = D
  }
}

#============================================================
# Edit variable names for TrainSurv and ValidSurv (baseline)
#============================================================
data_chg = c(1, 2)
for (j in 1:length(data_chg)) {
  if (data_chg[j] == 1) {
    D = TrainSurv
    TrainSurv_var_org  = names(D)
    TrainSurv_pred_org = names(D)[5:ncol(D)]
  } else {
    D = ValidSurv
    ValidSurv_var_org  = names(D)
    ValidSurv_pred_org = names(D)[5:ncol(D)]
  }
  k = ncol(D)
  colnames(D)[1:4] = c("case", "delta", "sx", "sc")
  colnames(D)[5:k] = paste0("base_pred", 1:(k-4))
  # names(D)
  if (data_chg[j] == 1) {
    TrainSurv = D
  } else {
    ValidSurv = D
  }
}

#==========================
# Resplit the data
#==========================
### Get codes and sample size
codes = names(TrainCode)[grep("pred", names(TrainCode))]
nn = nrow(TrainSurv)  ## labeled
NN = length(unique(TrainCode$case)) - nn ## unlabeled
nnv = nrow(ValidSurv)  ## validation
TrainPatNum = unique(TrainCode$case)
ValidPatNum = unique(ValidCode$case)

##### standardize follow up time or not
dirpath = paste0("All_", ifelse(StdFollowUp, "StdFollowUp/", "UnStdFollowUp/"))

if (dir.exists(paste0(datadir_base_func, dirpath)) == FALSE) {
  dir.create(path = paste0(datadir_base_func, dirpath))
}

if (StdFollowUp) {
  TrainCode$monthstd = TrainCode[,"month"] / TrainCode[,"analysisfu"] # standardize follow up time
  ValidCode$monthstd = ValidCode[,"month"] / ValidCode[,"analysisfu"] # standardize follow up time
  Tend = 1
  TrainFU = aggregate(TrainCode$analysisfu, list(TrainCode$case), max)
  TrainFU = TrainFU[match(TrainPatNum, TrainFU[, 1]), 2]
  ValidFU = aggregate(ValidCode$analysisfu, list(ValidCode$case), max)
  ValidFU = ValidFU[match(ValidPatNum, ValidFU[, 1]), 2]
} else {
  Tend = max(TrainCode$month)
}

##############################################
##############################################
# func02: Create (or read) base functions
##############################################
##############################################
if(read_base_func == FALSE) {
  
  #--TRAINING---
  K = NULL
  ft.e = ft.e.S = PKTS = NULL
  for(i in seq_along(codes)) {
    cat(i,"\n")
    tmp2 = TrainN[, i + 1] > 0
    TrainNP = TrainN[tmp2, i + 1]
    
    ### PKs from Two-step procedure
    txt = paste0("t=with(TrainCode,unlist(sapply(seq_along(month),function(j) rep(month[j],", codes[i], "[j]))))")
    eval(parse(text = txt))
    id = (1:{nn + NN})[tmp2]
    id = rep(id, TrainNP)
    PKTS = cbind(PKTS, GetPK(id = id, t = t, tseq = sort(unique(TrainCode$month)), nn = nrow(TrainN)))
    
    ### (i)
    txt  = paste0("t=with(TrainCode,unlist(sapply(seq_along(monthstd),function(j) rep(monthstd[j],", codes[i], "[j]))))")
    eval(parse(text = txt))
    tmp = TrainCode[, c("case", "month", codes[i])]
    tmp = tmp[tmp[, 3] > 0, -3]
    t1 = tapply(tmp[, 2], factor(tmp[, 1], levels = TrainPatNum), FUN = min)
    t1 = t1[!is.na(t1)]
    
    ### (ii)
    if(PPIC_K){
      tmp = PP_FPCA_CPP_Parallel(t, h1 = bw.nrd(t), h2 = bw.nrd(t)^{5/6}, TrainNP, bw = "nrd", ngrid = n.grid, Tend = Tend,
                                 K.select = "PPIC", derivatives = TRUE, nsubs = 4)
    } else{
      tmp = PP_FPCA_CPP_Parallel(t, h1 = bw.nrd(t), h2 = bw.nrd(t)^{5/6}, TrainNP, bw = "nrd", ngrid = n.grid, Tend = Tend, 
                                 K.select = "PropVar", propvar = propvar, derivatives = TRUE, nsubs = 4)
    }
    
    ft.e.tmp = cbind(matrix(TrainFU, nrow = length(TrainFU), ncol = 3), -tmp$baseline[1], log(1 + TrainN[, i + 1]))
    nns = sum(tmp2)
    ft.e.tmp[tmp2,1] = t1
    locm = apply(tmp$densities[, 1:nns + 1], 2, which.max)
    ft.e.tmp[tmp2, 2] = ft.e.tmp[tmp2, 2] * tmp$densities[locm, 1]
    ft.e.tmp[tmp2, 3] = ft.e.tmp[tmp2, 3] * tmp$derivatives[sapply(1:nns, function(i) {
      which.max(tmp$derivatives[1:locm[i], i + 1]) 
    }), 1]
    ft.e.tmp[tmp2, 4] = tmp$scores[, 2]
    ft.e = cbind(ft.e, ft.e.tmp)
    ft.e.S.tmp = cbind(ft.e.tmp[, 1], VTM(-tmp$baseline[1:4], length(TrainFU)), log(1 + TrainN[, i + 1]))
    ft.e.S.tmp[tmp2, 2:5] = as.matrix(tmp$scores[, 2:5])
    ft.e.S = cbind(ft.e.S, ft.e.S.tmp)
    K = c(K,tmp$K)
    write.table(tmp$scores,paste0(datadir_base_func, dirpath, codes[i], "_scores.dat"), row.names = FALSE)
    write.table(tmp$densities,paste0(datadir_base_func, dirpath, codes[i], "_dens.dat"), row.names = FALSE)
    write.table(tmp$derivatives,paste0(datadir_base_func, dirpath,  codes[i],"_deriv.dat"), row.names = FALSE)
    write.table(tmp$mean,paste0(datadir_base_func, dirpath, codes[i],"_mean.dat"), row.names = FALSE)
    write.table(tmp$basis,paste0(datadir_base_func, dirpath, codes[i], "_basis.dat"), row.names = FALSE)
    write.table(tmp$baseline,paste0(datadir_base_func, dirpath, codes[i], "_baseline.dat"), row.names = FALSE)
    txt = paste0(codes[i], "_FPCA=list(mean=tmp$mean,basis=tmp$basis);rm(tmp)")
    eval(parse(text = txt))
  }
  
  names(K) = codes
  write.table(K,paste0(datadir_base_func, dirpath, "FPCAnums.dat"))
  colnames(ft.e) = paste0(c("1stCode","Pk","ChP","1stScore","logN"), rep(c(seq_along(codes)), each = 5))
  colnames(ft.e.S) = paste0(c("1stCode","1stScore","2ndScore","3rdScore","4thScore","logN"), rep(c(seq_along(codes)), each = 6))
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
  
} else {
  
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
