#' @title Functional Principal Component Analysis (FPCA)
#' @description Performs FPCA to extract features from longitudinal encounter data.
#' @param data input data. See \code{data(data_org)} for example.
#' @param PPIC_K a logical indicating whether you want to use Pseudo-Poisson Information Criterion to choose 
#' the number of principal components K (K.select="PPIC") \code{TRUE} or another criterion to choose 
#' K (K.select="PropVar") \code{FALSE} in the PP_FPCA_Parallel function (hidden). Default is \code{FALSE}.
#' @param n.grid an integer value for grid points used in estimating covariance function g. Default is \code{401}.
#' @param propvar a proportion of variation used to select number of FPCs. Default is \code{0.85}.
#' @param n_core an integer to specify the number of core using for parallel computing. 
#' @export
mata.fpca <- function(data, PPIC_K = FALSE, n.grid = 401, propvar = 0.85, n_core = NULL) {
  if (is.null(n_core)) n_core <- parallel::detectCores()
  registerDoParallel(cores = n_core)  
  TrainSurv <- data.frame(data$TrainSurv)
  ValidSurv <- data.frame(data$ValidSurv)
  TrainCode <- data.frame(data$TrainCode)
  ValidCode <- data.frame(data$ValidCode)
  colnames(TrainSurv) <- colnames(ValidSurv) <- c(
    "case", "delta", "sx", "sc", paste0("base_pred", 1:(NCOL(TrainSurv) - 4)))
  colnames(TrainCode)[1:2] <- colnames(ValidCode)[1:2] <- c("case", "analysisfu")  
  TrainSurv_pred_org <- TrainSurv[-(1:4)]
  ValidSurv_pred_org <- ValidSurv[-(1:4)]  
  codes <- names(TrainCode[, -(1:3)])
  TrainN <- rowsum(TrainCode[, -(1:3)], group = TrainCode[, 1])
  ValidN <- rowsum(ValidCode[, -(1:3)], group = ValidCode[, 1])
  TrainN <- cbind(case = rownames(TrainN), TrainN)
  ValidN <- cbind(case = rownames(ValidN), ValidN)  
  colnames(TrainN)[-1] <- paste0(colnames(TrainCode)[-(1:3)], "_total")
  colnames(ValidN)[-1] <- paste0(colnames(ValidCode)[-(1:3)], "_total")
  TrainPatNum <- unique(TrainCode[, 1])
  ValidPatNum <- unique(ValidCode[, 1])
  nn <- nrow(TrainSurv) # labeled (training)
  nnv <- nrow(ValidSurv) # labeled (validation)
  NN <- length(TrainPatNum) - nn # unlabeled  
  Tend <- 1
  TrainCode$monthstd <- TrainCode$month / TrainCode$analysisfu # standardize follow up time
  ValidCode$monthstd <- ValidCode$month / ValidCode$analysisfu # standardize follow up time 
  TrainFU <- aggregate(TrainCode$analysisfu, list(TrainCode$case), max)  
  ValidFU <- aggregate(ValidCode$analysisfu, list(ValidCode$case), max)
  TrainFU <- TrainFU[match(TrainPatNum, TrainFU[, 1]), 2]
  ValidFU <- ValidFU[match(ValidPatNum, ValidFU[, 1]), 2]
  #--TRAINING---
  K <- NULL
  ft.e <- ft.e.S <- PKTS <- NULL
  FPCA <- vector("list", length(codes))
  names(FPCA) <- codes
  for(i in seq_along(codes)) {
    # cat(i,"\n")
    tmp2 <- TrainN[, i + 1] > 0
    TrainNP <- TrainN[tmp2, i + 1]
    ### PKs from Two-step procedure
    txt <- paste0("t=with(TrainCode,unlist(sapply(seq_along(month),function(j) rep(month[j],",codes[i],"[j]))))")
    eval(parse(text = txt))
    id <- (1:(nn+NN))[tmp2]
    id <- rep(id, TrainNP)
    PKTS <- cbind(PKTS, GetPK(id = id, t = t, tseq = sort(unique(TrainCode$month)), nn = nrow(TrainN)))    
    ### (i)
    txt <- paste0("t=with(TrainCode,unlist(sapply(seq_along(monthstd),function(j) rep(monthstd[j],",codes[i],"[j]))))")
    eval(parse(text = txt))
    tmp <- TrainCode[,c("case","month",codes[i])]
    tmp <- tmp[tmp[, 3] > 0, -3]
    t1 <- tapply(tmp[, 2], factor(tmp[, 1], levels = TrainPatNum), FUN = min)
    t1 <- t1[!is.na(t1)]    
    ### (ii)
    if (PPIC_K) {
      tmp <- PP_FPCA_Parallel(t, h1 = bw.nrd(t), h2 = bw.nrd(t)^{5/6}, TrainNP, bw = "nrd", ngrid = n.grid, Tend = Tend,
                              K.select = "PPIC", derivatives = TRUE, nsubs = 4)
    } else {
      tmp <- PP_FPCA_Parallel(t, h1 = bw.nrd(t), h2 = bw.nrd(t)^{5/6}, TrainNP, bw = "nrd", ngrid = n.grid, Tend = Tend, 
                              K.select = "PropVar", propvar = propvar, derivatives = TRUE, nsubs = 4)
    }    
    ft.e.tmp <- cbind(matrix(TrainFU, nrow = length(TrainFU), ncol = 3), -tmp$baseline[1], log(1 + TrainN[, i + 1]))
    nns <- sum(tmp2)
    ft.e.tmp[tmp2, 1] <- t1
    locm <- apply(tmp$densities[, 1:nns + 1], 2, which.max)
    ft.e.tmp[tmp2, 2] <- ft.e.tmp[tmp2, 2] * tmp$densities[locm, 1]
    ft.e.tmp[tmp2, 3] <- ft.e.tmp[tmp2, 3] * tmp$derivatives[sapply(1:nns, function(i) {
      which.max(tmp$derivatives[1:locm[i], i + 1])
    }), 1]
    ft.e.tmp[tmp2, 4] <- tmp$scores[, 2]
    ft.e <- cbind(ft.e, ft.e.tmp)
    ft.e.S.tmp <- cbind(ft.e.tmp[, 1], VTM(-tmp$baseline[1:4], length(TrainFU)), log(1 + TrainN[, i + 1]))
    ft.e.S.tmp[tmp2, 2:5] <- as.matrix(tmp$scores[, 2:5])
    ft.e.S <- cbind(ft.e.S, ft.e.S.tmp)
    K <- c(K, tmp$K)    
    FPCA[[i]] <- list(K = tmp$K,
                      scores = tmp$scores, 
                      dens = tmp$densities, 
                      deriv = tmp$derivatives,
                      mean = tmp$mean, 
                      basis = tmp$basis, 
                      baseline = tmp$baseline)
    txt <- paste0(codes[i],"_FPCA=list(mean=tmp$mean,basis=tmp$basis);rm(tmp)")
    eval(parse(text = txt))
  }
  names(K) <- codes
  colnames(ft.e) <- paste0(c("1stCode", "Pk", "ChP", "1stScore", "logN"), rep(seq_along(codes), each = 5))
  colnames(ft.e.S) <- paste0(c("1stCode", "1stScore", "2ndScore", "3rdScore", "4thScore", "logN"), rep(seq_along(codes), each = 6))
  colnames(PKTS) <- codes
  row.names(ft.e) <- row.names(ft.e.S) <- row.names(PKTS) <- TrainPatNum
  #--Validation---
  txt <- paste0("mean.fun=list(",paste0(paste0(codes,"_FPCA$mean"),collapse = ","),")")
  eval(parse(text = txt))
  txt <- paste0("basis.fun=list(",paste0(paste0(codes,"_FPCA$basis"),collapse = ","),")")
  eval(parse(text = txt))  
  ft.e2 <- ft.e.S2 <- PKTS2 <- NULL
  for (i in seq_along(codes)) {
    # cat(i, "\n")
    tmp2 <- ValidN[, i + 1] > 0
    ValidNP <- ValidN[tmp2, i + 1]
    names(ValidNP) <- rownames(ValidN[tmp2, ])
    ### PKs from Two-step procedure
    txt <- paste0("t=with(ValidCode,unlist(sapply(seq_along(month),function(j) rep(month[j],",codes[i],"[j]))))")
    eval(parse(text = txt))
    id <- (1:nnv)[tmp2]
    id <- rep(id, ValidNP)
    PKTS2 <- cbind(PKTS2, GetPK(id = id,t = t, tseq = sort(unique(ValidCode$month)), nn = nrow(ValidN)))
    txt <- paste0("t=with(ValidCode,unlist(sapply(seq_along(monthstd),function(j) rep(monthstd[j],",codes[i],"[j]))))")
    eval(parse(text = txt))
    tmp <- ValidCode[, c("case","month", codes[i])]
    tmp <- tmp[tmp[, 3] > 0, -3]
    t1 <- tapply(tmp[, 2], factor(tmp[, 1], levels = ValidPatNum), FUN = min)
    t1 <- t1[!is.na(t1)]
    tmp <- PP.FPCA.Pred(t, ValidNP, mean.fun[[i]], basis.fun[[i]], K[i])
    FPCA[[i]]$ValidPred <- tmp
    ft.e.tmp <- cbind(matrix(ValidFU, nrow = length(ValidFU), ncol = 3), -tmp$baseline[1], log(1 + ValidN[, i + 1]))
    nns <- sum(tmp2)
    ft.e.tmp[tmp2, 1] <- t1
    locm <- apply(tmp$densities[, 1:nns + 1], 2, which.max)
    ft.e.tmp[tmp2, 2] <- ft.e.tmp[tmp2, 2] * tmp$densities[locm, 1]
    ft.e.tmp[tmp2, 3] <- ft.e.tmp[tmp2, 3] * tmp$derivatives[sapply(1:nns, function(i) {
      which.max(tmp$derivatives[1:locm[i], i + 1])
    }), 1]
    ft.e.tmp[tmp2, 4] <- tmp$scores[, 2]
    ft.e2 <- cbind(ft.e2, ft.e.tmp)
    ft.e.S.tmp <- cbind(ft.e.tmp[, 1], VTM(-tmp$baseline[1:4], length(ValidFU)), log(1 + ValidN[, i + 1]))
    ft.e.S.tmp[tmp2, 2:5] <- as.matrix(tmp$scores[, 2:5])
    ft.e.S2 <- cbind(ft.e.S2, ft.e.S.tmp)
  }  
  colnames(ft.e2) <- paste0(c("1stCode", "Pk", "ChP", "1stScore", "logN"), rep(seq_along(codes), each = 5))
  colnames(ft.e.S2) <- paste0(c("1stCode", "1stScore", "2ndScore", "3rdScore", "4thScore", "logN"), rep(seq_along(codes), each = 6))
  colnames(PKTS2) <- codes
  row.names(ft.e2) <- row.names(ft.e.S2) <- row.names(PKTS2) <- ValidPatNum
  TrainFt <- ft.e[1:nn, ]
  TrainSc <- ft.e.S[1:nn, ]
  ValidFt <- ft.e2
  ValidSc <- ft.e.S2
  #---
  TrainPK <- PKTS[row.names(PKTS) %in% TrainSurv$case, ]
  ValidPK <- PKTS2
  list(nn = nn, 
       codes = codes, 
       Tend = Tend,
       TrainSurv = TrainSurv,
       ValidSurv = ValidSurv,     
       TrainSurv_pred_org = TrainSurv_pred_org, 
       ValidSurv_pred_org = ValidSurv_pred_org,       
       TrainN = TrainN,
       ValidN = ValidN, 
       TrainFt = TrainFt, 
       ValidFt = ValidFt,
       TrainPK = TrainPK,
       ValidPK = ValidPK,
       TrainSc = TrainSc,
       ValidSc = ValidSc,
       FPCA = FPCA)
}
