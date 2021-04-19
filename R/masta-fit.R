#' @title Main function implementing the MASTA algorithm
#' @description This function builds an algorithm to identify the occurrence of event outcome from trajectories of several predictors.
#' @param object The object returned by the \code{fpca.combine} function
#' @param survival the labeled data. The columns should be 1) id, 2) event indicator, 3) event time, followed by baseline predictors.
#' @param follow_up_time the follow-up data
#' @param cov_group a vector of consecutive integers describing the grouping only for covariates. When \code{NULL} is specified (default), each covariate will be in different group.
#' @param thresh a default is \code{0.7}, which means if there are codes with >70\% patients no codes, only use first code time.
#' @param PCAthresh a threshold value for PCA. Default is \code{0.9}.
#' @param seed random seed used for the sampling. Default is \code{100}.
#' @return A list with components:
#' @return \item{bgbbest_FromChengInit_BFGS}{Details of the fitted model}
#' @return \item{Cstat_BrierSc_ChengInit_BFGS}{Performance of the derived algorithm. C-statistics, etc.}
#' @return \item{group}{A vector of consecutive integers describing the grouping coefficients}
#' @export
masta.fit <- function(object, survival, follow_up_time, Tend=1, cov_group = NULL, thresh = 0.7, PCAthresh = 0.9, seed = 100) {

  aa = merge(survival, follow_up_time, by="id")
  
  TrainSurv=aa[aa$train_valid==1,-ncol(aa)]
  ValidSurv=aa[aa$train_valid==2,-ncol(aa)]
  npred = ncol(aa)-5
  codes <- paste0("pred", 1:npred)
  colnames(TrainSurv) <- colnames(ValidSurv) <- c("case", "delta", "sx",  paste0("base_pred", 1:npred), "sc")
  nn <- nrow(TrainSurv)
  TrainSurv_pred_org <- TrainSurv[, grep("base_pred",names(TrainSurv))] #--- baseline predictors
  ValidSurv_pred_org <- ValidSurv[, grep("base_pred",names(ValidSurv))] #--- baseline predictors

  TrainN <- object$TrainN
  ValidN <- object$ValidN
  TrainFt <- object$TrainFt[row.names(object$TrainFt) %in% as.character(TrainSurv$case), ]
  ValidFt <- object$ValidFt
  TrainPK <- object$TrainPK[row.names(object$TrainPK) %in% as.character(TrainSurv$case), ]
  ValidPK <- object$ValidPK
  
  
  for (i in seq_along(codes)) {
    idx0 <- 5 * (i - 1)
    #--switch peak & change point if later one is bigger---    
    tmp <- TrainFt[, idx0 + 2] < TrainFt[, idx0 + 3]
    if (sum(tmp) > 0) {
      aa <- TrainFt[tmp, idx0 + 2]
      TrainFt[tmp, idx0 + 2] <- TrainFt[tmp, idx0 + 3]
      TrainFt[tmp, idx0 + 3] <- aa
    }
    tmp <- ValidFt[, idx0 + 2] < ValidFt[, idx0 + 3]
    if (sum(tmp) > 0) {
      aa <- ValidFt[tmp, idx0 + 2]
      ValidFt[tmp, idx0 + 2] <- ValidFt[tmp, idx0 + 3]
      ValidFt[tmp, idx0 + 3] <- aa
    }
    # set pk & chp to first code time if 0
    tmp <- TrainFt[, idx0 + 2] == 0
    if (sum(tmp) > 0) TrainFt[tmp, idx0 + 2] <- TrainFt[tmp, idx0 + 1]
    tmp <- TrainFt[, idx0 + 3] == 0
    if (sum(tmp) > 0) TrainFt[tmp, idx0 + 3] <- TrainFt[tmp, idx0 + 1]
    tmp <- ValidFt[, idx0 + 2] == 0
    if (sum(tmp) > 0) ValidFt[tmp, idx0 + 2] <- ValidFt[tmp, idx0 + 1]
    tmp <- ValidFt[, idx0 + 3] == 0
    if (sum(tmp) > 0) ValidFt[tmp, idx0 + 3] <- ValidFt[tmp, idx0 + 1]
  }
  ##### log scale
  ncodes <- length(codes)
  tmp <- seq(5 * ncodes)
  aa <- tmp[tmp %% 5 %in% 1:3]
  TrainFt[, aa] <- log(TrainFt[, aa])
  ValidFt[, aa] <- log(ValidFt[, aa])
  FirstCode <- exp(ValidFt[, seq(1, 5 * ncodes, 5)])

  ### Regress out baseline ##
  tmp <- paste(paste0("base_pred", 1:length(TrainSurv_pred_org)), collapse = " + ")
  txt <- paste0("lm(x ~ ", tmp, ", data = TrainSurv)")
  func_1 <- function(x) {
    tmp <- eval(parse(text = txt))
    list(resid = tmp$residuals, coef = tmp$coefficients)
  }
  TrainFt <- data.frame(TrainFt)
  TrainFt <- lapply(TrainFt, func_1)
  RegCoef <- do.call(cbind, lapply(TrainFt, `[[`, 2))
  TrainFt <- do.call(cbind, lapply(TrainFt, `[[`, 1))

  tmp <- paste(paste0("base_pred", 1:length(ValidSurv_pred_org)), collapse = "', '")
  txt <- paste0("ValidFt[, l] - as.matrix(cbind(1, ValidSurv[, c('", tmp, "')])) %*% RegCoef[, l]")
  ValidFt <- sapply(1:ncol(ValidFt), function(l) eval(parse(text = txt)))
  colnames(ValidFt) <- colnames(TrainFt)


  ### get estimating equation for beta!!!!!!!!!!!!!!
  SX <- TrainSurv$sx
  SC <- TrainSurv$sc
  Delta <- TrainSurv$delta
  tseq <- sort(unique(SX[Delta == 1]))
  G_fit <- survfit(Surv(SX, 1 - Delta) ~ 1)
  G_tseq <- sapply(tseq, function(s) {
    tmp2 <- G_fit$time <= s
    if (sum(tmp2) == 0) {
      1
    } else {
      G_fit$surv[max(which(tmp2))]
    }
  })
  G_SX <- summary(G_fit, times = SX, extend = TRUE)$surv
  G_SX[sort(SX, index.return = TRUE)$ix] <- G_SX
  G_SX[G_SX == 0] <- min(G_SX[G_SX > 0])
  G_tseq[G_tseq == 0] <- min(G_tseq[G_tseq > 0])

  #--- for output ---
  G_SX_train <- G_SX
  G_tseq_train <- G_tseq

  #### get initial values using Cheng et al (1995,1997)
  ####### number of patients with no codes in the labeled set (training)+SEER
  TrainZC <- sapply(seq_along(codes), function(i) sum(TrainN[, i + 1] == 0))
  names(TrainZC) <- codes
  ####### number of patients with no codes in the labeled set (validation)
  ValidZC <- sapply(seq_along(codes), function(i) sum(ValidN[, i + 1] == 0))
  names(ValidZC) <- codes
  # for codes with > 70% patients no codes, only use first code time
  codes0 <- which((TrainZC + ValidZC) / (nrow(TrainN) + nrow(ValidN)) > thresh)
  codes0 <- c(sapply(codes0, function(x) (x - 1) * 5 + 2:4))
  codesK <- lapply(1:length(codes), function(x) {
    tmp <- (x - 1) * 5 + 1:5
    tmp <- tmp[!tmp %in% codes0]
  })
  # Z <- as.matrix(cbind(TrainSurv[, paste0("base_pred", 1:length(TrainSurv_pred_org))], TrainFt))
  Z <- as.matrix(cbind(TrainSurv[, grep("base_pred", names(TrainSurv))], TrainFt))
  ### standardize Z
  meanZ <- apply(Z, 2, mean)
  sdZ <- apply(Z, 2, sd)
  for (m in 1:length(TrainSurv_pred_org)) {
    aa <- unique(TrainSurv[, paste0("base_pred", m)])
    if (length(aa) == 2) {
      #--categorical -> do not standerdize
      meanZ[paste0("base_pred", m)] <- 0
      sdZ[paste0("base_pred", m)] <- 1
    }
  }
  Z <- (Z - VTM(meanZ, nrow(Z))) / VTM(sdZ, nrow(Z))

  # PCA within each code to pre-select features
  TrainFt_PCA <- lapply(1:length(codes), function(i) {
    if (length(codesK[[i]]) > 0) {
      tmp <- Z[, codesK[[i]] + length(TrainSurv_pred_org)]
      tmp2 <- prcomp(tmp)
      PCA_K <- sum(cumsum(tmp2$sdev) / sum(tmp2$sdev) < PCAthresh) + 1
      tmp3 <- tmp2$x[, 1:PCA_K]
      list(PCA_K = PCA_K, PCAs = tmp3, PCAobj = tmp2)
    } else {
      list(PCA_K = 0, PCAs = NULL, PCAobj = NULL)
    }
  })
  PCA_K <- do.call(c, lapply(TrainFt_PCA, `[[`, 1))
  PCAnames <- paste0(rep(codes, PCA_K), unlist(sapply(PCA_K, function(x) 1:x)), "_", rep(seq_along(codes), PCA_K))
  PCAobjs <- lapply(TrainFt_PCA, `[[`, 3)
  TrainFt_PCA <- do.call(cbind, lapply(TrainFt_PCA, `[[`, 2))
  colnames(TrainFt_PCA) <- PCAnames
  Z <- cbind(Z[, 1:length(TrainSurv_pred_org)], TrainFt_PCA)
  #--group: a vector of consecutive integers describing the grouping of the coefficients
  # group <- do.call(c, sapply(1:length(PCA_K), function(i) rep(i, PCA_K[i]))) #<--- this does not work if the elements of PCA_K are identical
  group <- rep(1:length(PCA_K), PCA_K)

  #--add grouping for covariaates and update the group vevtor--
  if (is.null(cov_group)) cov_group <- 1:length(TrainSurv_pred_org)
  group <- c(cov_group, group + cov_group[length(cov_group)])
  
  #corZExtreme(Z, 0.7)
  ### try delete ChP, feel like Pk is estimated better




  #### B-spline setting (9 knots plus 2 = 11 --- based on quantiles)
  degree <- 1
  knots <- quantile(SX[Delta == 1], prob = seq(0.1, 0.9, 0.1))
  Boundary.knots <- c(0, max(SX))
  allknots <- c(Boundary.knots[1], knots, Boundary.knots[2])
  q <- length(knots) + degree + 1

  #-- for output 
  nparm4Bspline <- q

#### ===================
#### Initial values
#### ===================
#  set.seed(seed)
#  if (0) {
#    b_ini <- runif(ncol(Z), -0.1, 0.1)
#    b_next <- b_ini
#    b_next <- estbb.data(bb = b_next, Z = Z, Delta = Delta, G_SX = G_SX, SX = SX)$est
#    #### --Option I: Cheng et al (1995,1997) Not as good as NPMLE initial
#    func1 <- function(bb, Z, Delta, G_SX, SX) estbb.data(bb = bb, Z = Z, Delta = Delta, G_SX = G_SX, SX = SX)$est
#    func2 <- function(bb, Z, Delta, G_SX, SX) estbb.data(bb = bb, Z = Z, Delta = Delta, G_SX = G_SX, SX = SX)$Jacob
#    bb <- multiroot(f = func1, start = b_ini, jacfunc = func2, Z = Z, Delta = Delta, G_SX = G_SX, SX = SX)
#    bb_Cheng <- bb.init <- bb$root
#  }

  bb_Cheng <- rep(0, ncol(Z))
  ht <- sapply(seq_along(tseq), function(i) {
    est <- function(a) mean((SX >= tseq[i]) / G_tseq[i] + g_fun(a + Z %*% bb_Cheng)) - 1
    if (est(-1e10) > 0) -1e10 else uniroot(est, lower = -1e10, upper = 1e10)$root
  })
  ht <- sort(ht)
  bbt_Cheng <- cbind(tseq, ht)

  ### add more weights to origin to get better estimation there ???
  ### on original scale
  # ht[1] <- -30  --- will give an error if ht[2] is less than 30
  ht[1] <- min(ht[1], ht[2] - 10)
  weit <- (2 / (tseq[-1] + tseq[-length(tseq)]))^(1 / 4)
  bgi <- function(bg) {
    tmpp <- h.fun(tseq[-1], knots, Boundary.knots, bg)
    tmp <- (ht[-1] - log(tmpp))^2
    dtmp <- -2 * (ht[-1] - log(tmpp)) * dh.fun(tseq[-1], knots, Boundary.knots, bg) / tmpp
    fn <- sum((tmp[-1] + tmp[-length(tmp)]) / 2 * diff(tseq[-1]) * weit[-1])
    gr <- apply((dtmp[-1, ] + dtmp[-length(tmp), ]) / 2 * diff(tseq[-1]) * weit[-1], 2, sum)
    return(list(fn = fn, gr = gr))
  }
  ini <- c(-12, 5, 6, rep(7, 7), -3.5) + runif(q, -0.5, 0.5)
  bg.init.Cheng <- optim(
    par = ini,
    fn = function(bg) bgi(bg)$fn,
    gr = function(bg) bgi(bg)$gr,
    method = "BFGS", control = list(maxit = 30000)
  )
  bg.init.Cheng <- bg.init.Cheng$par


  ###### Option III: arbitary initial but with some ridge
  ## optim with initial bg_bm
  bgbm.init <- c(bg.init.Cheng, bb_Cheng)
  lam <- 0
  bgbm.optim <- optim(
    par = bgbm.init,
    fn = function(x) {
      SurvLoglik(
        bg = x[1:q], bb = x[-(1:q)], knots, Boundary.knots, Delta, SX, Z
      )$likelihood
    },
    gr = function(x) {
      SurvLoglik(
        bg = x[1:q], bb = x[-(1:q)], knots, Boundary.knots, Delta, SX, Z, likelihood = FALSE, gradient = TRUE
      )$gradient
    },
    method = "BFGS",
    control = list(maxit = 100000)
  )


  #  inital.values = cbind(bgbm.init, c(bb_NPMLE, bg.init.NPMLE),bgbm.optim$par)
  #  colnames(inital.values)=c("Cheng (Option 1)","NPMLE (Option 2)","Ridge (Option 3)")
  #  print(inital.values)

  bgbm.optim
  #--- gamma ---
  cbind(bg.init.Cheng, bgbm.optim$par[1:q])
  #--- beta ---
  cbind(bb_Cheng, bgbm.optim$par[-c(1:q)])

  # ===================
  # Group Lasso
  # ===================
  ## optim with initial bgbm
  lam.glasso <- sort(seq(0, 0.003, 1e-6), decreasing = TRUE)

  bgbm.optim.glasso <- gglasso.Approx.BSNP(
    bgbm.optim$par[1:q], bgbm.optim$par[-c(1:q)], knots, Boundary.knots, Delta, SX, Z, 0, lam.glasso, group
  )

  mo21 <- which.min(bgbm.optim.glasso$AIC.LSA) # AIC
  mo22 <- which.min(bgbm.optim.glasso$BIC.LSA) # BIC
  mo23 <- which.min(bgbm.optim.glasso$AIC.Orig) # AIC on original model
  mo24 <- which.min(bgbm.optim.glasso$BIC.Orig) # BIC on original model
  nzpar <- cbind(bgbm.optim.glasso$beta[, c(mo21, mo22, mo23, mo24)]) != 0 # non-zero parameters
  nzpar <- rbind(matrix(TRUE, nrow = q + 4, ncol = 4), nzpar[-(1:4), ])


  bgbm.optim.glasso <- sapply(1:4, function(i) {
    tmp <- rep(0, q + ncol(Z))
    tmp[nzpar[, i]] <- optim(
      par = bgbm.optim$par[nzpar[, i]],
      fn = function(x) {
        SurvLoglik.nz(
          bg = x[1:q], bb = x[-(1:q)], nzpar[-(1:q), i], knots, Boundary.knots, Delta, SX, Z
        )$likelihood
      },
      gr = function(x) {
        SurvLoglik.nz(
          bg = x[1:q], bb = x[-c(1:q)], nzpar[-c(1:q), i], knots, Boundary.knots, Delta, SX, Z,
          likelihood = FALSE, gradient = TRUE
        )$gradient
      },
      method = "BFGS",
      control = list(maxit = 10000)
    )$par
    tmp
  })

  bgbm.optim.glasso <- list(
    bgbb = bgbm.optim.glasso,
    lam.glasso = lam.glasso[c(mo21, mo22, mo23, mo24)]
  )

  colnames(bgbm.optim.glasso$bgbb) <- names(bgbm.optim.glasso$lam.glasso) <- c("AIC", "BIC", "AIC.Orig", "BIC.Orig")
  bgbbest <- cbind(bgbm.init, bgbm.optim$par, bgbm.optim.glasso$bgbb)
  tree.fit <- treefit(Delta, Z[, -(1:length(TrainSurv_pred_org))]) # ??????


  set.seed(seed)
  train <- sample(1:nn, nn * 0.75)
  logi.fit <- logifit(Delta, Z, train, colnames(Z)[1:length(TrainSurv_pred_org)])
  #  vars <- logi.fit$vars
  #  vars <- sapply(vars, function(x) substr(x, 1, nchar(x) - 3), USE.NAMES = FALSE)
  #  vars <- unique(vars)
  vars_wk1 <- logi.fit$vars
  vars_wk2 <- substring(vars_wk1, regexpr("_", vars_wk1) + 1)
  vars_wk3 <- paste0("pred", vars_wk2)
  vars <- unique(vars_wk3)
  wei <- GetWei(TrainPK, vars, Delta, SX)

  # ================
  ### validation
  # ================
  SX <- ValidSurv$sx
  SC <- ValidSurv$sc
  Delta <- ValidSurv$delta

  ## the same PCA basis as training set (use PCAobjs from the training)
  tmp <- paste0("base_pred", 1:length(ValidSurv_pred_org))
  Z <- as.matrix(cbind(ValidSurv[, tmp], ValidFt))
  Z <- (Z - VTM(meanZ, nrow(Z))) / VTM(sdZ, nrow(Z))
  ValidFt_PCA <- do.call(cbind, sapply(1:length(codes), function(i) {
    if (length(codesK[[i]]) > 0) {
      tmp <- matrix(Z[, codesK[[i]] + length(ValidSurv_pred_org)], ncol = length(codesK[[i]]))
      colnames(tmp) <- colnames(Z)[codesK[[i]] + length(ValidSurv_pred_org)]
      tmp3 <- predict(PCAobjs[[i]], tmp)[, 1:PCA_K[i]]
      list(PCAs = tmp3)
    } else {
      list(PCAs = NULL)
    }
  }))
  colnames(ValidFt_PCA) <- PCAnames
  Z <- cbind(Z[, 1:length(ValidSurv_pred_org)], ValidFt_PCA)
  endF <- quantile(SX, 0.9)
  t.grid <- tseq[seq(9, length(tseq), 10)]

  #---- this uses the validation data to get G_fit
  G_fit <- survfit(Surv(SX, 1 - Delta) ~ 1)
  G_tseq <- sapply(tseq, function(s) {
    tmp2 <- G_fit$time <= s
    if (sum(tmp2) == 0) {
      1
    } else {
      G_fit$surv[max(which(tmp2))]
    }
  })
  G_SX <- summary(G_fit, times = SX, extend = TRUE)$surv
  G_SX[sort(SX, index.return = TRUE)$ix] <- G_SX
  G_SX[G_SX == 0] <- min(G_SX[G_SX > 0])
  G_tseq[G_tseq == 0] <- min(G_tseq[G_tseq > 0])

  eval <- GetPrecAll(
    bgbbest, SX, SX, SC, Delta, Z, tseq, t.grid,
    knots, Boundary.knots, G_SX, G_tseq, endF, Tend
  )
  tmp <- eval$CB
  tmp2 <- mean(BrierScore.KM2(tseq[tseq <= endF], SX, SX, SC, Delta, tseq, G_SX, G_tseq)[, 2])
  tmp[2, ] <- 1 - tmp[2, ] / tmp2
  cstats <- tmp

  # ================
  # Output
  # ================
  out <- list()
  #--- knots used for B-spline
  out$nparm4Bspline <- nparm4Bspline
  #--- regression parameter estimate (from Training) --
  out$bgbbest <- bgbbest
  out$gamma_est <- bgbbest[1:nparm4Bspline]
  out$beta_est <- bgbbest[-c(1:nparm4Bspline)]
  #--- PCA  (from Training) --
  out$PCAobjs <- PCAobjs
  #--- G_SX and G_tseq  (from Training) --
  out$GX_trian <- G_SX_train
  out$Gt_train <- G_tseq_train
  #--- group used for group lasso
  out$group <- group

  #--- result with the validation data 
  out$result_valid <- cstats
  #---- predicted score (validation data)
  out$risk_score_valid <- eval$risk_score
  # out$pred_surv_valid = eval$pred_surv

  #-- add class to the return object --
  class(out) <- "masta"

  return(out)
}
