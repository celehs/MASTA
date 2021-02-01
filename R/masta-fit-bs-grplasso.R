### All functions need in simulations
### Survival likelihood; Ridge, Glasso; Prediction Accuracy; Rule based algorithm (comparison)

### No penalty on BS at all in glasso and ridge
### do not optimize over BS at all in glasso and ridge

### change hessian to exact (not numerical derivative)

### Jul 31 Add SurvPred & SurvEst

### Aug 3 delete dstep in ridge AIC, all Norm simus code might need to be changed ??

### Aug 10 midnight, add the case that tree has 0 split

###########################################################################
### reparametrized intercept function h(t) = exp(beta_0(t))
### the intercept beta_0 is reparametrized to ensure monotone increasing
### BS intercept is constant (in bs function intercept=FALSE, add then add 1)
h.fun <- function(t, knots, Boundary.knots, bg, a = NULL, b = NULL, k = NULL) {
  allknots <- c(Boundary.knots[1], knots, Boundary.knots[2])
  q <- length(allknots)

  if (is.null(a)) {
    a <- {
      allknots[-1] * bg[-q] - allknots[-q] * bg[-1]
    } / (allknots[-1] - allknots[-q])
  }
  if (is.null(b)) {
    b <- {
      bg[-1] - bg[-q]
    } / (allknots[-1] - allknots[-q])
  }
  if (is.null(k)) k <- sapply(t, function(s) k <- sum(s > allknots))
  c <- ifelse(bg[-q] == bg[-1],
    exp(bg[-q]) * {
      allknots[-1] - allknots[-q]
    },
    {
      allknots[-1] - allknots[-q]
    } /
      {
        bg[-1] - bg[-q]
      } *
      {
        exp(a[-q] + b[-q] * allknots[-1]) - exp(bg[-q])
      }
  )
  c <- c * exp(bg[1])
  c[1] <- {
    allknots[2] - allknots[1]
  } / bg[2] * exp(bg[1]) * {
    exp(bg[2]) - 1
  }

  sapply(seq_along(t), function(i) {
    if (k[i] == 1) {
      tmp <- {
        allknots[2] - allknots[1]
      } / bg[2] * exp(bg[1]) * {
        exp(bg[2] * (t[i] - allknots[1]) / (allknots[2] - allknots[1])) - 1
      }
    } else {
      tmp <- exp(bg[1]) * ifelse(bg[k[i] + 1] == bg[k[i]], exp(bg[k[i]]) * (t[i] - allknots[k[i]]), {
        allknots[k[i] + 1] - allknots[k[i]]
      } /
        {
          bg[k[i] + 1] - bg[k[i]]
        } *
        {
          exp(a[k[i]] + b[k[i]] * t[i]) - exp(bg[k[i]])
        })
    }
    sum(c[1:k[i] - 1]) + tmp
  })
}


# derivative
# time: matrix calculation < numerical derivative < sapply
dh.fun <- function(t, knots, Boundary.knots, bg, a = NULL, b = NULL, k = NULL) {
  allknots <- c(Boundary.knots[1], knots, Boundary.knots[2])
  q <- length(allknots)

  if (is.null(a)) {
    a <- {
      allknots[-1] * bg[-q] - allknots[-q] * bg[-1]
    } / (allknots[-1] - allknots[-q])
  }
  if (is.null(b)) {
    b <- {
      bg[-1] - bg[-q]
    } / (allknots[-1] - allknots[-q])
  }
  if (is.null(k)) k <- sapply(t, function(s) k <- sum(s > allknots))

  dh <- matrix(0, nrow = length(t), ncol = length(bg))
  dh[, 1] <- h.fun(t, knots, Boundary.knots, bg, a = NULL, b = NULL, k = NULL)
  dh[k == 1, 2] <- exp(bg[1]) * {
    {
      bg[2] * (t[k == 1] - allknots[1]) - (allknots[2] - allknots[1])
    } / bg[2]^2 * exp(bg[2] *
      {
        t[k == 1] - allknots[1]
      } / (allknots[2] - allknots[1])) +
      (allknots[2] - allknots[1]) / bg[2]^2
  }
  dh[k == 2, 2] <- exp(bg[1]) * {
    (bg[2] - 1) * (allknots[2] - allknots[1]) / bg[2]^2 * exp(bg[2]) +
      (allknots[2] - allknots[1]) / bg[2]^2 +
      {
        (allknots[3] - allknots[2]) / (bg[3] - bg[2])^2 + (allknots[3] - t[k == 2]) / (bg[3] - bg[2])
      } * exp(a[2] + b[2] * t[k == 2]) -
      {
        (allknots[3] - allknots[2]) * (bg[3] - bg[2] + 1) / (bg[3] - bg[2])^2
      } * exp(bg[2])
  }
  dh[k > 2, 2] <- exp(bg[1]) * {
    (bg[2] - 1) * (allknots[2] - allknots[1]) / bg[2]^2 * exp(bg[2]) +
      (allknots[2] - allknots[1]) / bg[2]^2 +
      {
        (allknots[3] - allknots[2]) / (bg[3] - bg[2])^2
      } * exp(bg[3]) -
      {
        (allknots[3] - allknots[2]) * (bg[3] - bg[2] + 1) / (bg[3] - bg[2])^2
      } * exp(bg[2])
  }
  dh[cbind(which(k > 2), k[k > 2])] <- exp(bg[1]) * {
    {
      (allknots[k[k > 2]] - allknots[k[k > 2] - 1]) * (bg[k[k > 2]] - bg[k[k > 2] - 1] - 1) / (bg[k[k > 2]] - bg[k[k > 2] - 1])^2 -
        (allknots[k[k > 2] + 1] - allknots[k[k > 2]]) * (bg[k[k > 2] + 1] - bg[k[k > 2]] + 1) / (bg[k[k > 2] + 1] - bg[k[k > 2]])^2
    } * exp(bg[k[k > 2]]) +
      (allknots[k[k > 2]] - allknots[k[k > 2] - 1]) / (bg[k[k > 2]] - bg[k[k > 2] - 1])^2 * exp(bg[k[k > 2] - 1]) + {
        (allknots[k[k > 2] + 1] - allknots[k[k > 2]]) / (bg[k[k > 2] + 1] - bg[k[k > 2]])^2 + (allknots[k[k > 2] + 1] - t[k > 2]) / (bg[k[k > 2] + 1] - bg[k[k > 2]])
      } * exp(a[k[k > 2]] + b[k[k > 2]] * t[k > 2])
  }
  dh[cbind(which(k > 1), k[k > 1] + 1)] <- exp(bg[1]) * {
    {
      (t[k > 1] - allknots[k[k > 1]]) / (bg[k[k > 1] + 1] - bg[k[k > 1]]) - (allknots[k[k > 1] + 1] - allknots[k[k > 1]]) / (bg[k[k > 1] + 1] - bg[k[k > 1]])^2
    } * exp(a[k[k > 1]] + b[k[k > 1]] * t[k > 1]) +
      (allknots[k[k > 1] + 1] - allknots[k[k > 1]]) / (bg[k[k > 1] + 1] - bg[k[k > 1]])^2 * exp(bg[k[k > 1]])
  }
  for (i in 3:(q - 2)) {
    dh[k > i, i] <- exp(bg[1]) * {
      (allknots[i] - allknots[i - 1]) * (bg[i] - bg[i - 1] - 1) / (bg[i] - bg[i - 1])^2 * exp(bg[i]) +
        (allknots[i] - allknots[i - 1]) * exp(bg[i - 1]) / (bg[i] - bg[i - 1])^2 +
        (allknots[i + 1] - allknots[i]) * exp(bg[i + 1]) / (bg[i + 1] - bg[i])^2 -
        (allknots[i + 1] - allknots[i]) * (bg[i + 1] - bg[i] + 1) / (bg[i + 1] - bg[i])^2 * exp(bg[i])
    }
  }
  dh
}


## get initial value for bg
# bgbm<-function(tseq, beta0c, Boundary.knots){
#   bb0t   = 3*log(tseq)+beta0c
#   ht     = (tseq)^3*exp(beta0c)
#   dht    = log(3*tseq^2)+beta0c
#   B      = bs(tseq,degree = 1,knots = knots, Boundary.knots = Boundary.knots, intercept = FALSE)
#   B      = cbind(1,B)
#   bg_bm  = solve(t(B)%*%B)%*%t(B)%*%dht  ## based on linear regression
#   as.numeric(bg_bm)
# }


###########################################################################
### Survival loglikelihood and score; + ridge & group lasso
## bb = \bb_1
## -n^{-1} loglik
SurvLoglik <- function(bg, bb, knots, Boundary.knots, Delta, SX, Z,
                       likelihood = TRUE, gradient = FALSE) {
  allknots <- c(Boundary.knots[1], knots, Boundary.knots[2])
  q <- length(allknots)

  BX <- bs(SX, knots = knots, degree = 1, Boundary.knots = Boundary.knots, intercept = FALSE)
  BX <- cbind(1, BX)

  a <- {
    allknots[-1] * bg[-q] - allknots[-q] * bg[-1]
  } / (allknots[-1] - allknots[-q])
  b <- {
    bg[-1] - bg[-q]
  } / (allknots[-1] - allknots[-q])
  kk <- sapply(SX, function(s) k <- sum(s > allknots))

  hX <- h.fun(SX, knots, Boundary.knots, bg, a, b, kk)

  l <- NULL

  if (likelihood) {
    l$likelihood <- -mean(Delta * as.numeric({
      BX %*% bg + Z %*% bb
    }) - (1 + Delta) * log(1 + hX * exp(as.numeric(Z %*% bb))))
  }

  if (gradient) {
    dhX <- dh.fun(SX, knots, Boundary.knots, bg, a, b, kk)
    tmp <- exp(Z %*% bb) / (1 + hX * exp(Z %*% bb))
    lbg <- apply(t(VTM(Delta, ncol(dhX))) * BX - t(VTM(
      {
        1 + Delta
      } * tmp,
      ncol(dhX)
    )) * dhX, 2, mean)
    lbb <- apply(t(VTM(Delta - (1 + Delta) * hX * tmp, ncol(Z))) * Z, 2, mean)
    l$gradient <- -c(lbg, lbb)
  }

  return(l)
}


## nzpar is over bb
SurvLoglik.nz <- function(bg, bb, nzpar, knots, Boundary.knots, Delta, SX, Z,
                          likelihood = TRUE, gradient = FALSE) {
  SurvLoglik(
    bg = bg, bb = bb, knots = knots, Boundary.knots = Boundary.knots,
    Delta = Delta, SX = SX, Z = Z[, nzpar], likelihood = likelihood, gradient = gradient
  )
}

gglasso.Approx.BSNP <- function(bg, bb, knots, Boundary.knots, Delta, SX, Z, lam.rid, lam.glasso, group) {
  nn <- length(SX)

  hX <- h.fun(SX, knots, Boundary.knots, bg)
  tmp <- as.numeric(hX * exp(Z %*% bb))
  tmp <- tmp / (1 + tmp)^2
  d2l <- t(Z) %*% diag((1 + Delta) * tmp / nn) %*% Z + lam.rid * diag(rep(1, length(bb)))

  d2l <- svd(d2l)
  d2l.sqrt <- d2l$u %*% diag(sqrt(d2l$d)) %*% t(d2l$v)
  ynew <- d2l.sqrt %*% bb

  pf_glasso <- 1 / sqrt(aggregate(bb^2, by = list(group), FUN = sum)[, 2])
  pf_glasso[is.infinite(pf_glasso)] <- max(pf_glasso[!is.infinite(pf_glasso)])
  tmp <- gglasso(d2l.sqrt, ynew,
    group = group, loss = "ls",
    lambda = lam.glasso, pf = pf_glasso,
    intercept = FALSE
  )

  tmp2 <- apply((c(ynew) - d2l.sqrt %*% tmp$beta)^2, 2, sum)
  AIC.LSA <- tmp2 + tmp$df * 2 / nn
  BIC.LSA <- tmp2 + tmp$df * log(nn) / nn

  tmp2 <- unlist(sapply(seq_along(lam.glasso), function(i) {
    SurvLoglik(
      bg, tmp$beta[, i], knots, Boundary.knots,
      Delta, SX, Z
    )
  }))
  AIC.Orig <- tmp2 + tmp$df * 2 / nn
  BIC.Orig <- tmp2 + tmp$df * log(nn) / nn

  tmp2 <- list(
    AIC.LSA = AIC.LSA, BIC.LSA = BIC.LSA, AIC.Orig = AIC.Orig, BIC.Orig = BIC.Orig,
    LL = apply((c(ynew) - d2l.sqrt %*% tmp$beta)^2, 2, sum),
    Loglik = unlist(sapply(seq_along(lam.glasso), function(i) {
      SurvLoglik(
        bg, tmp$beta[, i], knots, Boundary.knots,
        Delta, SX, Z
      )
    }))
  )

  return(c(tmp, tmp2))
}


###
OrigSc <- function(bgbbest, mean, sd, bgbbbm = NULL) {
  bgbbest2 <- bgbbest[, -1]
  bgbbest2[1, ] <- bgbbest2[1, ] - apply(bgbbest2[-c(1:q), ] * mean / sd, 2, sum)
  bgbbest2[-c(1:q), ] <- bgbbest2[-c(1:q), ] / sd

  if (!is.null(bgbbbm)) {
    bgbbest2 <- cbind(bgbbbm, bgbbest2)
    colnames(bgbbest2)[1] <- "BM"
  }

  bgbbest2
}



### Get MLE, Ridge, GLASSO results
GetResultAll <- function(bgbm.init, knots, Boundary.knots, q,
                         Delta, SX, Z,
                         lam.rid, lam.glasso, group) {
  ## optim with initial bg_bm
  bgbm.optim <- optim(
    par = bgbm.init, fn = function(x) SurvLoglik(bg = x[1:q], bb = x[-c(1:q)], knots, Boundary.knots, Delta, SX, Z)$likelihood,
    gr = function(x) SurvLoglik(bg = x[1:q], bb = x[-c(1:q)], knots, Boundary.knots, Delta, SX, Z, likelihood = FALSE, gradient = TRUE)$gradient,
    method = "BFGS", control = list(maxit = 3000)
  )[c("par", "convergence")]

  ## group lasso with L2 approximation; no penalty on BS part
  bgbm.optim.glasso <- gglasso.Approx.BSNP(bgbm.optim$par[1:q], bgbm.optim$par[-c(1:q)], knots, Boundary.knots, Delta, SX, Z, 0, lam.glasso, group)
  # AIC
  mo21 <- which.min(bgbm.optim.glasso$AIC.LSA)
  # BIC
  mo22 <- which.min(bgbm.optim.glasso$BIC.LSA)
  # AIC on original model
  mo23 <- which.min(bgbm.optim.glasso$AIC.Orig)
  # BIC on original model
  mo24 <- which.min(bgbm.optim.glasso$BIC.Orig)
  # non-zero parameters
  nzpar <- cbind(bgbm.optim.glasso$beta[, c(mo21, mo22, mo23, mo24)]) != 0
  nzpar <- rbind(matrix(TRUE, nrow = q, ncol = 4), nzpar)
  bgbm.optim.glasso <- sapply(1:4, function(i) {
    tmp <- rep(0, q + ncol(Z))
    tmp[nzpar[, i]] <- optim(
      par = bgbm.optim$par[nzpar[, i]], fn = function(x) SurvLoglik.nz(bg = x[1:q], bb = x[-c(1:q)], nzpar[-c(1:q), i], knots, Boundary.knots, Delta, SX, Z)$likelihood,
      gr = function(x) SurvLoglik.nz(bg = x[1:q], bb = x[-c(1:q)], nzpar[-c(1:q), i], knots, Boundary.knots, Delta, SX, Z, likelihood = FALSE, gradient = TRUE)$gradient[nzpar[, i]],
      method = "BFGS", control = list(maxit = 3000)
    )$par
    tmp
  })

  colnames(bgbm.optim.glasso) <- c("AIC", "BIC", "AIC.Orig", "BIC.Orig")

  bgbbest <- cbind(bgbm.init, bgbm.optim$par, bgbm.optim.glasso)
  colnames(bgbbest) <- c(
    "BM", "MLE",
    paste0("Glasso.", c("AIC", "BIC", "AIC.Orig", "BIC.Orig"))
  )

  bgbbest
}

###########################################################################


###########################################################################
### Some basic functions
## logistic model
# distribution function of T given features Z
g_fun <- function(x) 1 / (1 + exp(-x))

###########################################################################


###########################################################################
### prediction accuracy: Brier Score
BrierScore2 <- function(tt = NULL, bg, bb, ST, SX, SC, Delta, Z, knots, Boundary.knots, tseq, GX = NULL, Gt = NULL) {
  # x=bgbbest[,3]
  #  tt=tseq[tseq<=endF]; bg=x[1:q]; bb=x[-c(1:q)] ; ST=SX;
  # == ST SX SC Delta Z knots Boundary.knots tseq ;
  #  GX=G_SX; Gt=G_tseq ;

  if (is.null(GX) & is.null(Gt)) {
    G_fit <- survfit(Surv(SX, 1 - Delta) ~ 1, type = "kaplan-meier")
    GX <- sapply(SX, function(s) {
      tmp2 <- G_fit$time <= s
      if (sum(tmp2) == 0) {
        1
      } else {
        G_fit$surv[max(which(tmp2))]
      }
    })
    Gt <- sapply(tseq, function(s) {
      tmp2 <- G_fit$time <= s
      if (sum(tmp2) == 0) {
        1
      } else {
        G_fit$surv[max(which(tmp2))]
      }
    })
  }

  if (is.null(tt)) tt <- tseq[Gt != 0]
  tmp <- sapply(tt, function(t) {
    {
      SX <= t
    } * Delta / GX
  })
  tmp[is.na(tmp)] <- 0
  wt <- tmp + sapply(tt, function(t) {
    {
      SX > t
    } / Gt[tseq == t]
  })

  tmp <- h.fun(tseq, knots, Boundary.knots, bg)
  tmp <- log(tmp)
  # tmp   = log(cumsum(exp(B%*%bg)))+log(diff(tseq)[1])

  #--- peak time point?
  tmpp <- apply(outer(SC, tseq, FUN = function(x, y) abs(x - y)), 1, which.min)
  tmpp <- tmp[tmpp]
  tmpp <- outer(tmpp, tmp, FUN = "pmin")
  tmp2 <- outer(ST, tseq, FUN = "<=") * Delta # 1(T<= min(t,C))
  aa <- apply(
    {
      tmp2[, tseq %in% tt] - g_fun(as.numeric(Z %*% bb) + tmpp[, tseq %in% tt])
    }^2 * wt,
    2,
    mean
  )
  #--- returning ave[I(T>t) - S(t)]^2 for all t.
  #--- predicted survival function: 
  pred_surv <- g_fun(as.numeric(Z %*% bb) + tmpp[, tseq %in% tt])

  # return(cbind(tt,aa))
  out <- list()
  out$BrierSore_tt <- cbind(tt, aa)
  out$pred_surv <- pred_surv
  return(out)
}



GetPrecAll <- function(bgbbest, ST, SX, SC, Delta, Z, tseq, t.grid,
                       knots, Boundary.knots, GX, Gt, endF, Tend) {
  # bgbbest,SX,SX,SC,Delta,Z,tseq,t.grid,knots,Boundary.knots,G_SX,G_tseq,endF,Tend
  # ST=SX; GX=G_SX; Gt=G_tseq ;

  q <- length(knots) + 2

  #--- this does not use the b-spline part. Using only predictors (does not use intercept) 
  Cstat.CV <- apply(bgbbest, 2, function(x) {
    Est.Cval(cbind(SX, Delta, Z %*% x[-c(1:q)]), tau = endF, nofit = TRUE)$Dhat
  })

  #--- risk score -- 
  risk_score <- apply(bgbbest, 2, function(x) {
    Z %*% x[-c(1:q)]
  })

  ### Brier Score ###
  #-- for each parameter, get the brier score from 0 to endF ;
  BrierS.CV <- apply(bgbbest, 2, function(x) {
    BrierScore2(
      tseq[tseq <= endF], x[1:q], x[-c(1:q)],
      ST, SX, SC, Delta, Z, knots, Boundary.knots, tseq, GX, Gt
    )$BrierSore_tt[, 2]
  })
  #-- average over the time from 0 to endF ??
  BrierS.CV <- apply(BrierS.CV, 2, mean)

  CB <- rbind(Cstat.CV, BrierS.CV)
  row.names(CB) <- c("Cstat", "BrierSc.Adj")

  #--- retrun predicted survival function ---
  # chk=BrierScore2(tseq[tseq<=endF],bgbbest[1:q,3],bgbbest[-c(1:q),3],ST,SX,SC,Delta,Z,knots,Boundary.knots,tseq,GX,Gt)$pred_surv
  pred_surv <- list()
  for (kk in 1:ncol(bgbbest)) {
    tmp <- BrierScore2(
      tseq[tseq <= endF], bgbbest[1:q, kk], bgbbest[-c(1:q), kk],
      ST, SX, SC, Delta, Z, knots, Boundary.knots, tseq, GX, Gt
    )
    pred_surv[[kk]] <- tmp$pred_surv
    rownames(pred_surv[[kk]]) <- paste0("case", 1:nrow(Z))
    colnames(pred_surv[[kk]]) <- tmp$BrierSore_tt[, 1]
  }


  # return(CB)
  out <- list()
  out$CB <- CB
  out$risk_score <- risk_score
  out$pred_surv <- pred_surv
  return(out)
}


### Decision tree
BrierScore.tree2 <- function(tt = NULL, ST, SX, Delta, Z,
                             SX.tree, Delta.tree, tseq, GX, Gt) {
  if (is.null(GX) & is.null(Gt)) {
    G_fit <- survfit(Surv(SX, 1 - Delta) ~ 1, type = "kaplan-meier")
    GX <- sapply(SX, function(s) {
      tmp2 <- G_fit$time <= s
      if (sum(tmp2) == 0) {
        1
      } else {
        G_fit$surv[max(which(tmp2))]
      }
    })
    Gt <- sapply(tseq, function(s) {
      tmp2 <- G_fit$time <= s
      if (sum(tmp2) == 0) {
        1
      } else {
        G_fit$surv[max(which(tmp2))]
      }
    })
  }
  if (is.null(tt)) tt <- tseq[Gt != 0]
  tmp <- sapply(tt, function(t) {
    {
      SX <= t
    } * Delta / GX
  })
  tmp[is.na(tmp)] <- 0
  wt <- tmp + sapply(tt, function(t) {
    {
      SX > t
    } / Gt[tseq == t]
  })

  tmp <- outer(SX.tree, tseq, FUN = "<=") * Delta.tree
  tmp2 <- outer(ST, tseq, FUN = "<=") * Delta
  aa <- apply(
    {
      tmp2[, tseq %in% tt] - tmp[, tseq %in% tt]
    }^2 * wt,
    2,
    mean
  )
  return(cbind(tt, aa))
}


GetPrecAll.tree <- function(tree.fit, ST, SX, SC, Delta, Z, FirstCode, tseq, t.grid, GX, Gt, endF) {
  ## Terminal nodes (TN)
  treepred.fit <- treepred(tree.fit, SC, Z, FirstCode, type = "terminal nodes")
  ## All vars included in the tree model (All)
  treepred.fit2 <- treepred(tree.fit, SC, Z, FirstCode, type = "all vars")
  ## Cstat
  Cstat.CV <- c(
    TreeTN = Est.Cval(cbind(SX, Delta, treepred.fit$SX.tree), tau = endF)$Dhat,
    TreeAll = Est.Cval(cbind(SX, Delta, treepred.fit2$SX.tree), tau = endF)$Dhat
  )
  ## Brier
  BrierS.CV <- mean(BrierScore.tree2(
    tseq[tseq <= endF], ST, SX, Delta, Z,
    treepred.fit$SX.tree, treepred.fit$Delta.tree, tseq, GX, Gt
  )[, 2])
  BrierS.CV <- c(
    BrierS.CV,
    mean(BrierScore.tree2(
      tseq[tseq <= endF], ST, SX, Delta, Z,
      treepred.fit2$SX.tree, treepred.fit2$Delta.tree, tseq, GX, Gt
    )[, 2])
  )
  names(BrierS.CV) <- c("TreeTN", "TreeAll")

  CB <- rbind(Cstat.CV, BrierS.CV)
  row.names(CB) <- c("Cstat", "BrierSc.Adj")

  return(CB)
}

###########################################################################



###########################################################################
##### decision tree (rule based algorithm)
treefit <- function(Delta, Z) {
  df <- data.frame(Delta = Delta, Z)
  fit <- rpart(Delta ~ ., data = df, method = "class")
  opt <- which.min(fit$cptable[, "xerror"])
  fit <- prune(fit, cp = fit$cptable[opt, "CP"])
  return(fit)
}


treepred <- function(fit, SC, Z, FirstCode, type = c("terminal nodes", "all vars")) {
  type <- match.arg(type)

  Delta.tree <- as.numeric(predict(fit, data.frame(Z), type = "class")) - 1
  SX.tree <- SC

  if (nrow(fit$frame) > 1) {
    tmp <- rpart.subrules.table(fit)
    if (type == "all vars") {
      vars <- unique(as.character(tmp$Variable))
      vars <- as.numeric(sapply(vars, function(x) substr(x, nchar(x), nchar(x))))
      vars[vars == 0] <- as.numeric(sapply(unique(as.character(tmp$Variable))[vars == 0], function(x) substr(x, nchar(x) - 1, nchar(x))))
      vars <- unique(vars)
      if (length(vars) == 1) {
        SX.tree[Delta.tree == 1] <- FirstCode[Delta.tree == 1, vars]
      } else {
        SX.tree[Delta.tree == 1] <- apply(FirstCode[Delta.tree == 1, vars], 1, min, na.rm = TRUE)
      }
    } else if (type == "terminal nodes") {
      # tmp2 = as.integer(row.names(fit$frame))[rpart:::pred.rpart(fit,Z)]
      tmp2 <- as.integer(row.names(fit$frame))[predict(fit, data.frame(Z))] # added by MY
      tmp3 <- rpart.rules.table(fit)
      tmp3$Var <- tmp$Variable[match(tmp3$Subrule, tmp$Subrule)]
      tmp3$Level <- as.numeric(sapply(as.character(tmp3$Subrule), function(x) {
        ifelse(x == "NULL", 0, substr(x, 2, 2))
      }))
      tmp3 <- data.frame(num = unique(tmp3$Rule), var = sapply(unique(tmp3$Rule), function(x) {
        with(tmp3[tmp3$Rule %in% x, ], Var[which.max(Level)])
      }))
      vars <- as.character(tmp3$var[match(tmp2, tmp3$num)])
      vars <- as.numeric(sapply(vars, function(x) substr(x, nchar(x), nchar(x))))
      vars[vars == 0] <- as.numeric(sapply(as.character(tmp3$var[match(tmp2, tmp3$num)])[vars == 0], function(x) substr(x, nchar(x) - 1, nchar(x))))
      # 10
      indx <- Delta.tree == 1
      indx <- cbind(which(indx), vars[indx])
      SX.tree[Delta.tree == 1] <- FirstCode[indx]
      SX.tree[is.na(SX.tree)] <- SC[is.na(SX.tree)] / 2
    }
  }

  return(data.frame(Delta.tree, SX.tree))
}

###########################################################################


###########################################################################
### Logistic regression, two-step
logifit <- function(Delta, Z, train, baseline) {
  df <- data.frame(Delta, Z)
  tmp <- glm(Delta ~ ., data = df[train, ], family = binomial)
  fit <- cv.glmnet(Z[train, ], Delta[train],
    family = "binomial", alpha = 1,
    penalty.factor = abs(1 / tmp$coefficient), standardize = FALSE
  )
  cc <- as.numeric(coef(fit, s = "lambda.min"))
  fm <- as.formula(paste0("Delta~", paste0(colnames(Z)[cc[-1] != 0], collapse = "+")))
  fit <- glm(fm, data = df[train, ], family = binomial)

  fit$vars <- colnames(Z)[cc[-1] != 0 & {
    !colnames(Z) %in% baseline
  }]

  pp <- predict(fit, newdata = df[-train, ], type = "response")

  fit$thresh <- as.numeric(optimize(f = function(c) {
    sum(abs(Delta[-train] - {
      pp >= c
    }))
  }, interval = c(0, 1))$minimum)

  fit$Delta <- predict(fit, newdata = df, type = "response") > fit$thresh

  return(fit)
}


GetWei <- function(PK, vars, Delta, SX) {
  tmp <- Delta == 1
  lv <- length(vars)

  if (lv == 1) {
    wei <- 1
    adj <- median(PK[tmp, vars] - SX[tmp], na.rm = TRUE)
  } else {
    wei <- 1 / apply(PK[tmp, vars], 2, function(x) {
      rg <- quantile(x, prob = c(0.03, 0.97), na.rm = TRUE)
      var(x[x <= rg[2] & x >= rg[1]], na.rm = TRUE)
    })
    wei <- wei / sum(wei)
    adj <- apply(PK[tmp, vars] - SX[tmp], 2, median, na.rm = TRUE)
  }

  list(wei = wei, adj = adj) # add adj
}


GetSX <- function(PK, wei, adj, logi.fit, vars, Z, SC, Tend) { ## add adj
  df <- data.frame(Z)
  nn <- nrow(Z)
  Delta.fit <- predict(logi.fit, newdata = df, type = "response") >= logi.fit$thresh

  PK[, vars] <- PK[, vars] - VTM(adj, nn)

  PK[is.na(PK)] <- 0
  SX.fit <- SC
  # if(any(Delta.fit==1)){
  #
  # }
  if (length(logi.fit$vars) == 1) {
    SX.fit[Delta.fit] <- PK[Delta.fit, vars]
  } else {
    SX.fit[Delta.fit] <- {
      PK[Delta.fit, vars] %*% wei
    } / {
      matrix(PK[Delta.fit, vars] != 0, ncol = length(vars)) %*% wei
    }
  }

  SX.fit[is.na(SX.fit)] <- Tend

  return(data.frame(Delta = as.numeric(Delta.fit), SX = SX.fit))
}


GetPrecAll.logi <- function(logi.fit, vars, wei, adj, PKTS, ST, SX, SC, Delta, Z,
                            tseq, t.grid, GX, Gt, endF, Tend) {
  tmp <- GetSX(PKTS, wei, adj, logi.fit, vars, Z, SC, Tend)
  ## Cstat
  Cstat.CV <- Est.Cval(cbind(SX, Delta, tmp$SX), tau = endF)$Dhat
  ## Brier
  BrierS.CV <- mean(BrierScore.tree2(
    tseq[tseq <= endF], ST, SX, Delta, Z,
    tmp$SX, tmp$Delta, tseq, GX, Gt
  )[, 2])
  CB <- c(Cstat.CV, BrierS.CV)
  names(CB) <- c("Cstat", "BrierSc.Adj")

  return(CB)
}

###########################################################################


GetPrecAll.Comb <- function(bgbbest, bb_Cheng, bbt_Cheng,
                            bb, bbt, ST, SX, SC, Delta, Z, tseq, t.grid,
                            knots, Boundary.knots, GX, Gt, endF, Tend,
                            tree.fit, FirstCode,
                            logi.fit, vars, wei, adj, PKTS) {
  ### ours
  CB.CV <- GetPrecAll(
    bgbbest, ST, SX, SC, Delta, Z, tseq, t.grid,
    knots, Boundary.knots, GX, Gt, endF, Tend
  )

  ### Cheng
  tmp <- NPMLE.Pred(bb_Cheng, bbt_Cheng, ST, SX, SC, Delta, Z, tseq, t.grid, GX, Gt, endF)
  CB.CV <- cbind(CB.CV, NPMLE = tmp)

  ### NPMLE
  tmp <- NPMLE.Pred(bb, bbt, ST, SX, SC, Delta, Z, tseq, t.grid, GX, Gt, endF)
  CB.CV <- cbind(CB.CV, NPMLE = tmp)

  ### decision tree
  tmp <- GetPrecAll.tree(tree.fit, ST, SX, SC, Delta, Z, FirstCode, tseq, t.grid, GX, Gt, endF)
  CB.CV <- cbind(CB.CV, tmp)

  ### logistic reg
  tmp <- GetPrecAll.logi(logi.fit, vars, wei, adj, PKTS, ST, SX, SC, Delta, Z, tseq, t.grid, GX, Gt, endF, Tend)
  CB.CV <- cbind(CB.CV, logi = tmp)

  return(CB.CV)
}


BrierScore.KM2 <- function(tt = NULL, ST, SX, SC, Delta, tseq, GX = NULL, Gt = NULL) {
  if (is.null(GX) & is.null(Gt)) {
    G_fit <- survfit(Surv(SX, 1 - Delta) ~ 1, type = "kaplan-meier")
    Gt <- summary(G_fit, times = tseq, extend = TRUE)$surv
    GX <- summary(G_fit, times = SX, extend = TRUE)$surv
    GX[sort(SX, index.return = TRUE)$ix] <- GX
    GX[GX == 0] <- min(GX[GX > 0])
    Gt[Gt == 0] <- min(Gt[Gt > 0])
  }
  if (is.null(tt)) tt <- tseq[Gt != 0]
  tmp <- sapply(tt, function(t) {
    {
      SX <= t
    } * Delta / GX
  })
  tmp[is.na(tmp)] <- 0
  wt <- tmp + sapply(tt, function(t) {
    {
      SX > t
    } / Gt[tseq == t]
  })

  KM <- survfit(Surv(SX, Delta) ~ 1, type = "kaplan-meier")
  tmp <- 1 - summary(KM, times = tseq, extend = TRUE)$surv
  tmpp <- apply(outer(SC, tseq, FUN = function(x, y) abs(x - y)), 1, which.min)
  tmpp <- tmp[tmpp]
  tmpp <- outer(tmpp, tmp, FUN = "pmin")

  # ## exact result, too slow
  # tmpp  = t(sapply(SC,function(x) pmin(x,tseq)))
  # tmpp  = 1-t(apply(tmpp,1,function(x) summary(G_fit,times=x,extend=TRUE)$surv))

  tmp2 <- outer(ST, tseq, FUN = "<=") * Delta # 1(T<= min(t,C))
  aa <- apply(
    {
      tmp2[, tseq %in% tt] - tmpp[, tseq %in% tt]
    }^2 * wt,
    2,
    mean
  )
  return(cbind(tt, aa))
}
