NPMLE <- function(bb, HZ, SX, Z, Delta) {
  nn <- length(Delta)

  SX2 <- sort(unique(SX[Delta == 1]))
  indx <- match(SX, SX2)
  HZ <- exp(HZ)
  HZN <- HZ[indx]
  HZN[is.na(HZN)] <- 1
  HZsum <- cumsum(HZ)
  HZsum <- c(0, HZsum)
  tmp <- outer(SX, SX2, FUN = ">=")
  indx <- apply(tmp, 1, sum)
  HZsumN <- HZsum[indx + 1]
  Zbb <- as.numeric(Z %*% bb)

  l <- NULL
  l$likelihood <- mean(Zbb - Delta * log(HZN) + (1 + Delta) * log(HZsumN + exp(-Zbb)))
  lbb <- apply(Z - matrix((1 + Delta) * exp(-Zbb) / {
    HZsumN + exp(-Zbb)
  }, nrow = nn, ncol = length(bb)) * Z, 2, mean)
  lHZ <- 1 / nn - apply(tmp * VTM(HZ, nn) * (1 + Delta) / {
    HZsumN + exp(-Zbb)
  }, 2, mean)
  l$gradient <- c(lbb, lHZ)

  return(l)
}
NULL

NPMLE.est <- function(bb.init, HZ.init, SX, Z, Delta) {
  lHZ <- length(unique(SX[Delta == 1]))
  SX.grid <- sort(unique(SX[Delta == 1]))
  parinit <- c(bb.init, HZ.init)

  lbb <- length(bb.init)
  aa <- optim(
    par = parinit, fn = function(x) NPMLE(x[1:lbb], x[-c(1:lbb)], SX, Z, Delta)$likelihood,
    gr = function(x) NPMLE(x[1:lbb], x[-c(1:lbb)], SX, Z, Delta)$gradient,
    method = "CG", control = list(maxit = 10000)
  )
  bb <- aa$par[1:lbb]
  HZ <- log(cumsum(exp(aa$par[-c(1:lbb)])))
  HZ <- cbind(SX.grid, HZ)

  list(bb = bb, bbt = HZ, convergence = aa$convergence)
}
NULL

BrierScore.NPMLE2 <- function(tt = NULL, bb, bbt, ST, SX, SC, Delta, Z, tseq, GX = NULL, Gt = NULL) {
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

  tmp <- stepfun(bbt[-1, 1], bbt[, 2])
  tmp <- tmp(tseq)
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
  return(cbind(tt, aa))
}
NULL

NPMLE.Pred <- function(bb, bbt, ST, SX, SC, Delta, Z, tseq, t.grid, GX, Gt, endF) {
  Cstat.CV <- Est.Cval(cbind(SX, Delta, Z %*% bb), tau = endF, nofit = TRUE)$Dhat

  ## new definition
  BrierS.CV <- BrierScore.NPMLE2(
    tseq[tseq <= endF], bb, bbt,
    ST, SX, SC, Delta, Z, tseq, GX, Gt
  )[, 2]
  BrierS.CV <- mean(BrierS.CV)

  CB <- c(Cstat.CV, BrierS.CV)
  names(CB) <- c("Cstat", "BrierSc.Adj")

  return(CB)
}
NULL
