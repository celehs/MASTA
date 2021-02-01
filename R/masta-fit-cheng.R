xi.fun <- function(s) {
  tmp <- {
    exp(s) * (s - 1) + 1
  } / (exp(s) - 1)^2
  tmp[s == 0] <- 0.5
  tmp[s > 300] <- 0
  tmp
}


xi.deriv.fun <- function(s) {
  tmp <- {
    2 *
      {
        exp(s) - 1
      } - s * {
        exp(s) + 1
      }
  } / (exp(s) - 1)^2 / (1 - exp(-s))
  tmp[s == 0] <- -1 / 6
  tmp[s > 300] <- 0
  tmp
}


estbb <- function(bb, Delta, SX, Z) {
  G_fit <- survfit(Surv(SX, 1 - Delta) ~ 1)
  G_SX <- summary(G_fit, times = SX, extend = TRUE)$surv
  G_SX[sort(SX, index.return = TRUE)$ix] <- G_SX
  G_SX[G_SX == 0] <- min(G_SX[G_SX > 0])
  tmp <- VTM(Delta / G_SX^2, length(SX)) * outer(SX, SX, FUN = ">=")
  Zb <- as.numeric(Z %*% bb)
  Zb.outer <- outer(Zb, Zb, FUN = "-")
  A <- (tmp - xi.fun(Zb.outer))
  A <- apply(A, 1, sum) - apply(A, 2, sum)
  est <- apply(Z * A, 2, sum)
  A <- -xi.deriv.fun(Zb.outer)
  Jacob <- -t(Z) %*%
    {
      A + t(A)
    } %*% Z + t(Z) %*% (apply(A, 1, sum) * Z) + t(Z) %*% (apply(A, 2, sum) * Z)
  return(list(est = est, Jacob = Jacob))
}


estbb.data <- function(bb, Z, Delta, G_SX, SX) {
  tmp <- VTM(Delta / G_SX^2, length(SX)) * outer(SX, SX, FUN = ">=") # used in estbb.data()


  Zb <- as.numeric(Z %*% bb)
  Zb.outer <- outer(Zb, Zb, FUN = "-")
  A <- (tmp - xi.fun(Zb.outer))
  A <- apply(A, 1, sum) - apply(A, 2, sum)
  est <- apply(Z * A, 2, sum)
  A <- -xi.deriv.fun(Zb.outer)
  Jacob <- -t(Z) %*%
    {
      A + t(A)
    } %*% Z + t(Z) %*% (apply(A, 1, sum) * Z) + t(Z) %*% (apply(A, 2, sum) * Z)
  return(list(est = est, Jacob = Jacob))
}


Cheng.est <- function(Delta, SX, Z) {
  bb <- multiroot(
    f = function(bb) estbb(bb, Delta, SX, Z)$est, start = runif(ncol(Z), -0.1, 0.1),
    jacfunc = function(bb) estbb(bb, Delta, SX, Z)$Jacob, maxiter = 3000
  )
  iters <- 1
  while (bb$estim.precis > 1 & iters < 10) {
    bb <- multiroot(
      f = function(bb) estbb(bb, Delta, SX, Z)$est, start = runif(ncol(Z), -0.1, 0.1),
      jacfunc = function(bb) estbb(bb, Delta, SX, Z)$Jacob, maxiter = 3000
    )
    iters <- iters + 1
  }
  bb <- bb$root
  if (iters >= 10) {
    cat("Cheng's method might not converge\n")
    bb <- optim(
      par = runif(ncol(Z), -0.1, 0.1), fn = function(bb) mean(estbb(bb, Delta, SX, Z)$est^2),
      gr = function(bb) 2 * apply(estbb(bb, Delta, SX, Z)$Jacob * estbb(bb, Delta, SX, Z)$est, 2, mean),
      method = "BFGS", control = list(maxit = 6000)
    )
    bb <- bb$par
  }
  G_fit <- survfit(Surv(SX, 1 - Delta) ~ 1)
  SX.grid <- sort(unique(SX[Delta == 1]))
  G_SX <- summary(G_fit, times = SX.grid, extend = TRUE)$surv

  ht <- sapply(seq_along(SX.grid), function(i) {
    est <- function(a) {
      mean({
        SX >= SX.grid[i]
      } / G_SX[i] + g_fun(a + Z %*% bb)) - 1
    }
    if (est(-30) > 0) -30 else uniroot(est, lower = -30, upper = 1e10)$root
  })

  ht <- sort(ht)
  ht <- cbind(SX.grid, ht)

  return(list(bb = bb, HZ = ht))
}


corZExtreme <- function(Z, thresh) {
  corZ <- cor(Z)
  ind <- which(abs(corZ) > thresh & abs(corZ) != 1)
  ind1 <- ind %% ncol(Z)
  ind1[ind1 == 0] <- ncol(Z)
  ind2 <- ceiling(ind / ncol(Z))
  tmpp <- ind1 <= ind2
  cbind(colnames(Z)[ind1[tmpp]], colnames(Z)[ind2[tmpp]], corZ[ind[tmpp]])
}
